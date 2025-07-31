from __future__ import annotations
import os, time, pickle, sqlite3, re, tempfile
import logging
from pathlib import Path
from datetime import datetime
from functools import lru_cache
from typing import List, Sequence

import pandas as pd
import numpy as np
import networkx as nx
from bioservices import KEGG
from mygene import MyGeneInfo
from sqlalchemy import create_engine, text
from tqdm.auto import tqdm

from proximity_util import build_matrix, DistMatrix, compute_network_distances_GPU

log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')


class HerbDatabase:
    def __init__(self, db_file: str):
        self.db_file = db_file
        self.db_url = f"sqlite:///{db_file}"
        self.engine = create_engine(self.db_url)

    @lru_cache(maxsize=None)
    def get_standard_name(self, alias: str) -> str | None:
        sql = text("SELECT herb_cn_name FROM Herb_Aliases WHERE alias = :v OR herb_cn_name = :v")
        with self.engine.connect() as conn:
            row = conn.execute(sql, {"v": alias}).fetchone()
        if not row:
            log.warning("⚠️ Alias not found: %s", alias)
        return row[0] if row else None

    @lru_cache(maxsize=None)
    def get_ensp_ids(self, herb_name: str) -> List[str]:
        std_name = self.get_standard_name(herb_name)
        if not std_name:
            return []
        sql = text("SELECT ENSP_ID FROM Herb_targets_filtered WHERE herb = :name")
        with self.engine.connect() as conn:
            ids = [str(r[0]) for r in conn.execute(sql, {"name": std_name})]
        return ids


class TargetResolver:
    def __init__(self):
        self.mg = MyGeneInfo()
        self.kegg = KEGG()
        self.kegg.http_header = {"user-agent": "Mozilla/5.0"}

    def resolve(self, target_type: str, identifier: str) -> List[str]:
        if target_type == "disease":
            return self._from_disease(identifier)
        elif target_type == "pathway":
            return self._from_pathway(identifier)
        elif target_type == "genes":
            return self._from_gene(identifier)
        else:
            raise ValueError("Unknown target type")

    def _from_disease(self, disease_id: str, cutoff: float = 0.30) -> List[str]:
        df = pd.read_parquet("./Data/Disease_data")
        genes = df[(df.diseaseId == disease_id) & (df.score >= cutoff)]["targetId"].dropna().unique().tolist()
        hits = self.mg.querymany(genes, scopes="ensembl.gene", fields="ensembl.protein", species="human")
        return self._extract_proteins(hits)

    def _from_pathway(self, path_id: str) -> List[str]:
        if path_id.startswith("map"):
            path_id = "ko" + path_id[3:]

        if path_id.startswith("ko"):
            ko_genes_txt = self._kegg_rest_api(f"link/ko/{path_id}")
            ko_ids = re.findall(r"(K\d{5})", ko_genes_txt)
            hsa_ids = []
            for i in range(0, len(ko_ids), 100):
                batch = ",".join(ko_ids[i:i+100])
                txt = self._kegg_rest_api(f"link/hsa/{batch}")
                hsa_ids += re.findall(r"hsa:(\d+)", txt)
            hits = self.mg.querymany(hsa_ids, scopes="entrezgene", fields="ensembl.protein", species="human")
        elif path_id.startswith("hsa"):
            data = self.kegg.get(path_id)
            parsed = self.kegg.parse(data)
            genes = parsed.get("GENE", {})
            symbols = {genes[k].split(";")[0].split()[0] for k in genes if k.isdigit()}
            hits = self.mg.querymany(list(symbols), scopes="symbol", fields="ensembl.protein", species="human")
        else:
            raise ValueError("Invalid pathway ID")

        return self._extract_proteins(hits)

    def _from_gene(self, symbol: str) -> List[str]:
        if symbol.startswith("9606."):
            return [symbol]
        hits = self.mg.querymany([symbol], scopes="symbol,alias,entrezgene,ensembl.gene,uniprot", fields="ensembl.protein", species="human")
        return self._extract_proteins(hits)

    def _extract_proteins(self, hits) -> List[str]:
        prots = []
        for h in hits:
            ens = h.get("ensembl", {})
            if isinstance(ens, dict): ens = [ens]
            for e in ens:
                ps = e.get("protein")
                if isinstance(ps, str): prots.append(ps)
                elif isinstance(ps, list): prots.extend(ps)
        return [f"9606.{p}" for p in prots if p]

    def _kegg_rest_api(self, endpoint: str, retries: int = 3) -> str:
        from urllib.request import urlopen, Request
        url = f"https://rest.kegg.jp/{endpoint}"
        for i in range(retries):
            try:
                req = Request(url, headers=self.kegg.http_header)
                with urlopen(req) as resp:
                    return resp.read().decode("utf-8")
            except Exception as e:
                log.warning("[KEGG REST] %s → retry %d/%d", e, i+1, retries)
                time.sleep(2)
        raise RuntimeError(f"KEGG REST failed → {endpoint}")

class HerbGraphManager:
    def __init__(self, db_file: str = "./Data/herb_network_Result.db", cache_dir: str = ".cache"):
        self.db_file = db_file
        self.db_url = f"sqlite:///{db_file}"
        self.tmp_gpath = Path(os.path.join(Path(tempfile.gettempdir()), "ppi_tmp_big.gpickle"))
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.g = None
        self.dist = None

    def prepare_graph(self):
        if not self.tmp_gpath.exists():
            with sqlite3.connect(self.db_file) as conn:
                df = pd.read_sql_query("SELECT protein1, protein2 FROM Human_PPI", conn)
            g = nx.Graph()
            g.add_edges_from(df.itertuples(index=False, name=None))
            with open(self.tmp_gpath, "wb") as f:
                pickle.dump(g, f)
            log.info("✅ Graph built and saved.")
        else:
            log.info("♻️ Reusing graph pickle.")
        with open(self.tmp_gpath, "rb") as f:
            self.g = pickle.load(f)
        self._build_distance_matrix()

    def _build_distance_matrix(self):
        cache = self.cache_dir / "human_ppi_dist.npy"
        if cache.exists():
            mat = np.load(cache, mmap_mode="r")
            mapping = {n: i for i, n in enumerate(self.g.nodes())}
        else:
            import igraph as ig
            g_ig = ig.Graph.TupleList(self.g.edges(), directed=False)
            g_ig.vs["name"] = list(self.g.nodes())
            mapping = {n: i for i, n in enumerate(g_ig.vs["name"])}
            mat = build_matrix(g_ig, cache)
        self.dist = DistMatrix(mat, mapping)

    def get_graph(self):
        return self.g

    def get_distance_matrix(self):
        return self.dist

class HerbScreeningPipeline:
    def __init__(self, graph_manager: HerbGraphManager, herb_db: HerbDatabase):
        self.g = graph_manager.get_graph()
        self.dist = graph_manager.get_distance_matrix()
        self.db = herb_db

    def _lcc_nodes(self, nodes: Sequence[str]) -> list[str]:
        present = [n for n in nodes if n in self.g]
        if not present:
            return []
        sub = self.g.subgraph(present)
        biggest = max(nx.connected_components(sub), key=len)
        return list(biggest)

    def _calc_proximity(self, drug_nodes, target_nodes) -> dict[str, float]:
        res = compute_network_distances_GPU(self.g, self.dist, drug_nodes, target_nodes, random_time=500)
        return {"shortest": res["shortest"], "Z_score": res["Z_score"]["z"], "p": res["Z_score"]["p"]}

    def screen_herbs(self, target_nodes: List[str], herb_list: List[tuple[str, str]]) -> pd.DataFrame:
        results = []
        kor_map = {cn: kor for cn, kor in herb_list}
        for cn, _ in tqdm(herb_list, desc="Screening Herbs"):
            try:
                drug_nodes = self._lcc_nodes(self.db.get_ensp_ids(cn))
                print(len(self.db.get_ensp_ids(cn)))
                t_nodes = self._lcc_nodes(target_nodes)
                print(len(t_nodes))
                if not drug_nodes or not t_nodes:
                    raise ValueError("Missing valid nodes")
                out = self._calc_proximity(drug_nodes, t_nodes)
                out["error"] = ""
            except Exception as e:
                out = {"shortest": np.nan, "Z_score": np.nan, "p": np.nan, "error": str(e)}
            out["herb_cn_name"] = cn
            out["herb_kor_name"] = kor_map.get(cn, "")
            results.append(out)

        df = pd.DataFrame(results).dropna(subset=["Z_score"])
        df = df.query("p > 0.05").sort_values("Z_score").reset_index(drop=True)
        return df

    def save_screening_result(self, df: pd.DataFrame, target_val: str, target_type: str,
                              out_dir: str = "Screening_proximity_result") -> str:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        tag = str(target_val).replace(":", "_")

        target_dir = Path(out_dir) / tag / target_type
        target_dir.mkdir(parents=True, exist_ok=True)

        out_path = target_dir / f"{tag}_{target_type}_screen_{ts}.csv"
        df.to_csv(out_path, index=False, encoding="utf-8-sig")
        log.info("✅ Saved screening result → %s", out_path)
        return str(out_path)

class EvaluationEngine:
    def __init__(self, base_dir: str, p: float = 0.05, z_score: float = -1.00):
        self.base_dir = Path(base_dir)
        self.p = p
        self.z_score = z_score

    def evaluate(self) -> pd.DataFrame:
        level_dirs = ["disease", "pathway", "genes"]
        count_all: list[tuple[str, str]] = []

        out_path = self.base_dir / "screened_herbs.xlsx"
        with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
            for level in level_dirs:
                result_cols: dict[str, list[tuple[str, str]]] = {}
                dir_path = self.base_dir / level
                if not dir_path.is_dir():
                    continue

                files = list(dir_path.glob("*.csv")) + list(dir_path.glob("*.xlsx")) + list(dir_path.glob("*.xls"))

                for fp in files:
                    try:
                        df = pd.read_excel(fp) if fp.suffix in [".xls", ".xlsx"] else pd.read_csv(fp)
                    except Exception:
                        continue

                    if not {"p", "Z_score", "shortest", "herb_cn_name", "herb_kor_name"}.issubset(df.columns):
                        continue

                    for col in ["p", "Z_score", "shortest"]:
                        df[col] = pd.to_numeric(df[col], errors="coerce")
                    df = df.dropna(subset=["p", "Z_score", "shortest"])

                    if level == "Gene Level":
                        filtered = df.sort_values("shortest").head(50)
                    else:
                        filtered = df.query("Z_score < @self.z_score").sort_values("Z_score").head(50)

                    if filtered.empty:
                        continue

                    col_name = fp.stem.split("_")[0]
                    while col_name in result_cols:
                        col_name += "_dup"

                    herb_pairs = list(zip(filtered["herb_cn_name"], filtered["herb_kor_name"]))
                    result_cols[col_name] = herb_pairs

                if not result_cols:
                    print(f"⚠ No valid screening results found in {level}")
                    continue

                # Raw 시트 만들기
                max_len = max(map(len, result_cols.values()))
                raw_df_dict = {}
                for col, herb_list in result_cols.items():
                    padded = herb_list + [("", "")] * (max_len - len(herb_list))
                    cn_names = [x[0] for x in padded]
                    kr_names = [x[1] for x in padded]

                    raw_df_dict[f"{col}"] = cn_names
                    raw_df_dict[f"{col} (KR)"] = kr_names

                    count_all.extend(herb_list)

                raw_df = pd.DataFrame(raw_df_dict)
                raw_df.to_excel(writer, sheet_name=level.capitalize(), index=False)

            # Count 시트 만들기 (모든 level 통합)
            count_df = pd.DataFrame([pair for pair in count_all if all(pair)],
                                    columns=["herb_cn_name", "herb_kor_name"])
            count_df = count_df.value_counts().reset_index(name="Count")
            count_df.to_excel(writer, sheet_name="Count", index=False)

        print(f"[SAVE] Evaluation results saved → {out_path}")
        return None
