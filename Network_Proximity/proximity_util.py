#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
(⬇ 2025‑07‑04 디버그 친화 버전)
=====================================================
* **python‑igraph** + **NPY 거리행렬 캐시**로 초고속 proximity 계산
* **디버깅/로깅** 편의성을 위해
  - logging 모듈 기본 설정
  - 단계별 타이머 + 사이즈/노드 정보 출력
  - 캐시 hit / miss 명확 표기
  - `--rebuild-cache` CLI 플래그 추가
"""
from __future__ import annotations
from pathlib import Path
import logging
import time
import igraph as ig
from concurrent.futures import ThreadPoolExecutor, as_completed
import os, random, numpy as np
from typing import Sequence, Dict, Any
import networkx as nx

###############################################################################
# 로깅 설정
###############################################################################
logging.basicConfig(
    format="[%(levelname)s] %(message)s",
    level=logging.INFO,
)
log = logging.getLogger(__name__)

###############################################################################
# 상수 & 헬퍼
###############################################################################
CACHE_DIR = Path(".cache"); CACHE_DIR.mkdir(exist_ok=True)
SENTINEL   = np.uint16(0xFFFF)   # unreachable 거리 표시 (65 535)

###############################################################################
# 0. NetworkX → igraph, 거리행렬 캐시
###############################################################################

def nx_to_igraph(G: nx.Graph):
    g = ig.Graph.TupleList(G.edges(), directed=False)
    g.vs["name"] = list(G.nodes())
    name2idx = {n: i for i, n in enumerate(g.vs["name"])}
    return g, name2idx
SENT = 65535          # DistMatrix에서 쓰는 sentinel 값

def _shortest_mean(mat: np.ndarray, rows: np.ndarray, cols: np.ndarray) -> float:
    """
    A×B 부분행렬 평균. sentinel(SENT) 값은 제외.
    """
    sub = mat[rows][:, cols]          # (|A|, |B|) 뷰
    valid = sub != SENT               # True = 유한 거리
    if not valid.any():               # 완전히 단절된 경우
        return float("inf")
    return float(sub[valid].mean())

def build_matrix(g: ig.Graph, out_file: Path) -> np.ndarray:
    log.info("⚙️  building all‑pairs shortest‑path matrix … (one‑time)")
    t0 = time.perf_counter()
    dist = np.array(g.shortest_paths_dijkstra(), dtype=np.float64)
    dist[np.isinf(dist)] = SENTINEL               # ∞ → sentinel
    dist_u16 = dist.astype(np.uint16)
    np.save(out_file, dist_u16)
    log.info("✅ matrix saved → %s  [%.1f MB, %.1fs]", out_file, dist_u16.nbytes/1e6, time.perf_counter()-t0)
    return dist_u16


def load_graph_assets(
    edgelist: str | Path,
    *,
    matrix_path: str | Path | None = None,
    rebuild: bool = False,
):
    """(G_nx, DistMatrix) 반환. 캐시가 있으면 메모리‑맵, 없으면 생성."""
    edgelist = Path(edgelist)
    log.info("📂 reading edgelist: %s", edgelist)
    G = nx.read_edgelist(edgelist, data=False)
    log.info("   ↳ %d nodes / %d edges", G.number_of_nodes(), G.number_of_edges())

    if matrix_path is None:
        matrix_path = CACHE_DIR / f"{edgelist.stem}_dist.npy"
    else:
        matrix_path = Path(matrix_path)
        matrix_path.parent.mkdir(parents=True, exist_ok=True)

    if matrix_path.exists() and not rebuild:
        log.info("🗄  using cached matrix %s", matrix_path)
        mat = np.load(matrix_path, mmap_mode="r")
        name2idx = {n: i for i, n in enumerate(G.nodes())}
    else:
        log.info("🚧 cache miss → calculating matrix")
        g_ig, name2idx = nx_to_igraph(G)
        mat = build_matrix(g_ig, matrix_path)

    return G, DistMatrix(mat, name2idx)

###############################################################################
# 1. 거리 조회 래퍼
###############################################################################
class DistMatrix:
    """NPY 행렬 기반 최단거리 look‑up (µs)"""

    def __init__(self, mat: np.ndarray, name2idx: Dict[str, int]):
        self.mat = mat
        self.idx = name2idx

    def get(self, u: Any, v: Any) -> float:
        try:
            i, j = self.idx[str(u)], self.idx[str(v)]
        except KeyError:
            return float("inf")
        d = self.mat[i, j]
        return float("inf") if d == SENTINEL else float(d)

###############################################################################
# 2. PairwiseLengths
###############################################################################
class PairwiseLengths:
    def __init__(self, dist: DistMatrix, A: Sequence[Any], B: Sequence[Any]):
        self.d = dist; self.A = list(A); self.B = list(B)

    def build(self):
        len_AB, len_BA, len_AA, len_BB = {}, {}, {}, {}
        for a in self.A:
            len_AB[a] = {b: self.d.get(a, b) for b in self.B}
            len_AA[a] = {x: self.d.get(a, x) for x in self.A if x != a}
        for b in self.B:
            len_BA[b] = {a: len_AB[a][b] for a in self.A}
            len_BB[b] = {y: self.d.get(b, y) for y in self.B if y != b}
        return len_AB, len_BA, len_AA, len_BB

###############################################################################
# 3. NetworkDistance (필요 지표만 유지)
###############################################################################
class NetworkDistance:
    def __init__(self, tables, A, B):
        self.len_AB, self.len_BA, self.len_AA, self.len_BB = tables
        self.A, self.B = A, B

    def _closest(self):
        dA = [min(self.len_AB[a].values()) for a in self.A]
        dB = [min(self.len_BA[b].values()) for b in self.B]
        return float(np.mean(dA + dB))

    def _separation_inner(self, group, table):
        if len(group) <= 1: return 0.0
        return float(np.mean([min(table[n].values()) for n in group]))

    def separation(self):
        return self._closest() - (
            self._separation_inner(self.A, self.len_AA) +
            self._separation_inner(self.B, self.len_BB)
        ) / 2

    def shortest(self):
        return float(np.mean([
            *[d for a in self.A for d in self.len_AB[a].values()],
            *[d for b in self.B for d in self.len_BA[b].values()],
        ]))

###############################################################################
# 4. high‑level API
###############################################################################

# Z-score 계산을 위한 사용───────────────────────────────────────────────────────────

# def compute_network_distances_GPU(G: nx.Graph,
#                                   dist: DistMatrix,
#                                   A: Sequence[Any],
#                                   B: Sequence[Any],
#                                   *, random_time=100,
#                                   seed=42,
#                                   max_workers=None
#                                   ):
#     import os, random
#     from concurrent.futures import ThreadPoolExecutor, as_completed
#
#     dist_mat_np = dist.mat
#     name2idx = dist.idx
#     A = [n for n in A if n in G]
#     B = [n for n in B if n in G]
#     if not A or not B: raise ValueError("노드 없음")
#
#     A_idx = np.asarray([name2idx[n] for n in A], dtype=np.int32)
#     B_idx = np.asarray([name2idx[n] for n in B], dtype=np.int32)
#
#     # ── 1) 원본 거리 ─────────────────────────────────────
#     d_orig = _shortest_mean(dist_mat_np, A_idx, B_idx)
#
#     # ── 2) degree bucket 준비 ───────────────────────────
#     deg2 = {}
#     for n in G:
#         deg2.setdefault(G.degree[n], []).append(name2idx[n])
#
#     rng = random.Random(seed)
#     def one_sample(k):
#         rA = np.fromiter((rng.choice(deg2[G.degree[a]]) for a in A), int, len(A))
#         rB = np.fromiter((rng.choice(deg2[G.degree[b]]) for b in B), int, len(B))
#         return _shortest_mean(dist_mat_np, rA, rB)
#
#     # ── 3) 병렬 실행 (스레드) ───────────────────────────
#     max_workers = max_workers or min(32, (os.cpu_count() or 1)*4)
#     with ThreadPoolExecutor(max_workers=max_workers) as ex:
#         rnd_vals = [f.result()
#                     for f in as_completed(ex.submit(one_sample, i)
#                                           for i in range(random_time))]
#
#     rnd_vals = [v for v in rnd_vals if np.isfinite(v)]  # ← 추가: ∞ 샘플 제거
#     if not rnd_vals:
#         raise RuntimeError("모든 무작위 샘플이 연결되지 않았습니다.")
#
#     mean, std = float(np.mean(rnd_vals)), float(np.std(rnd_vals))
#     z = 0.0 if std == 0 else (d_orig - mean) / std
#     p = sum(v <= d_orig for v in rnd_vals) / len(rnd_vals)  # 표본 수 변경
#
#     return {"shortest": d_orig,
#             "Z_score": {"d": d_orig,"z":z,"mean":mean,"std":std,"p":p}}

def compute_network_distances_GPU(
    G: nx.Graph,
    dist: DistMatrix,
    A: Sequence[Any],
    B: Sequence[Any],
    *,
    min_bin_size: int = 50,
    random_time: int = 100,
    seed: int = 42,
    max_workers: int = None
):
    import random, os
    from concurrent.futures import ThreadPoolExecutor, as_completed

    dist_mat_np = dist.mat
    name2idx = dist.idx
    A = [n for n in A if n in G]
    B = [n for n in B if n in G]
    if not A or not B:
        raise ValueError("노드 없음")

    A_idx = np.asarray([name2idx[n] for n in A], dtype=np.int32)
    B_idx = np.asarray([name2idx[n] for n in B], dtype=np.int32)

    # ── 1) 원본 거리 ─────────────────────────────────────
    d_orig = _shortest_mean(dist_mat_np, A_idx, B_idx)

    # ── 2) min_bin_size 기반 bin 생성 ───────────────────
    deg_to_nodes = {}
    for node, deg in G.degree():
        deg_to_nodes.setdefault(deg, []).append(node)

    sorted_degrees = sorted(deg_to_nodes)
    bins = []
    i = 0
    while i < len(sorted_degrees):
        bin_nodes = []
        low = sorted_degrees[i]
        while len(bin_nodes) < min_bin_size and i < len(sorted_degrees):
            bin_nodes.extend(deg_to_nodes[sorted_degrees[i]])
            i += 1
        high = sorted_degrees[i-1]
        if bins and len(bin_nodes) < min_bin_size:
            # 마지막 bin과 병합
            low_prev, high_prev, nodes_prev = bins.pop()
            bins.append((low_prev, high, nodes_prev + bin_nodes))
        else:
            bins.append((low, high, bin_nodes))

    def get_equivalent_nodes(node):
        node_deg = G.degree[node]
        for low, high, members in bins:
            if low <= node_deg <= high:
                return [n for n in members if n != node]
        return []

    rng = random.Random(seed)

    def one_sample(_):
        try:
            rA = [rng.choice(get_equivalent_nodes(a)) for a in A]
            rB = [rng.choice(get_equivalent_nodes(b)) for b in B]
            rA_idx = [name2idx[n] for n in rA if n in name2idx]
            rB_idx = [name2idx[n] for n in rB if n in name2idx]
            if not rA_idx or not rB_idx:
                return float("inf")
            return _shortest_mean(dist_mat_np, rA_idx, rB_idx)
        except:
            return float("inf")

    # ── 3) 병렬 실행 ────────────────────────────────────
    max_workers = max_workers or min(32, (os.cpu_count() or 1) * 4)
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        rnd_vals = [f.result() for f in as_completed(
            ex.submit(one_sample, i) for i in range(random_time)
        )]

    rnd_vals = [v for v in rnd_vals if np.isfinite(v)]
    if not rnd_vals:
        raise RuntimeError("모든 무작위 샘플이 연결되지 않았습니다.")

    mean, std = float(np.mean(rnd_vals)), float(np.std(rnd_vals))
    z = 0.0 if std == 0 else (d_orig - mean) / std
    p = sum(v <= d_orig for v in rnd_vals) / len(rnd_vals)

    return {
        "shortest": d_orig,
        "Z_score": {
            "d": d_orig,
            "z": z,
            "mean": mean,
            "std": std,
            "p": p
        }
    }


