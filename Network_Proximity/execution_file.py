from prox_screeninig_FInal import (
    HerbGraphManager,
    HerbDatabase,
    TargetResolver,
    HerbScreeningPipeline,
    EvaluationEngine
)
import pandas as pd

target_val = "MONDO_0021187"

def run_screening(target_val: str, target_type: str,
                  db_path: str = "./Data/herb_network_Result.db",
                  herb_list_path: str = "./Data/herbs_list.csv",
                  out_dir: str = "Screening_proximity_result") -> str:
    """
    단일 target_val을 기준으로 스크리닝을 수행하고 결과를 저장합니다.
    예: run_screening("MONDO_0021187", "disease")
    target_val = "MONDO_0021187
    target_type = ["disease","pathway","genes"]
    """
    # 1. 그래프 준비
    graph_mgr = HerbGraphManager()
    graph_mgr.prepare_graph()

    # 2. DB 및 로더
    herb_db = HerbDatabase(db_path)
    resolver = TargetResolver()

    # 3. 한약 리스트 불러오기
    def _load_herb_LISTS(csv_path) -> list[tuple[str, str]]:
        df = pd.read_csv(csv_path)
        return list(df[["herb_cn_name", "herb_kor_name"]].itertuples(index=False, name=None))

    herbs = _load_herb_LISTS(herb_list_path)
    # print(herbs)
    # 4. 타겟 ID → ENSP 변환
    target_nodes = resolver.resolve(target_type, target_val)
    # print(target_nodes[:3])
    # 5. 스크리닝 실행 및 저장
    pipeline = HerbScreeningPipeline(graph_mgr, herb_db)
    df_result = pipeline.screen_herbs(target_nodes, herbs)

    print(df_result.head(3))
    return pipeline.save_screening_result(df_result, target_val, target_type, out_dir=out_dir)

# run_screening(target_val, "disease")
#7 Evaluating results
base_dir = f"Screening_proximity_result/{target_val}"
evaluator = EvaluationEngine(base_dir)
df_result = evaluator.evaluate()