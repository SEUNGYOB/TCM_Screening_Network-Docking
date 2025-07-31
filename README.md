Prerequisite




1. https://drive.google.com/file/d/19pyxIsQkW7GVuG1Tff0JwPs8fq3pb3pN/view?usp=sharing 에서 "data.zip" download
  ./Data 안에 Disease_data folder 및 내용물, herb_network_Result.db, herbs_detailed_database.db, herbs_list.csv

2. .cache -> human_ppi_dist.npy (Human PPI 안의 노드 간의 모든 거리 계산된 캐시) 코드안에서 구동 or data.zip 파일 안에서 폴더 이동(.cache로)

    herb_graph = HerbGraphManager()
    herb_graph.save_distance_cache_from_db()

****3. execution_file.py 사용방법
**   **targe**t_val = "MONDO_0021187"
   ** run_screening(target_val, "disease")****
********
