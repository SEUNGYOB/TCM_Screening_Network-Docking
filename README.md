# How to use this code

## Prerequisite

### 1. https://drive.google.com/file/d/19pyxIsQkW7GVuG1Tff0JwPs8fq3pb3pN/view?usp=sharing 에서 "data.zip" download
  ./Data 안에 Disease_data folder 및 내용물, herb_network_Result.db, herbs_detailed_database.db, herbs_list.csv

### 2. Human PPI 안의 노드 간의 모든 거리 계산된 캐시 확보(다운로드 or 직접 계산 코드)
  .cache 디렉토리 안에 human_ppi_dist.npy 파일 필요
- 2-1. 코드안에서 구동

```
from prox_screeninig_FInal import HerbGraphManager

  herb_graph = HerbGraphManager()
  herb_graph.save_distance_cache_from_db()
```
- 2-2. 1에서 다운받은 data.zip 활용
  
  data 디렉토리에서 .cache로 이동

 

# 3. execution_file.py 
  ## 3-1 "disease" 레벨에서 screening 원하는 경우
    target_val = "MONDO_0021187"
    run_screening(target_val, "disease")

  ## 3-2 "pathway" 레벨에서 screening 원하는 경우
    pathway_list =["hsa04979"]
    for pathway in pathway_list:
      run_screening(pathway, "pathway")

  ## 3-3 "genes" 레벨에서 screening 원하는 경우
    gene_list = ["HMGCR"]
    for gene in gene_list:
      run_screening(gene, "genes")
