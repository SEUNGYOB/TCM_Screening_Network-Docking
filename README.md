# How to Use This Code

## Prerequisites

### 1. Download Required Data Files
Download the `data.zip` file from the following link:  
[📦 Download data.zip](https://drive.google.com/file/d/19pyxIsQkW7GVuG1Tff0JwPs8fq3pb3pN/view?usp=sharing)

The archive contains:
- `Disease_data/` folder
- `herb_network_Result.db`
- `herbs_detailed_database.db`
- `herbs_list.csv` (list of herb names)

### 2. Prepare Cached Human PPI Distance Matrix
A precomputed all-pairs shortest path matrix is required:  
`.cache/human_ppi_dist.npy`

You can prepare it in two ways:

#### 2-1. Generate the cache using code

```python
from prox_screeninig_FInal import HerbGraphManager

herb_graph = HerbGraphManager()
herb_graph.save_distance_cache_from_db()
```

#### 2-2. Use the downloaded file

Copy the human_ppi_dist.npy from the downloaded data/ directory into your local .cache/ directory.

# 3. execution_file.py 
  ## 3-1 Run screening at the "disease" level
    target_val = "MONDO_0021187"
    run_screening(target_val, "disease")

  ## 3-2 Run screening at the "pathway" level (from KEGG)
    pathway_list =["hsa04979"]
    for pathway in pathway_list:
      run_screening(pathway, "pathway")

  ## 3-3 Run screening at the "genes" level (symbol,alias,entrezgene,ensembl.gene,uniprot)
    gene_list = ["HMGCR"]
    for gene in gene_list:
      run_screening(gene, "genes")
