# Cheng et al. *Nat Commun* 2021

This repo contains all the data and scripts required to generate the TRAP-seq figure as presented in Cheng et al. *Nat Commun* 2021 (see below for full citation). The count data and metadata is hosted on GEO ([GSE176202](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176202)). The summary data from GEO is already found in the **data** folder.

## Files
```
+-- data
    +-- GSE176202_counts.tsv.gz
    +-- GSE176202_metadata.csv.gz
+-- scripts
   +-- analysis.R
+-- environment.yml
+-- README.md
```

## Software
This pipeline is best run with [Anaconda](https://www.anaconda.com/products/individual) to replicate the original environment. Alternatively, you can assemble your own packages and specific versions based on the **environment.yml** file.

## Pipeline
1. Create the conda environment: `conda create -n cheng-natcommun-2021 --file environment.yml`
2. Activate the conda environment: `conda activate cheng-natcommun-2021`
3. Generate the figure from the count data: `R -e "source('scripts/analysis.R')"`
4. Deactivate the conda environment: `conda deactivate`

## Output
* **figures/heatmap.pdf**: plot of enriched genes
* **figures/zscores.xlsx**: raw data from **figures/heatmap.pdf**

## Citation
(To come)
