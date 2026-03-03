# Analysis of PBX1 ChIPmentation

This analysis was computed using a local [galaxy](https://doi.org/10.1093/nar/gkac247) server.

The workflow has been exported [here](./Galaxy-Workflow-ChIPmentation_PE.ga). It has been run with the following inputs:

- fastq R1: will be available on SRA/GEO
- fastq R2: will be available on SRA/GEO
- reference_genome: mm10
- effective_genome_size: 1870000000


## Heatmaps and Venn diagram generation

Plots were obtained with:
```
Rscript ChIPm/analysis_ChIPm.R
```

The output is [here](../outputs/ChIPm/).
