# RNA-seq analysis of Chiron and Activin A gastruloids

The samples are described [here](./samplesPlan_ActAvsChiron_wtonly_all.txt).

The scripts used are available [here](https://github.com/lldelisle/rnaseq_rscripts/tree/39049c96ae9934985372ab49ef4b3791fbc8fe0f).

## Merge tables

Individual counts and FPKM tables were merged in order to have a single table with all samples.

```
Rscript /RNAseq/RNAseq_ActA/ActAvsChrion/configFileRNAseq_step1.R
```

The tables are available [here](../../../outputs/RNAseq/RNAseq_ActA/ActAvsChrion/).

## Heatmaps

Heatmaps plots were obtained with:

```
Rscript /RNAseq/RNAseq_ActA/ActAvsChrion/heatmap_ActAvsChiron_wtOnly.R
```

The output is [here](../../../outputs/RNAseq/RNAseq_ActA/ActAvsChrion/).
