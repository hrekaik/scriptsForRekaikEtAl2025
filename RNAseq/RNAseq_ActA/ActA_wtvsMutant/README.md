# RNA-seq analysis of Activin A gastruloids

The samples are described [here](./samplesPlan_ActA_wtvsBADC.txt).

The scripts used are available [here](https://github.com/lldelisle/rnaseq_rscripts/tree/39049c96ae9934985372ab49ef4b3791fbc8fe0f).

## Merge tables

Individual counts and FPKM tables were merged in order to have a single table with all samples.

```
Rscript /RNAseq/RNAseq_ActA/ActA_wtvsMutant/configFileRNAseq_step1.R
```

The tables are available [here](../../../outputs/RNAseq/RNAseq_ActA/ActA_wtvsMutant/).

## DESeq2 analysis

DESeq2 and volcano analysis was run on 96h and 120h counts separately:

```
For 96h
Rscript /RNAseq/RNAseq_ActA/ActA_wtvsMutant/volcano_96h.R

For 120h
Rscript /RNAseq/RNAseq_ActA/ActA_wtvsMutant/volcano_120h.R
```

## PCA analysis

PCA plot was obtained with:

```
Rscript /RNAseq/RNAseq_ActA/ActA_wtvsMutant/PCA_ActivinA.R
```

## Heatmaps

Heatmaps plots were obtained with:

```
Rscript /RNAseq/RNAseq_ActA/ActA_wtvsMutant/heatmap_ActA_wtvsBADC.R
```

The output is [here](../../../outputs/RNAseq/RNAseq_ActA/ActA_wtvsMutant/).
