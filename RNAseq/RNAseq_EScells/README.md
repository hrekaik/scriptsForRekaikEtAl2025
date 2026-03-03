# RNA-seq analysis of ES cells

The samples are described [here](./samplesPlan_HoxLess_ES.txt).

The scripts used are available [here](https://github.com/lldelisle/rnaseq_rscripts/tree/39049c96ae9934985372ab49ef4b3791fbc8fe0f).

## Merge tables

Individual counts and FPKM tables were merged in order to have a single table with all samples.

```
Rscript ~/Documents/mygit/rnaseq_rscripts/step1-generateTables.R RNAseq/RNAseq_EScells/configFileRNAseq_step1.R
```

The tables have been gzipped [here](../../outputs/RNAseq/RNAseq_EScells/).

## DESeq2 analysis

DESeq2 analysis was run on all counts:

```
Rscript ~/Documents/mygit/rnaseq_rscripts/step2-DESeq2.R RNAseq/RNAseq_EScells/configFileRNAseq_step2.R
```

And the volcano was obtained with:

```
Rscript RNAseq/RNAseq_EScells/volcano.R
```

The output is [here](../../outputs/RNAseq/RNAseq_EScells/).
