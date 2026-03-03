# Multiome at 96 and 120h

## Fastqs to matrices/loom and fragments/peaks

This analysis was computed using a local [galaxy](https://doi.org/10.1093/nar/gkac247) server.

scRNAseq and scATACseq were processed independently:

### scRNAseq

2 workflows have been used:
1. The workflow for the count matrices has been exported [here](./Galaxy-Workflow-scRNA-seq_preprocessing_10X_v3_Bundle.ga). It has been run with the following inputs:

- fastq PE collection : will be available on SRA/GEO
- gtf: available on [zenodo](https://zenodo.org/records/10079673)
- cellranger_barcodes_3M-february-2018.txt: [barcodes specific to multiome](./737K-arc-v1_GEX.txt)
- reference genome: mm10
- Barcode Size is same size of the Read: False

2. The workflow to get the loom file from velocyto has been exported [here](../scRNAseq/Galaxy-Workflow-Velocyto-on10X.ga). It has been launched with the following inputs:

- BAM files with CB and UB: the BAM from the previous workflow (RNA STARSolo on collection 13: Alignments)
- filtered matrices in bundle: the filtered matrice from the previous workflow
- gtf file: same as above

All the command lines used have been written [here](./get_loom.sh).

### scATACseq

The workflow used has been exported [here](./Galaxy-Workflow-MultiomeATAC-seq_part.ga). It has been run with the following inputs:

- fastq file containing barcode sequences (24 nt, I2): will be available on SRA/GEO
- fastq R1: will be available on SRA/GEO
- fastq R2: will be available on SRA/GEO
- Genome: mm10
- effective_genome_size: 1870000000
- RNA_737K-arc-v1.txt.gz: [RNA_737K-arc-v1.txt.gz](./737K-arc-v1_GEX.txt)
- ATAC_737K-arc-v1.txt.gz: [ATAC_737K-arc-v1.txt.gz](./737K-arc-v1_ATAC.txt)

All the command lines used have been written [here](./get_atac_peaks.sh).

## Matrices to figures

This part was run on SCITAS (Scientific IT and Application Support) at EPFL using the docker image lldelisle/verse_with_more_packages:4.4.1_8.

First, the data were retrieved from galaxy and moved to good directories and indexed see [here](./Multiome.Step0.sh) for the command lines.

For 96h analysis
```bash
sbatch Multiome/96h/Multiome.Step1.SCITAS.sh
sbatch Multiome/96h/Multiome.Step2.SCITAS.sh
sbatch Multiome/96h/Multiome.Step2.5SCITAS.sh
sbatch Multiome/96h/Multiome.Step3SCITAS.sh
```
