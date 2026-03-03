# scRNA-seq

## Fastq to matrices/loom

The first steps of single-cell RNA-seq analysis was performed using a local [galaxy](https://doi.org/10.1093/nar/gkac247) server.

2 workflows have been used:
1. The workflow for the count matrices has been exported [here](./scRNA-seq_preprocessing_10X_cellPlex_UPDUB.ga). They have been run with the following inputs:

- fastqPE collection GEX: fastqs are available on GEO/SRA.
- reference genome: `mm10` (from UCSC)
- gtf: available on [zenodo](https://zenodo.org/records/10079673)
- cellranger_barcodes_3M-february-2018.txt: downloaded from [zenodo](https://zenodo.org/record/3457880/files/3M-february-2018.txt.gz).
- fastqPE collection CMO: fastqs are available on GEO/SRA
- cmo_10X_seq.txt: available [here](./cmo_10X_seq.txt)
- sample name and CMO sequence collection: see [CMO_samples](./CMO_samples) directory
- Number of expected cells (used by CITE-seq-Count): `24000`

2. The workflow to get the loom file from velocyto has been exported [here](./Galaxy-Workflow-Velocyto-on10X.ga). It has been launched with the following inputs:

- BAM files with CB and UB: the BAM from the previous workflow (RNA STARSolo on collection 13: Alignments)
- filtered matrices in bundle: the filtered matrice from the previous workflow
- gtf file: same as above

All the command lines used have been written [here](./get_loom.sh).

## Matrices to figures

This part was run on SCITAS (Scientific IT and Application Support) at EPFL using the docker image lldelisle/verse_with_more_packages:4.4.1_8.

First, the data were retrieved from galaxy and moved to good directories and indexed see [here](./Multiome.Step0.sh) for the command lines.

```bash
sbatch scRNAseq/RNAseq.demultiplex.SCITAS.sh
sbatch scRNAseq/RNAseq.step1.SCITAS.sh
sbatch scRNAseq/RNAseq.step2.SCITAS.sh
sbatch scRNAseq/RNAseq.step3.SCITAS.sh
```

The steps to integrate the scRNAseq with the available mouse Atlas are available [here](./integration/README.md).
