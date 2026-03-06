# scriptsForRekaikEtAl2025

All scripts necessary to build figures from raw data in Rekaik et al. 2025.

Here is description of each directory:

- [annotations](./annotations/) contains 2 bed files, one with the coordinates of the deletion and one with the Hox genes coordinates from 1 (or 4) to 13.
- [Sangerseq](./Sangerseq/) contains the table and script to generate Sanger sequencing figures of mutant mES cells.
- [ATACseq](./ATACseq/) contains all information about the processing of ATAC-seq datasets as well as how to generate the figure.
- [general](./general/) contains a script that has been used to generate the list of genes affected by the deletion and the genes in chrY.
- [Multiome](./Multiome/) contains all information about the processing of Multiome datasets from fastqs to figures.
- [quantifications](./quantifications/) contains the table and the script used to generate the plot for image-based quantification.
- [ChIPmentation](./ChIPm/) contains all information about the processing of PBX1 ChIPmentation datasets from fastqs to figures.
- [RNAseq](./RNAseq/) contains all information about the processing of RNA-seq datasets (ES cells and ES cells differentiated to endoderm cells) from fastqs to figures.
- [RT-qPCR](./RT-qPCR/) contains the tables with the quantification of the gene expression by RT-qPCR and the script to run the plot.
- [scRNAseq](./scRNAseq/) contains all information about the processing of single-cell RNA-seq datasets from fastqs to figures including the integration with the mouse atlas.

All outputs are in the [outputs directory](./outputs/).
