# scRNA-seq integration

The scRNA-seq dataset is integrated with E6.5 to E9.5 from TOME https://tome.gs.washington.edu.

This part was run on SCITAS (Scientific IT and Application Support) at EPFL using the docker image lldelisle/verse_with_more_packages:4.4.1_8.

```bash
sbatch scRNAseq/integration/Atlas.step1.SCITAS.sh
sbatch scRNAseq/integration/Atlas.step2.SCITAS.sh
```
