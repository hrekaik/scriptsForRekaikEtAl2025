#!/bin/sh
#SBATCH --time=08:00:00
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH -o scriptsForRekaikEtAl2025/scRNAseq/Demultiplex.out
#SBATCH -e scriptsForRekaikEtAl2025/scRNAseq/Demultiplex.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --job-name demultiplex
#SBATCH --chdir /scratch/ldelisle/rstudio_test
# customize --output path as appropriate (to a directory readable only by the user!)

sif_version=4.4.1_8
SIF=$PWD/verse_with_more_packages_${sif_version}.sif

if [ ! -e $SIF ]; then
   export APPTAINER_CACHEDIR=$PWD/.cache
   singularity pull docker://lldelisle/verse_with_more_packages:${sif_version}
fi
export APPTAINER_BIND="/scratch/$(id -un)/"
singularity exec --cleanenv $SIF \
    Rscript $PWD/scriptsForRekaikEtAl2025/scRNAseq/Step1.Seurat.Demultiplexing.Analysis.R
