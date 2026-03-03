#!/bin/sh
#SBATCH --time=08:00:00
#SBATCH --ntasks=5
#SBATCH --mem=20G
#SBATCH -o scriptsForRekaikEtAl2025/Multiome/Multiome.Step4.out
#SBATCH -e scriptsForRekaikEtAl2025/Multiome/Multiome.Step4.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --job-name step4
#SBATCH --chdir /scratch/ldelisle/rstudio_test

wget -nc "http://datacache.galaxyproject.org/singularity/h/o/homer:4.11--pl5321h9f5acd7_7"

# Homer version 4.11
export APPTAINER_BIND="/scratch/$(id -un)/"

for file in scriptsForRekaikEtAl2025/outputs/Multiome/*_for_Homer.bed; do
    singularity exec --cleanenv 'homer:4.11--pl5321h9f5acd7_7' findMotifsGenome.pl ${file} inputs_galaxy/mm10_UCSC.fa ${file/.bed/}_motifs -size given -len 8,10,12 &> ${file/.bed/}_motifs.log &
done

wait