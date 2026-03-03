# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/1.16.8
# command_version:1.16
ln -f -s 'Endodiff_BADC_Endo_rep3_S34_R1_001.fastq.gz.fastqsanger.gz' 'Endodiff_BADC_Endo_rep3_S34_R1_001_fastq_gz.fq.gz' 
 ln -f -s 'Endodiff_BADC_Endo_rep3_S34_R2_001.fastq.gz.fastqsanger.gz' 'Endodiff_BADC_Endo_rep3_S34_R2_001_fastq_gz.fq.gz' 
 cutadapt -j ${GALAXY_SLOTS:-1} -a 'TruSeqR1'='GATCGGAAGAGCACACGTCTGAACTCCAGTCAC' -A 'TruSeq'='GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT' --output='out1.fq.gz' --paired-output='out2.fq.gz' --minimum-length=15 --quality-cutoff=30 'Endodiff_BADC_Endo_rep3_S34_R1_001_fastq_gz.fq.gz' 'Endodiff_BADC_Endo_rep3_S34_R2_001_fastq_gz.fq.gz' > report.txt

# toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.7a
# command_version:
STAR --runThreadN ${GALAXY_SLOTS:-4} --genomeDir '/data/galaxy/galaxy/var/tool-data/rnastar/2.7.4a/mm10_UCSC/mm10_UCSC/dataset_163252_files' --sjdbOverhang 99 --sjdbGTFfile 'mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsOnly_UCSC.gtf.gz.gtf' --readFilesIn 'cutadapt of Endodiff_BADC_Endo_rep3_S34_R1_001.fastq.gz.fastqsanger.gz' 'cutadapt of Endodiff_BADC_Endo_rep3_S34_R1_001.fastq.gz.fastqsanger.gz' --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate '' --quantMode GeneCounts --outSAMattributes NH HI AS nM --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outBAMsortingThreadN ${GALAXY_SLOTS:-4} --limitBAMsortRAM $((${GALAXY_MEMORY_MB:-0}*1000000)) 
 samtools view -b -o 'RNA STAR on cutadapt of Endodiff_BADC_Endo_rep3_S34_R1_001.fastq.gz.bam' Aligned.sortedByCoord.out.bam

# toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_awk_tool/1.1.0
# awk program is  NR>4{print $1,$4} 
# command_version:GNU Awk 4.1.3, API: 1.1
awk --sandbox -v FS='	' -v OFS='	' --re-interval -f "/data/galaxy/galaxy/jobs/000/149/149800/configs/tmp0n7zn2pc" "RNA STAR on data 49, data 727, and data 726: reads per gene.tabular" > "htseqCountFormat for RNA STAR on data 49, data 727, and data 726: reads per gene.tabular"

# toolshed.g2.bx.psu.edu/repos/devteam/bamtools_filter/bamFilter/2.4.1
# The filter is tag NH=1 
# command_version:
cp '/data/galaxy/galaxy/jobs/000/149/149799/configs/tmp8yeaqcs_' 'Filter on data 731: JSON filter rules.txt' 
 ln -s 'RNA STAR on cutadapt of Endodiff_BADC_Endo_rep3_S34_R1_001.fastq.gz.bam' localbam.bam 
 ln -s '/data/galaxy/data/_metadata_files/027/metadata_27276.dat' localbam.bam.bai 
 cat '/data/galaxy/galaxy/jobs/000/149/149799/configs/tmp8yeaqcs_' 
 bamtools filter -script '/data/galaxy/galaxy/jobs/000/149/149799/configs/tmp8yeaqcs_' -in localbam.bam -out 'uniquely mapped of RNA STAR on cutadapt of Endodiff_BADC_Endo_rep3_S34_R1_001.fastq.gz.bam'

# toolshed.g2.bx.psu.edu/repos/devteam/cufflinks/cufflinks/2.2.1.3
# command_version:cufflinks v2.2.1
cufflinks -q --no-update-check 'uniquely mapped of RNA STAR on cutadapt of Endodiff_BADC_Endo_rep3_S34_R1_001.fastq.gz.bam' --num-threads "${GALAXY_SLOTS:-4}" -G 'mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.102_ExonsOnly_UCSC.gtf.gz.gtf' 2> stderr 
 python '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/devteam/cufflinks/d080005cffe1/cufflinks/mass.py' stderr 'None' "transcripts.gtf"
