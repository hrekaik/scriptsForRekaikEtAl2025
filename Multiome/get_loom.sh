Visiting 08b713d785d3147b
Visiting 850b9cbdafed5d93
Visiting fc0e0a5540e818e1
Visiting f00012b3bd1fdad9
Visiting 053d8a605e137ee2
Visiting bbe033cb40ea3d12
Visiting e1b88a4b7e0ece9b
Visiting 9fb7870cf0e2f9b8
Visiting d498d96e4d1e33d5
Visiting bea836864dae108f


Command-lines:


# toolshed.g2.bx.psu.edu/repos/iuc/rna_starsolo/rna_starsolo/2.7.11a+galaxy0
# command_version:
STAR  --runThreadN ${GALAXY_SLOTS:-4} --genomeLoad NoSharedMemory --genomeDir '/data/galaxy/galaxy/var/tool-data/rnastar/2.7.4a/mm10_UCSC/mm10_UCSC/dataset_163252_files' --sjdbOverhang 100 --sjdbGTFfile 'bf5a119_mm10_allGastruloids_min10_extended.gtf.gz.gtf' --sjdbGTFfeatureExon 'exon'   --soloType CB_UMI_Simple   --readFilesIn Hox_Del_7_multi_wt_GEX_S3_R2_001.fastq.gz.fastqsanger.gz Hox_Del_7_multi_wt_GEX_S3_R1_001.fastq.gz.fastqsanger.gz --soloCBmatchWLtype 1MM_multi  --readFilesCommand zcat   --soloCBwhitelist '737K-arc-v1_GEX.txt.txt' --soloBarcodeReadLength 0 --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12   --soloStrand Forward --soloFeatures Gene --soloUMIdedup 1MM_CR --soloUMIfiltering - --quantMode TranscriptomeSAM GeneCounts --outSAMattributes NH HI AS nM GX GN CB UB --outSAMtype BAM SortedByCoordinate  --soloCellFilter None  --soloOutFormatFeaturesGeneField3 'Gene Expression'  --outSAMunmapped None --outSAMmapqUnique 60  --limitOutSJoneRead 1000 --limitOutSJcollapsed 1000000 --limitSjdbInsertNsj 1000000     
 mv Solo.out/Gene Solo.out/soloFeatures 
 cat <(echo "Barcodes:") Solo.out/Barcodes.stats <(echo "Genes:") Solo.out/soloFeatures/Features.stats > 'RNA STARSolo on data 45, data 12, and others: Barcode/Feature Statistic Summaries.txt'  
 samtools view -b -o 'RNA STARSolo on data 45, data 12, and others: Alignments.bam' Aligned.sortedByCoord.out.bam

# toolshed.g2.bx.psu.edu/repos/iuc/dropletutils/dropletutils/1.10.0+galaxy2
# command_version:
mkdir 'tenx.input' 
 ln -s 'RNA STARSolo on data 45, data 12, and others: Matrix Gene Counts raw.mtx' 'tenx.input/matrix.mtx' 
 ln -s 'RNA STARSolo on data 45, data 12, and others: Genes raw.tsv' 'tenx.input/genes.tsv' 
 ln -s 'RNA STARSolo on data 45, data 12, and others: Barcodes raw.tsv' 'tenx.input/barcodes.tsv' 
  mkdir 'tenx.output' 
  Rscript '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/iuc/dropletutils/a9caad671439/dropletutils/scripts/dropletutils.Rscript' '/data/galaxy/galaxy/jobs/000/141/141575/configs/tmpghxocfu5'

# toolshed.g2.bx.psu.edu/repos/iuc/velocyto_cli/velocyto_cli/0.17.17+galaxy2
# command_version:Matplotlib created a temporary config/cache directory at /data/galaxy/galaxy/jobs/000/148/148313/tmp/matplotlib-9qih9z05 because the default path (/data/galaxy/galaxy/.config/matplotlib) is not a writable directory; it is highly recommended to set the MPLCONFIGDIR environment variable to a writable directory, in particular to speed up the import of Matplotlib and to better support multiprocessing.
velocyto, version 0.17.17
export NUMBA_CACHE_DIR="${TEMP:-/tmp}";  mkdir -p 'Hox_Del_7_multi_wt_GEX_S3_001_fastq/outs/filtered_gene_bc_matrices/whatever/' 
 ln -s 'RNA STARSolo on data 45, data 12, and others: Alignments.bam' 'Hox_Del_7_multi_wt_GEX_S3_001_fastq/outs/possorted_genome_bam.bam' 
 ln -s 'DropletUtils 10X Barcodes on data 55, data 54, and data 56.tsv' 'Hox_Del_7_multi_wt_GEX_S3_001_fastq/outs/filtered_gene_bc_matrices/whatever/barcodes.tsv' 
 velocyto  run10x  -t 'uint16'  --samtools-threads ${GALAXY_SLOTS:-1} --samtools-memory ${GALAXY_MEMORY_MB:-100}  '-vv' 'Hox_Del_7_multi_wt_GEX_S3_001_fastq' 'bf5a119_mm10_allGastruloids_min10_extended.gtf.gz.gtf' 
 mv 'Hox_Del_7_multi_wt_GEX_S3_001_fastq/velocyto/'*.loom 'output.loom'
