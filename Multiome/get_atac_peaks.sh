Visiting 5acf6090577ce616
Visiting 1f8a44effea0981f
Visiting fcd8b3720d73561a
Visiting 3835a2cb7caf96b8
Visiting b2309f1057b41940
Visiting c0194da2036ddfc1
Visiting 54fc69243cbfd28f
Visiting 79c28865ca5e36bd
Visiting f1ec0e16e38e27d3
Visiting 42e8e0e9600e8bf6
Visiting 85c417f66a4a27a8
Visiting c92f8ff16f868d46
Visiting d6cbdfde4cd543b3
Visiting d81556b8c300598a
Visiting 14b4b95793897346
Visiting 7646057bb3541973
Visiting 8cf47a6ba6361c9b
Visiting 0894d9448313404c
Visiting d2479f860420bd72
Visiting 374d0a85c985d028
Visiting ba159b172f403113
Visiting 84d3633657770443
Visiting 45d9ad7b3b3a5ff9
Visiting 8b0d1adf0f91e787
Visiting 4d3832a7faa71fce
Visiting 983ae705fa45d06b
Visiting 480a7961b3073fa5
Visiting c7a54acf828dc950
Visiting 4d3281afb870c16b
Visiting 41ac9d26aca899fa
Visiting 5f1143f6e0fffcb4
Visiting 0ce6add3abe7e5d6
Visiting def63c377a6188c9
Visiting a23f54e41a690d60
Visiting 4cf9cf4ae4105632
Visiting 2df24a6abea397ee
Visiting d5dee84b34101153
Visiting 97f0e93c657f60d6
Visiting 61daf7906a0e2afa
Visiting 0f700938f7ac81d7
Visiting 4faac172bfb874c1
Visiting 76c57ca8938cca8f


Command-lines:


# addValue
# command_version:
perl '/data/galaxy/galaxy/server/tools/filters/fixedValueColumn.pl' 'Sorted fragments.tabular' 'Add column on data 129.tabular' '-' no

# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/4.9+galaxy1
# command_version:4.9
ln -f -s 'Sinto barcode on data 15, data 12, and data 91: barcoded read 1.fastqsanger.gz' 'BADC_1.fq.gz' 
 ln -f -s 'Sinto barcode on data 15, data 12, and data 91: barcoded read 2.fastqsanger.gz' 'BADC_2.fq.gz' 
  cutadapt  -j=${GALAXY_SLOTS:-4}   -a 'Nextera R1'='CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'    -A 'Nextera R2'='CTGTCTCTTATACACATCTGACGCTGCCGACGA'    --error-rate=0.1 --times=1 --overlap=3    --action=trim   --quality-cutoff=30       --minimum-length=20      -o 'out1.fq.gz' -p 'out2.fq.gz'  'BADC_1.fq.gz' 'BADC_2.fq.gz'  > report.txt

# toolshed.g2.bx.psu.edu/repos/iuc/sinto_barcode/sinto_barcode/0.10.1+galaxy0
# command_version:sinto 0.10.1
ln -s 'seqtk_seq on data 9.fastqsanger.gz' barcodes.fastq.gz 
 ln -s 'Hox_Del_7_multi_BADC_ATAC_S2_R1_001.fastq.gz.fastqsanger.gz' read1.fastq.gz 
 ln -s 'Hox_Del_7_multi_BADC_ATAC_S2_R2_001.fastq.gz.fastqsanger.gz' read2.fastq.gz 
 sinto barcode --barcode_fastq barcodes.fastq.gz --read1 read1.fastq.gz --read2 read2.fastq.gz --bases 16

# toolshed.g2.bx.psu.edu/repos/iuc/seqtk/seqtk_seq/1.3.3
# command_version:
seqtk seq -q 0 -X 255 -n '0' -l 0 -Q 33 -s 11 -f 1.0 -L 0  -r      'Hox_Del_7_multi_BADC_ATAC_S2_I2_001.fastq.gz.fastqsanger.gz'  | pigz -p ${GALAXY_SLOTS:-1} --no-name --no-time > 'seqtk_seq on data 9.fastqsanger.gz'

# toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.18
# command_version:
set -o | grep -q pipefail 
 set -o pipefail;     bwa mem  -t "${GALAXY_SLOTS:-1}" -v 1                 '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/bwa_mem_index/mm10_UCSC/mm10_UCSC.fa' 'Cutadapt on data 97 and data 96: Read 1 Output.fastqsanger.gz' 'Cutadapt on data 97 and data 96: Read 2 Output.fastqsanger.gz'  | samtools sort -@${GALAXY_SLOTS:-2} -T "${TMPDIR:-.}" -O bam -o 'Map with BWA-MEM on data 105 and data 104 (mapped reads in BAM format).bam'

# toolshed.g2.bx.psu.edu/repos/iuc/sinto_fragments/sinto_fragments/0.9.0+galaxy1
# command_version:sinto 0.9.0
ln -s 'Map with BWA-MEM on data 105 and data 104 (mapped reads in BAM format).bam' 'input.bam' 
 ln -s '/data/galaxy/data/_metadata_files/026/metadata_26948.dat' 'input.bam.bai' 
 sinto fragments --bam 'input.bam' --min_mapq 30 --barcodetag 'CB' --barcode_regex '[^:]*' --max_distance 5000 --min_distance 10 --shift_plus 4 --shift_minus -5  --fragments Sinto fragments on data 109: fragments BED.bed --nproc "${GALAXY_SLOTS:-1}"

# toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_easyjoin_tool/9.3+galaxy1
# command_version:join (GNU coreutils) 9.3
cp '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/86755160afbf/text_processing/sort-header' ./ 
 chmod +x sort-header 
 perl '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/86755160afbf/text_processing/easyjoin'  -t $'\t'  -e '0' -o auto  -1 '1' -2 '4' 'Paste on data 114 and data 113.tabular' 'Sinto fragments on data 109: fragments BED.bed' > 'Join on data 112 and data 115.tabular'

# Cut1
# command_version:
perl '/data/galaxy/galaxy/server/tools/filters/cutWrapper.pl' 'Join on data 112 and data 115.tabular' 'c3-c5,c2,c6' T 'Cut on data 120.tabular'

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_sortbed/2.30.0+galaxy2
# command_version:bedtools v2.30.0
sortBed -i 'Cut on data 120.tabular'   -g '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len'  > 'Sorted fragments.tabular'

# addValue
# command_version:
perl '/data/galaxy/galaxy/server/tools/filters/fixedValueColumn.pl' 'Sorted fragments.tabular' 'Add column on data 129.tabular' '+' no

# cat1
# command_version:
python /data/galaxy/galaxy/server/tools/filters/catWrapper.py 'Concatenate datasets on data 132 and data 135.tabular' 'Add column on data 129.tabular' 'Add column on data 129.tabular'

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_sortbed/2.30.0+galaxy2
# command_version:bedtools v2.30.0
sortBed -i 'Concatenate datasets on data 132 and data 135.tabular'   -g '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len'  > 'SortBed on Concatenate datasets on data 132 and data 135.tabular'

# toolshed.g2.bx.psu.edu/repos/iuc/macs2/macs2_callpeak/2.2.9.1+galaxy0
# command_version:macs2 2.2.9.1
export PYTHON_EGG_CACHE=`pwd` 
   (macs2 callpeak   -t 'SortBed on Concatenate datasets on data 132 and data 135.tabular'  --name BADC    --format BED   --gsize '1870000000'             --keep-dup 'all'  --d-min 20 --buffer-size 100000  --bdg  --qvalue '0.05'  --nomodel --extsize '200' --shift '-100'  2>&1 > macs2_stderr) 
 cp BADC_peaks.xls 'MACS2 callpeak on data 141 (Peaks in tabular format).tabular'   
 exit_code_for_galaxy=$? 
 cat macs2_stderr 2>&1 
 (exit $exit_code_for_galaxy)

# addValue
# command_version:
perl '/data/galaxy/galaxy/server/tools/filters/fixedValueColumn.pl' 'Sorted fragments.tabular' 'Add column on data 128.tabular' '-' no

# toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/4.9+galaxy1
# command_version:4.9
ln -f -s 'Sinto barcode on data 14, data 11, and data 90: barcoded read 1.fastqsanger.gz' 'wt_1.fq.gz' 
 ln -f -s 'Sinto barcode on data 14, data 11, and data 90: barcoded read 2.fastqsanger.gz' 'wt_2.fq.gz' 
  cutadapt  -j=${GALAXY_SLOTS:-4}   -a 'Nextera R1'='CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'    -A 'Nextera R2'='CTGTCTCTTATACACATCTGACGCTGCCGACGA'    --error-rate=0.1 --times=1 --overlap=3    --action=trim   --quality-cutoff=30       --minimum-length=20      -o 'out1.fq.gz' -p 'out2.fq.gz'  'wt_1.fq.gz' 'wt_2.fq.gz'  > report.txt

# toolshed.g2.bx.psu.edu/repos/iuc/sinto_barcode/sinto_barcode/0.10.1+galaxy0
# command_version:sinto 0.10.1
ln -s 'seqtk_seq on data 8.fastqsanger.gz' barcodes.fastq.gz 
 ln -s 'Hox_Del_7_multi_wt_ATAC_S1_R1_001.fastq.gz.fastqsanger.gz' read1.fastq.gz 
 ln -s 'Hox_Del_7_multi_wt_ATAC_S1_R2_001.fastq.gz.fastqsanger.gz' read2.fastq.gz 
 sinto barcode --barcode_fastq barcodes.fastq.gz --read1 read1.fastq.gz --read2 read2.fastq.gz --bases 16

# toolshed.g2.bx.psu.edu/repos/iuc/seqtk/seqtk_seq/1.3.3
# command_version:
seqtk seq -q 0 -X 255 -n '0' -l 0 -Q 33 -s 11 -f 1.0 -L 0  -r      'Hox_Del_7_multi_wt_ATAC_S1_I2_001.fastq.gz.fastqsanger.gz'  | pigz -p ${GALAXY_SLOTS:-1} --no-name --no-time > 'seqtk_seq on data 8.fastqsanger.gz'

# toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.18
# command_version:
set -o | grep -q pipefail 
 set -o pipefail;     bwa mem  -t "${GALAXY_SLOTS:-1}" -v 1                 '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/bwa_mem_index/mm10_UCSC/mm10_UCSC.fa' 'Cutadapt on data 95 and data 94: Read 1 Output.fastqsanger.gz' 'Cutadapt on data 95 and data 94: Read 2 Output.fastqsanger.gz'  | samtools sort -@${GALAXY_SLOTS:-2} -T "${TMPDIR:-.}" -O bam -o 'Map with BWA-MEM on data 102 and data 101 (mapped reads in BAM format).bam'

# toolshed.g2.bx.psu.edu/repos/iuc/sinto_fragments/sinto_fragments/0.9.0+galaxy1
# command_version:sinto 0.9.0
ln -s 'Map with BWA-MEM on data 102 and data 101 (mapped reads in BAM format).bam' 'input.bam' 
 ln -s '/data/galaxy/data/_metadata_files/026/metadata_26949.dat' 'input.bam.bai' 
 sinto fragments --bam 'input.bam' --min_mapq 30 --barcodetag 'CB' --barcode_regex '[^:]*' --max_distance 5000 --min_distance 10 --shift_plus 4 --shift_minus -5  --fragments Sinto fragments on data 108: fragments BED.bed --nproc "${GALAXY_SLOTS:-1}"

# Paste1
# command_version:
perl '/data/galaxy/galaxy/server/tools/filters/pasteWrapper.pl' 'ATAC_737K-arc-v1.txt.gz.txt' 'RNA_737K-arc-v1.txt.gz.txt' T 'Paste on data 114 and data 113.tabular'

# toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_easyjoin_tool/9.3+galaxy1
# command_version:join (GNU coreutils) 9.3
cp '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/86755160afbf/text_processing/sort-header' ./ 
 chmod +x sort-header 
 perl '/data/galaxy/galaxy/var/shed_tools/toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/86755160afbf/text_processing/easyjoin'  -t $'\t'  -e '0' -o auto  -1 '1' -2 '4' 'Paste on data 114 and data 113.tabular' 'Sinto fragments on data 108: fragments BED.bed' > 'Join on data 111 and data 115.tabular'

# Cut1
# command_version:
perl '/data/galaxy/galaxy/server/tools/filters/cutWrapper.pl' 'Join on data 111 and data 115.tabular' 'c3-c5,c2,c6' T 'Cut on data 119.tabular'

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_sortbed/2.30.0+galaxy2
# command_version:bedtools v2.30.0
sortBed -i 'Cut on data 119.tabular'   -g '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len'  > 'Sorted fragments.tabular'

# addValue
# command_version:
perl '/data/galaxy/galaxy/server/tools/filters/fixedValueColumn.pl' 'Sorted fragments.tabular' 'Add column on data 128.tabular' '+' no

# cat1
# command_version:
python /data/galaxy/galaxy/server/tools/filters/catWrapper.py 'Concatenate datasets on data 131 and data 134.tabular' 'Add column on data 128.tabular' 'Add column on data 128.tabular'

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_sortbed/2.30.0+galaxy2
# command_version:bedtools v2.30.0
sortBed -i 'Concatenate datasets on data 131 and data 134.tabular'   -g '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len'  > 'SortBed on Concatenate datasets on data 131 and data 134.tabular'

# toolshed.g2.bx.psu.edu/repos/iuc/macs2/macs2_callpeak/2.2.9.1+galaxy0
# command_version:macs2 2.2.9.1
export PYTHON_EGG_CACHE=`pwd` 
   (macs2 callpeak   -t 'SortBed on Concatenate datasets on data 131 and data 134.tabular'  --name wt    --format BED   --gsize '1870000000'             --keep-dup 'all'  --d-min 20 --buffer-size 100000  --bdg  --qvalue '0.05'  --nomodel --extsize '200' --shift '-100'  2>&1 > macs2_stderr) 
 cp wt_peaks.xls 'MACS2 callpeak on data 140 (Peaks in tabular format).tabular'   
 exit_code_for_galaxy=$? 
 cat macs2_stderr 2>&1 
 (exit $exit_code_for_galaxy)

# toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_cat/9.3+galaxy1
# command_version:cat (GNU coreutils) 9.3
cat 'Peaks from fragments.bed' >> 'Concatenate datasets on data 151 and data 147.bed' 
 cat 'Peaks from fragments.bed' >> 'Concatenate datasets on data 151 and data 147.bed' 
 exit 0

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_sortbed/2.31.1+galaxy0
# command_version:bedtools v2.31.1
sortBed -i 'Concatenate datasets on data 151 and data 147.bed'   -g '/data/galaxy/galaxy/var/tool-data/mm10_UCSC/len/mm10_UCSC.len'  > 'SortBed on Concatenate datasets on data 151 and data 147.bed'

# toolshed.g2.bx.psu.edu/repos/iuc/bedtools/bedtools_mergebed/2.31.1
# command_version:bedtools v2.31.1
mergeBed -i 'SortBed on Concatenate datasets on data 151 and data 147.bed'  -d 0    > 'Merged SortBed on Concatenate datasets on data 151 and data 147.bed'
