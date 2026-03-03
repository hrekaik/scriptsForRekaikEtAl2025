# The command lines to obtain inputs are in commandlines.
# The preprocessing was run on Galaxy

cd /home/ldelisle/mountDuboule/Hocine/Papier/HoxLess/Multiome_v3/inputs_galaxy
# Get RNA part:
unzip BADC.zip
unzip wt.zip 
mv All\ matrices\ in\ cell-ranger\ format/Hox_Del_7_multi_BADC_GEX_S4_001.fastq/ BADC
mv All\ matrices\ in\ cell-ranger\ format/Hox_Del_7_multi_wt_GEX_S3_001.fastq/ wt
rm -r All\ matrices\ in\ cell-ranger\ format

# Index fragments with htslib 1.21
unzip Sorted_BGZIP.zip
tabix -p bed "Convert Interval to BGZIP on collection 127/wt.bgzip"
tabix -p bed "Convert Interval to BGZIP on collection 127/BADC.bgzip"

mv "Convert Interval to BGZIP on collection 127" "Sorted_Fragments"
