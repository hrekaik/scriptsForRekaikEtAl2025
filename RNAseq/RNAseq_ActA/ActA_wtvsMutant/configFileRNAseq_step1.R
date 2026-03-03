gitHubDir <- "/Users/hocine.rekaik/Desktop/EXP/Papers/DelHox/Code_figures/RNAseq/RNAseq_ActA/ActA_wtvsMutant"

### Required for all steps ###
RNAseqFunctionPath<-"/Users/hocine.rekaik/Desktop/EXP/Tools/R/RNAseq_analysis/RNAseqFunctions.R"
samplesPlan<-file.path(gitHubDir,"samplesPlan_ActA_wtvsBADC.txt") #This file should be a tabulated file with at least one column called "sample". Optionnaly, the paths to the counts tables and FPKM tables can be provided under the column called: htseq_count_file and cufflinks_file.


#### STEP 1 - MERGE TABLES ### 
#If the merged tables are not already generated:
outputFolderForStep1<-file.path(gitHubDir, "../../../outputs/RNAseq/RNAseq_ActA/ActA_wtvsMutant")
#Needed for DESeq2:
mergeCounts<-T #Do you want to merge counts? T=yes F or commented=no
#Optional: subset the count table
subsetCounts<-F #Do you want to remove some genes from the count table
#Optional:
mergeFPKM<-T
oneLinePerEnsemblID<-T #By default cufflinks split the transcripts which do not overlap in different locus and so different lines, put T if you want to sum the FPKM for non overlapping transcripts (put F if not).
normFPKMWithAnoukMethod<-F #Anouk method: Genes that have the less variable rank should have the same expression.