#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Oleg Moskvin; info@scienceforever.com 
# May 29, 2013 
#########################################################
# 
# This is the default attribute file
# (these attributes are used if no database data exchange is chosen; otherwise, they are overwritten by 
# the values retrieved from the database) 
#
#########################################################
#
# Usage: called from the GLSeq.top.R script;
# multiple versions of local attribute files may exist; the name of the particular version
# may be supplied as an option to the GLSeq.top.R script 
#
#########################################################
#
# directory containing raw files
# (may be non-writable!)
raw.dir <- "/your-directory-with-raw-compressed-files"
#
# Files in the raw dir are normally compressed but may be not:
unzipped <- FALSE
#
# directory contining ready-to-go (split+QC-processed) fq files (Oct 17, 2013):
readyData.dir <- "/your-directory-with-ready-to-go-files"
#
# the directory containing the GLSeq scripts:
base.dir <- "/path-to-GLSeq"
#
# raw file names: 
raw.fNames <- NULL
#
# strain
strain <- ""
#
# reference genome - the index directory for the respective method (RSEM or BWA)
# (must match the name of the  subfolder under base.dir):
rGenome <- "K12"
#
# name of the reference fasta file (may differ from the base name of the reference -
# still, it should be located in the folder for the selected feference): 
refFASTAname <- "E.coli.fa" 
#
# name of the reference genomic features file (may differ from the base name of the reference -
# this will be eventually taken from the database); it should be located in the folder for the selected feference)
# gtf files may be used instead of gff where applicable; the same object is used for both cases:
refGFFname <- "E.coli.gtf"
# refGFFname <- "Novosphingobium_3replicons.Clean.gff"
#
# number of the column in GTF file with the gene / other IDs charachter string (9, unless the file is non-standard for some reason):
gtfFeatureColumn <- 9
#
# GFF attribute to be used as feature ID (HTSeq):
idAttr <- "gene_id" # default
# idAttr <- "locus_tag" # useful if "gene_id" has duplicated entries
#
# single / paired end
paired.end <- FALSE
# 
# quality scores format
qScores <- "phred33" 
#
# number of cores to use
nCores <- 6
#
# number of parallel computation streams for expression computation
nStreams <- 8
#
# number of parallel computation streams for data preparation
# (may differ from the number of streams for expression computation because of particular software demands) 
nStreamsDataPrep <- 4
#
# compute confidence intervals? 
compConf <- TRUE
#
# Start from raw genome? (yes/no); this attribute is de facto unused at the moment
rawGen <- FALSE
#
# quantification algorithm
# qAlgor <- "bwaHTSeq" # Dana pipeline of BWA alignment AND HTseq counting 
qAlgor <- "RSEM" 
#
# Maximal length of fragment (for paired-end libraries)
fragMaxLength <- 1000
#
# Maximal size (MB) of the auxiliary buffer used for computing credibility intervals (CI) - for RSEM (+extra 2Gb per stream)
ciMem <- 4096
#
# Output genome bam
genobam <- TRUE
#
# Strandness of the library (NULL, F, R) 
libstrand <- "R"
#
# Base for the strandness representation ("transcriptome", "genome") - makes no difference at the moment:
strandBase <- "genome" 
#
# Extract explicit forward and reverse coverage from the original BAM file? 
strandExtract <- TRUE
#
# Base of the destination directory (added May 9, 2013)
# This should be located on a FAST volume (SCSI is recommended)
# a particular subfolder names after the run ID will be created by GLSeq below this folder
dest.dir.base <-  "/your-base-path-to-expected-results"
#
# the default run attempt
runAttempt <- formatC(1, width=2, flag="0")
#
# the actual unique run ID - 
# text.add <- paste(expID, runAttempt, sep=".")
# now is being generated inside GLSeq.top.R
#
# Number of unique characters in the beginning of the each file (library ID length):
libNchar <- 36
#
# Subset of the libraries to process (optional; normally the list wil be generated from the actual directory content)
libList <- NULL
#
# date of the run
runDate <-  paste(strsplit(date()," ")[[1]], collapse="_")
runDate <- gsub(":", "_", runDate)
#
# path to Trimmomatic: 
trimPath <- '/opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar'
#
# path to PicardTools jar directory:
# picardToolsPath <- '/soft/picard-tools-1.98/picard-tools-1.98/'
picardToolsPath <- '/opt/bifxapps/picard-tools/'
#
# path to fastqc:
fastqcPath <- '/opt/bifxapps/bin/fastqc' 
#
# path to BWA:
bwaPath <- '/opt/bifxapps/bin/bwa'
#
# trim the reads and generate QC reports for before- and after-trimming FASTQ files? 
readTrim <- TRUE
#
# trimmomatic parameter values for HEADCROP  
trimHead <- 12
#
# name of the FASTA file with artificail sequences (adapters, primers etc) - must be located in the base.dir
artificial.fq <- "JGI.clean2.fa"
#
# path to the shell script that converts bam to wig
bam2wigPath <- "/home/GLBRCORG/omoskvin/run/bam2wig.sh"
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Getting input values:
# defaults - End
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
