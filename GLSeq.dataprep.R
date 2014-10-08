#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Oleg Moskvin; info@scienceforever.com 
# April 2013 
#########################################################
# 
# Preparing the raw data
#
#########################################################
#
# ID of the computation run (text.add object) should be supplied to this script directly (by GLSeq.top.R) as the first and only argument
args <- commandArgs(trailingOnly = TRUE)
text.add <- as.character(args[1])
source("GLSeq.Util.R")
#
# loading variables for the current run (from the current directory) 
currentRun.dataFile <- paste("GLSeq.vars.", text.add, ".rda", sep="")
load(currentRun.dataFile) # loads all parameters for the run with the runID indicated in the text.add object
#
################################################################################
# routine for compressed (".gz") fastq files, paired-end sequencing and concatenated (1-st and 2-nd reads pooled) FASTQ files,  
# i.e. current situation in GLBRC/JGI for both E.coli and yeast data (April 2013);
# with appearence of new data types, respective data preparation blocks will be added as needed
################################################################################
# 
### assembling Trimmomatic system commands: 
if (paired.end) {
trimAssemble <- function(leftDirtyName, rightDirtyName, trimPath, qScore, headcrop=12, artifactsFile) {
dashqScore <- paste("-", qScore, sep="")
logName <- paste(leftDirtyName, "pairedtrim.log", sep=".")
pairedTrimmed.1 <- paste("p", leftDirtyName, sep=".")
unpairedTrimmed.1 <- paste("u",  leftDirtyName, sep=".")
pairedTrimmed.2 <- paste("p",  rightDirtyName, sep=".")
unpairedTrimmed.2 <- paste("u", rightDirtyName, sep=".") 
# trimParam <- paste("HEADCROP:", headcrop, " SLIDINGWINDOW:3:30 MINLEN:36 CROP:", crop, sep="")
trimParam <- paste("ILLUMINACLIP:", artifactsFile, ":2:30:10", " HEADCROP:", headcrop, " SLIDINGWINDOW:3:30 MINLEN:", trimMin, sep="")
tcomm <- paste("java -jar", trimPath, "PE -threads 20", dashqScore, "-trimlog", logName, leftDirtyName, rightDirtyName, pairedTrimmed.1, unpairedTrimmed.1, pairedTrimmed.2, unpairedTrimmed.2, trimParam)
tcomm }}
#
# with single-end libraries, Trimmomatic commands are shorter; let's keep 'p' prefix for the trimmed .fq files for consistency
if (!(paired.end)) {
trimAssemble <- function(leftDirtyName, trimPath, qScore, headcrop=12, artifactsFile) {
dashqScore <- paste("-", qScore, sep="")
logName <- paste(leftDirtyName, "pairedtrim.log", sep=".")
pairedTrimmed.1 <- paste("p", leftDirtyName, "fq", sep=".") # '.fq' is added here for single end libraries (it is added at the splitting stage for the paired-end)
# trimParam <- paste("HEADCROP:", headcrop, " SLIDINGWINDOW:3:30 MINLEN:36 CROP:", crop, sep="")
trimParam <- paste("ILLUMINACLIP:", artifactsFile, ":2:30:10", " HEADCROP:", headcrop, " SLIDINGWINDOW:3:30 MINLEN:", trimMin, sep="")
tcomm <- paste("java -jar", trimPath, "SE -threads 20", dashqScore, "-trimlog", logName, leftDirtyName,  pairedTrimmed.1, trimParam)
tcomm }}
###
setwd(dest.dir)
# files indicating readiness of the libraries (i.e. split FQ files):
files2watch.dataprep <- NULL
fqFiles.gz <- NULL
# here we have 2 ways of getting the list of compressed FASTQ files:
# 1) using supplied list of raw files (planned to be the standard way): 
if (!(is.null(libList))) fqFiles.gz <- libList
# 2) all the raw files are in the same directory 
# and the list of files is not explicitly supplied (the situation we had before intoducing a route of communication with CBDB): 
if (length(unique(raw.dir)) == 1 & is.null(libList))  fqFiles.gz <- dir(raw.dir[1])[grep(".gz", dir(raw.dir[1]))]
# if fqFiles.gz vector is still empty for some reason:
if (is.null(fqFiles.gz)) stop("Please check the list of raw FASTQ files and the contents of the raw directory \n")
# 
# names of the respective uncompressed FASTQ files (to be generated): 
fqFiles <- substr(fqFiles.gz, 1, nchar(fqFiles.gz)-3)
# 
# Ranges for data preparation: 
rangelist.Dataprep <- chunk(1:length(fqFiles), nStreamsDataPrep)
# data prep ranges (for fqFiles and fqFiles.gz) are in rangelist.Dataprep[1:nSreamsDataPrep]
for (zz in 1:nStreamsDataPrep) {
comm.pool <- "date " # separate command pool for every stream
for (j in rangelist.Dataprep[[zz]]) {
# fixing the "cannot create regular file `.....JGI.clean2.fa' cp: : File exists" error that started in Aug 2014 with all the same code 
# (preventing second to last raw files from being copied)  
if (j == 1) {
if (length(unique(raw.dir)) == 1) copy.comm <- paste("cp ", base.dir, artificial.fq, " ", dest.dir, " && cp ", raw.dir, fqFiles.gz[j], " ", dest.dir, sep="")
if (length(unique(raw.dir)) > 1) copy.comm <- paste("cp ", base.dir, artificial.fq, " ", dest.dir, " && cp ", raw.dir[j], fqFiles.gz[j], " ", dest.dir, sep="") }
if (j != 1) {
if (length(unique(raw.dir)) == 1) copy.comm <- paste("cp ", raw.dir, fqFiles.gz[j], " ", dest.dir, sep="")
if (length(unique(raw.dir)) > 1) copy.comm <- paste("cp ", raw.dir[j], fqFiles.gz[j], " ", dest.dir, sep="") }
#
gunzip.comm <- paste("gunzip ", fqFiles.gz[j], sep="")
# handling unzipped (but not split) files:
if (unzipped) gunzip.comm <- ""
# fqFile.base <-  substr(fqFiles[j], 1, nchar(fqFiles[j])-6)
# since files may be named either .fq or .fastq, we'll leave this extension alone for safety and add our extention on top of that:
fqFile.base <- fqFiles # May 29, 2013 
first.read.filename <- paste(fqFile.base[j], ".1.fq", sep="")
second.read.filename <- paste(fqFile.base[j], ".2.fq", sep="")
split.comm <- paste("cat ", dest.dir, fqFiles[j], " | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[1234]/; @i = 0 if @i == 8' > ", first.read.filename, " && cat ",   dest.dir, fqFiles[j], " | ruby -ne 'BEGIN{@i=0} ; @i+=1; puts $_  if @i.to_s =~ /[5678]/; @i = 0 if @i == 8' > ", second.read.filename, sep="")
# (the ruby code incorporated into the line above was adapted from SeqAnsweres forum (seqanswers.com))
#
# QC report for the untrimmed reads (always run): 
if (!(readTrim) & paired.end) preQC <- paste(fastqcPath, first.read.filename, second.read.filename)
if (!(readTrim) & !(paired.end)) preQC <- paste(fastqcPath, fqFile.base)
#
# trimming-specific commands: 
if (readTrim) {
leftDirtyFname <-  paste("dirty.", fqFile.base[j], ".1.fq", sep="")
rightDirtyFname <-  paste("dirty.", fqFile.base[j], ".2.fq", sep="")
pairedTrimmed.1 <-  paste("p.", fqFile.base[j], ".1.fq", sep="")
pairedTrimmed.2 <- paste("p.", fqFile.base[j], ".2.fq", sep="")
# 
# special case of SE libraries:
unpaired.fq <- paste(fqFile.base[j], ".fq", sep="") #
SE.dirtyFname <-  paste("dirty.", fqFile.base[j], ".fq", sep="")
SE.trimmedFname <-  paste("p.", fqFile.base[j], ".fq", sep="")	
#
if (paired.end) trimCommand <- trimAssemble(first.read.filename, second.read.filename, trimPath, qScores, trimHead, artificial.fq)
if (!(paired.end)) trimCommand <- trimAssemble(fqFile.base[j], trimPath, qScores, trimHead, artificial.fq)
if (paired.end) fileShuffle <- paste("mv", first.read.filename, leftDirtyFname, " && ", "mv", second.read.filename, rightDirtyFname, " && ", "mv", pairedTrimmed.1, first.read.filename, " && ", "mv", pairedTrimmed.2, second.read.filename)
if (!(paired.end)) fileShuffle <- paste("mv", fqFile.base[j], SE.dirtyFname, " && ",  "mv", SE.trimmedFname, unpaired.fq)
if (paired.end) preQC <- paste(fastqcPath, leftDirtyFname, rightDirtyFname)
if (!(paired.end)) preQC <- paste(fastqcPath, fqFile.base[j])
# (the above results in trimmed fq files with the names referenced elsewhere for further processing)
# with the name shiffling, command for QC after trimming will look the same:
if (paired.end) postQC <- paste(fastqcPath, first.read.filename, second.read.filename)
if (!(paired.end)) postQC <- paste(fastqcPath, unpaired.fq) # unpaired.fq is the FINAL TRIMMED file in the case of single-end library!
} # if readTrim
#
rm.comm <- paste("rm ", fqFiles[j])
readyfile.j <- paste(fqFiles[j], ".ready", sep="")
files2watch.dataprep <- c(files2watch.dataprep, readyfile.j)
ready.comm <- paste("echo > ", readyfile.j, sep="")
#
# final command stacks for no-trimming and trimming versions: 
if (!(readTrim) & paired.end) comm.pool.j <- paste(copy.comm, " && ", gunzip.comm, " && ", split.comm, " && ",  preQC, " && ", rm.comm, " && ", ready.comm)
if (!(readTrim) & !(paired.end)) comm.pool.j <- paste(copy.comm, " && ", gunzip.comm, " && ",  preQC, " && ", rm.comm, " && ", ready.comm)
if (readTrim & paired.end) comm.pool.j <- paste(copy.comm, " && ", gunzip.comm, " && ", split.comm, " && ",  trimCommand, " && ", fileShuffle, " && ", preQC, " && ", postQC, " && ",  ready.comm)
if (readTrim & !(paired.end)) comm.pool.j <- paste(copy.comm, " && ", gunzip.comm, " && ",  trimCommand, " && ", fileShuffle, " && ", preQC, " && ", postQC, " && ", ready.comm)
if (j != rangelist.Dataprep[[zz]][length(rangelist.Dataprep[[zz]])]) comm.pool <- paste(comm.pool, " && ", comm.pool.j)
if (j == rangelist.Dataprep[[zz]][length(rangelist.Dataprep[[zz]])]) {
comm.pool <- paste(comm.pool, " && ", comm.pool.j, " &")	
system(comm.pool) }
} # for j
} # for zz
#
# I.e., the script is run from base.dir, raw data are taken from raw.dir, processed in dest.dir and the ready-to-go (split) FQ files are placed into dest.dir 
#
# Checking for completeness: 
dataReady <- FALSE 
while(!(dataReady)) {
Sys.sleep(20) 
nProcessed.data <- sum(files2watch.dataprep %in% dir(dest.dir))
dataReady <-  nProcessed.data  == length(files2watch.dataprep) 
# number of libraries already processed:
if (!(dataReady)) cat(timestamp(), " Data preparation is not completed yet", " \n", nProcessed.data, " libraries are prepared so far", "\n")
if (dataReady) cat(timestamp(), " Data preparation is completed!", " \n",  nProcessed.data, " libraries are ready to go", "\n") }
#
if (dataReady) {
dataReady.signal <- paste("echo > ", text.add, ".DataReady", sep="") 
system(dataReady.signal) }
#################
