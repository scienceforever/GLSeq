#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Oleg Moskvin; info@scienceforever.com 
# April 2013 
#########################################################
# 
# Collection of the transcript abundance results, generation of experiemnt-wide summaries, 
# generation of starnd-specific visualization files and file sorting / houskeeping  
#
#########################################################
#
# ID of the computation run (text.add object) should be supplied to this script directly (by GLSeq.top.R) as the first and only argument
args <- commandArgs(trailingOnly = TRUE)
text.add <- as.character(args[1])
newRunDir <- as.character(args[2])
newDestDir <- as.character(args[3])
forceStart <- as.logical(as.numeric(args[4]))
source("GLSeq.Util.R")
#
# loading variables for the current run (from the current directory) 
currentRun.dataFile <- paste("GLSeq.vars.", text.add, ".rda", sep="")
load(currentRun.dataFile) # loads all parameters for the run with the runID indicated in the text.add object
#
# updating the directories for post-hoc runs (Jan 16, 2014):
dirUpdate(newRunDir, newDestDir)
cat('the new destination directory is', dest.dir, "\n")
#
############################
# watching for the completeness of the computation results
############################
#
setwd(dest.dir)
# the last step is indexing the bam file -> appearence of a .bai file
bai.pull <- function(dest.dir, text.add) {
allfiles <- dir(dest.dir)
bai.files <- allfiles[grep("bam.bai", allfiles)]
# just in case, let's restrict the list to the files with the current "text.add" mark:
bai.files.thisRun <- bai.files[grep(text.add, bai.files)]
# those are mix of transcript- and genome-level bai files;
# restricting to genome-level files: 
bai.files.thisRun.genome <- bai.files.thisRun[grep("genome.sorted.bam.bai", bai.files.thisRun)]
# if genome bam file is not requested, then transcriptome bai (the only indexed bam's captured by bai.pull() )is being used instead:
if (!(genobam)) bai.files.thisRun.genome <- bai.files
bai.files.thisRun.genome
}
##
# expression computation script is expected to save 
# ("text.add").completeExpression.1 .. to .nStreams files
run.ready.check <- function(runID, nStreams) {
completed <- TRUE 
for(stream in 1:nStreams) {
file2seek <- paste(runID, "completeExpression", stream, sep=".") 
completed <- completed & file2seek %in% dir(dest.dir) }
completed 
} 
#
# checking for expression calculation completeness every 10 minutes:
completeExpression.boo <- run.ready.check(text.add, nStreams)
if (forceStart) completeExpression.boo <- TRUE
while(!(completeExpression.boo)) {
Sys.sleep(600) 
completeExpression.boo <- run.ready.check(text.add, nStreams) 
# number of libraries already processed:
if (qAlgor == "RSEM") nProcessed <- length(bai.pull(dest.dir, text.add))
if (qAlgor == "bwaHTSeq") nProcessed <- length(dir(dest.dir, pattern=paste(".",text.add,".counts",sep="")))
if (!(completeExpression.boo)) cat(timestamp(), " Expression computation is not completed yet", " \n", nProcessed, " libraries are processed so far", "\n")
if (completeExpression.boo) cat(timestamp(), " Expression computation is completed!", " \n",  nProcessed, " libraries are processed", "\n") }
#
############################
# If the files are ready, 
# start extracting expression values
# and preparing visualization sets
############################
# 
if (completeExpression.boo) {
if (qAlgor == "RSEM") {
# 
############################
# generating matrices of computation results
# for all the libraries in the run
############################
#
# name signature of the results files:
results.fileSig <- paste(text.add, "genes.results", sep=".")
# file names with library-centric results:
result.fnames <- dir(pattern=results.fileSig) 
# library name + text add:
result.names <- substr(result.fnames, 1, nchar(result.fnames)-14)  
# just library names, in the respective order:
lib.names <- substr(result.names, 1,libNchar)
# collecting expected count and other processed data:
counts <- NULL
counts_pme <- NULL
FPKM <- NULL
FPKM_pme <- NULL
FPKM_lower <- NULL
FPKM_upper <- NULL
TPM <- NULL
TPM_pme <- NULL
TPM_lower <- NULL
TPM_upper <- NULL
for (i in 1:length(result.fnames)) {
i.data <- read.table(result.fnames[i], header=TRUE, sep="\t", row.names=1, as.is=TRUE)
counts <- cbind(counts, i.data[,"expected_count"])
# counts_pme <- cbind(counts_pme, i.data[,"pme_expected_count"]) # field name changed in RSEM v.1.2.15 (Jun 16, 2014)
counts_pme <- cbind(counts_pme, i.data[,"posterior_mean_count"])
FPKM <- cbind(FPKM, i.data[,"FPKM"])
FPKM_pme <- cbind(FPKM_pme, i.data[,"pme_FPKM"])
FPKM_lower <- cbind(FPKM_lower, i.data[,"FPKM_ci_lower_bound"])
FPKM_upper <- cbind(FPKM_upper, i.data[,"FPKM_ci_upper_bound"])
TPM <- cbind(TPM, i.data[,"TPM"])
TPM_pme <-  cbind(TPM_pme, i.data[,"pme_TPM"])
TPM_lower <-  cbind(TPM_lower, i.data[,"TPM_ci_lower_bound"])
TPM_upper <- cbind(TPM_upper, i.data[,"TPM_ci_upper_bound"])
} # for i
colnames(counts) <- result.names 
colnames(counts_pme) <- result.names 
colnames(FPKM) <- result.names 
colnames(FPKM_pme) <- result.names 
colnames(FPKM_lower) <- result.names 
colnames(FPKM_upper) <- result.names 
colnames(TPM) <- result.names 
colnames(TPM_pme) <- result.names 
colnames(TPM_lower) <- result.names 
colnames(TPM_upper) <- result.names 
#
rownames(counts) <- rownames(i.data)
rownames(counts_pme) <- rownames(i.data)
rownames(FPKM) <- rownames(i.data) 
rownames(FPKM_pme) <- rownames(i.data) 
rownames(FPKM_lower) <- rownames(i.data) 
rownames(FPKM_upper) <- rownames(i.data)
rownames(TPM) <- rownames(i.data) 
rownames(TPM_pme) <- rownames(i.data)
rownames(TPM_lower) <- rownames(i.data) 
rownames(TPM_upper) <- rownames(i.data)
#
counts <- round(counts, 0)
counts_pme <- round(counts_pme, 0)
#
counts.fName <- paste(destDirCount, text.add, ".counts.csv", sep="")
counts_pme.fName <-  paste(destDirCount, text.add, ".counts_pme.csv", sep="")
FPKM.fName <-  paste(destDirCount, text.add, ".FPKM.csv", sep="")
FPKM_pme.fName <-  paste(destDirCount, text.add, ".FPKM_pme.csv", sep="")
FPKM_lower.fName <-  paste(destDirCount, text.add, ".FPKM_lower.csv", sep="")
FPKM_upper.fName <-  paste(destDirCount, text.add, ".FPKM_upper.csv", sep="")
TPM.fName <-  paste(destDirCount, text.add, ".TPM.csv", sep="")
TPM_pme.fName <-  paste(destDirCount, text.add, ".TPM_pme.csv", sep="")
TPM_lower.fName <-  paste(destDirCount, text.add, ".TPM_lower.csv", sep="")
TPM_upper.fName <-  paste(destDirCount, text.add, ".TPM_upper.csv", sep="")
#
write.csv(counts, file=counts.fName)
write.csv(counts_pme, file=counts_pme.fName)
write.csv(FPKM, file=FPKM.fName)
write.csv(FPKM_pme, file=FPKM_pme.fName)
write.csv(FPKM_lower, file=FPKM_lower.fName)
write.csv(FPKM_upper, file=FPKM_upper.fName)
write.csv(TPM, file=TPM.fName)
write.csv(TPM_pme, file=TPM_pme.fName)
write.csv(TPM_lower, file=TPM_lower.fName)
write.csv(TPM_upper, file=TPM_upper.fName)
#
#@@@@@@@@@@@@@@@@@@@@@@@
# generating matrices END
#@@@@@@@@@@@@@@@@@@@@@@@
#
#
############################
# collecting genome-level BAM and BAI files,
# generating wiggle files
# and moving all of them (.bam, .bai, .wig)
# to a chosen folder (destDirBam)
############################
# if genome file is not requested (genobam == FALSE), then all the files below will indeed represent reanscriptome visualization files
# because bai.files.thisRun.genome will actually contain transcriptome bam's (see bai.pull function for details) 
bai.files.thisRun.genome <- bai.pull(dest.dir, text.add)
# respective names of the bam and wig files: 
bam.files.thisRun.genome <- substr(bai.files.thisRun.genome, 1, nchar(bai.files.thisRun.genome)-4)
wig.files.thisRun.genome <- paste(substr(bam.files.thisRun.genome, 1, nchar(bam.files.thisRun.genome)-4),"wig", sep=".")
# 
# The case for generating wigle files without explicit separating forward- and reverse-strand reads
# (implies using of a standard gtf file without addition of reverse-strand features):
if(!(strandExtract)) { 
# Generating .wig files: 
for (bam in 1:length(bam.files.thisRun.genome)) {
current.bam <- bam.files.thisRun.genome[bam]
current.bai <- bai.files.thisRun.genome[bam]
# name of the wig file to generate: 
current.wig <- wig.files.thisRun.genome[bam]
# base name (before ".genome.sorted.bam"): 
current.base <- substr(current.bam, 1, nchar(current.bam)-18)
wig.command <- paste("rsem-bam2wig ", current.bam, current.wig, current.base, sep=" ")
system(wig.command) 
}
cat(timestamp(), " Wiggle files are generatd!", "\n") 
} # if not strandExtract
# strand extract version is common to RSEM and BWA-HTSeq (below)
} # if qAlgor == RSEM 
###############################
# Collection in the case of BWA-HTSeq route
###############################
if (qAlgor == "bwaHTSeq") {
setwd(dest.dir)
cfiles <- dir(pattern="counts")
countDirName <- paste(text.add, "counts", sep=".")
cfiles <- cfiles[cfiles != countDirName]
count.mtrx <- NULL
for (i in 1:length(cfiles)) {
data.i <- read.table(cfiles[i], row.names=1)
colnames(data.i) <- substr(cfiles[i],1,libNchar)
if (is.null(count.mtrx)) count.mtrx <- data.i
if (!(is.null(count.mtrx)) & i != 1) count.mtrx <- cbind(count.mtrx, data.i) } 
# collecting library-specific counts end
#
# generating file with report on all abnormal events during alignment:
droplines <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned" , "alignment_not_unique")
exceptionReport.boo <- rownames(count.mtrx) %in% droplines
exceptionReport <- count.mtrx[exceptionReport.boo, ]
exrep.fName <-  paste(destDirCount, text.add, ".exceptionReport.csv", sep="")
write.csv(exceptionReport, file=exrep.fName)
#
# saving raw counts to disk:
count.mtrx <- count.mtrx[!(exceptionReport.boo), ]
counts.fName <-  paste(destDirCount, text.add, ".counts.csv", sep="")
write.csv(count.mtrx, file=counts.fName)
#
# Generating normalized counts (FPKM) and saving the file to disk: 
# function to pull GeneID from the first record of a particular element of the long character string ("mess") in the column 9 of the gtf file: 
genePull <- function(gtf.gene.record) { 
geneID.record <- strsplit(gtf.gene.record, split=";")[[1]][1]
geneID <- substr(geneID.record, 9, nchar(geneID.record))
geneID
}
# More general solution - 
# Extracting an arbitrary IDs from the annotation column of a gtf/gff file:
idPull <- function(annotColumn, ID2look) {
annotColumn.split <- strsplit(annotColumn, split=";")
ID2look.positions <- lapply(annotColumn.split, grep, pattern=ID2look)	
IDvec <- rep(NA, length(annotColumn.split))
for (vv in 1:length(IDvec)) {
if (length(ID2look.positions[[vv]]) > 0) { 
	IDcopyNum <- length(ID2look.positions[[vv]]) # the LAST copy of the same ID is actullay used in the count computation! 
	IDvec[vv] <- annotColumn.split[[vv]][ID2look.positions[[vv]]][IDcopyNum]
	IDvec[vv] <- substr(IDvec[vv], nchar(ID2look)+2, nchar(IDvec[vv]))
}
} # for vv
IDvec # the length of the original annotColumn is preserved (records not containing the ID2look now contain NAs)
}
#
# loading GTF:
gtf <-  read.table(refGFFname, sep="\t", header=FALSE, as.is=TRUE)
mess <- gtf[,gtfFeatureColumn] # 9
#@ dim(mess) <- c(length(mess), 1)
#@ geneIDs.full <- apply(mess,1,genePull)
# geneIDs.full now may stand for a vector of arbitrary IDs extracted with idPull()
geneIDs.full <- idPull(mess, idAttr)
lengthData <- cbind(gtf[,3:5], geneIDs.full) # columns: 1) feature, 2) start, 3) end, 4) ID
lengthData <- lengthData[lengthData[,1] == "exon",]
# we may see situations when "exon" record does not have  desireable ID record in the gtfFeatureColumn:
lengthData <- lengthData[!(is.na(lengthData[,4])),]
lengthData <- cbind(lengthData, lengthData[,3] - lengthData[,2])
colnames(lengthData)[5] <- "length"
#
## considering multiple exons per gene: 
duplicated.IDs.summary <- table(lengthData[,4])[table(lengthData[,4]) > 1]
lengthData.singleExons <- lengthData[!(lengthData[,4] %in% names(duplicated.IDs.summary)),]
lengthData.multipleExons <- lengthData[lengthData[,4] %in% names(duplicated.IDs.summary),]
lengthData <- lengthData.singleExons # will be extended to summarized exon length below
for (multiexon in 1:length(names(duplicated.IDs.summary))) {
data.multiexon <- lengthData.multipleExons[lengthData.multipleExons[,4] == names(duplicated.IDs.summary)[multiexon],]
# recording the sum of all exon's lengths for a multi-exon gene in the first line of the length data for this gene: 
data.multiexon[1,5] <- sum(data.multiexon[,5])
lengthData <- rbind(lengthData, data.multiexon[1,])
} # for multiexon
## end of multiple exons per gene summarization
#
lengthData <- lengthData[lengthData[,4] %in% rownames(count.mtrx),] # avoiding trouble with occasional NAs in the column 4 of the lengthData
rownames(lengthData) <- lengthData[,4]
# Sorting / restricting the count matrix based on the non-redundant ID vector: 
lengthData <- lengthData[rownames(count.mtrx),]
# count.mtrx <- count.mtrx[rownames(lengthData),]
RPKM.mtrx <- count.mtrx 
for (normCol in 1:ncol(RPKM.mtrx)) {
mil.mappedR <- sum(RPKM.mtrx[,normCol]) / 1000000
kbase.transcr <- lengthData[,5] / 1000
RPKM.mtrx[,normCol] <- RPKM.mtrx[,normCol] / (kbase.transcr * mil.mappedR)
} # for normCol
RPKM.fName <-  paste(destDirCount, text.add, ".RPKM.csv", sep="")
write.csv(RPKM.mtrx, file=RPKM.fName)
#
# visualization files (the name is common with RSEM to unify the code for strandExtract version:
bam.files.thisRun.genome <- dir(pattern="*.bam")
bai.ind <- grep("bam.bai", bam.files.thisRun.genome)
bam.files.thisRun.genome <- bam.files.thisRun.genome[-bai.ind]
bai.files.thisRun.genome <- bam.files.thisRun.genome[bai.ind]
#
if(!(strandExtract)) { 
vizfiles.base <- substr(bam.files.thisRun.genome, 1, nchar(bam.files.thisRun.genome)-4)
wigfiles <- paste(vizfiles.base, "wig", sep=".")
bamMove <- paste("mv *.bam*", destDirBam)
# setVizDir <- paste("cd", destDirBam) # non-stable (files are referred to as "not found" even when the access of them should happen AFTER the move and directory change (controlled with "&&")!!!
for (ii in 1:length(bam.files.thisRun.genome)) {
# wigGen.ii <- paste("rsem-bam2wig",  bamfiles[ii], wigfiles[ii], vizfiles.base[ii])
wigGen.ii <- paste(bam2wigPath, bam.files.thisRun.genome[ii], wigfiles[ii])
if (ii == 1) wigGen <- wigGen.ii
if (ii != 1) wigGen <- paste(wigGen, "&", wigGen.ii)
} # for bamfiles
try(system(bamMove)) 
Sys.sleep(5)
setwd(destDirBam)
system(paste(wigGen, "&"))
Sys.sleep(2)
} # if not strand extract
# setwd(base.dir)
} # if qAlgor == "bwaHTSeq" 
#
# July 2014
#
# The case for explicit extraction of the F and R strand coverage from the bam file 
# (implies using either raw genome alignments or "double" gtf file with antisense features added): 
#  generated strand-specific "S" and "A" files are treated differently depending on the direction of the library! - see 'libstrand' checkpoints
if(strandExtract) { 
wig.files.thisRun.genome <- NULL
bam.files.thisRun.genome.split <- NULL
bai.files.thisRun.genome.split <- NULL
# names of "ready files" (visualization completeness indicators) for all libraries: 
allReadyFiles <- NULL
# 
# Visualization file preparation: 
for (bam in 1:length(bam.files.thisRun.genome)) {
current.bam <- bam.files.thisRun.genome[bam]
current.bai <- bai.files.thisRun.genome[bam]
# names of the bam files representing "antisense" pairs (in fact, FORWARD GENOMIC strand for  libstrand =="R" case)
# please note the difference between read alignment strandness (relative to actual transcript)
# and strandness recorded in the FLAG of the alignment file (relative to genome) 
A1.bamName <- paste(substr(current.bam, 1, nchar(current.bam)-4),"A1.bam", sep=".")
A2.bamName <- paste(substr(current.bam, 1, nchar(current.bam)-4),"A2.bam", sep=".")
A.bamName <- paste(substr(current.bam, 1, nchar(current.bam)-4),"A.bam", sep=".")
# names of the bam files representing "sense" pairs (in fact, REVERSE GENOMIC strand for  libstrand =="R" case)
# please note the difference between read alignment strandness (relative to actual transcript)
# and strandness recorded in the FLAG of the alignment file (relative to genome) 
S1.bamName <- paste(substr(current.bam, 1, nchar(current.bam)-4),"S1.bam", sep=".")
S2.bamName <- paste(substr(current.bam, 1, nchar(current.bam)-4),"S2.bam", sep=".")
S.bamName <- paste(substr(current.bam, 1, nchar(current.bam)-4),"S.bam", sep=".")
#
# system commands, depending on the direction of the library, F and R bam files: 
# 
# for libraries representing R strand of the transcript
if (libstrand =="R") {
F.bamName <-   paste(substr(A.bamName, 1, nchar(A.bamName)-6),"F.bam", sep=".") 
R.bamName <-   paste(substr(S.bamName, 1, nchar(S.bamName)-6),"R.bam", sep=".") 
# bai file names (not needed for generation but neded for moving around): 
F.baiName <- paste(F.bamName, "bai", sep=".")
R.baiName <- paste(R.bamName, "bai", sep=".")
# 
	if (paired.end) {
bamF.command <- paste("samtools view -bf 83 ", current.bam, " > ", A1.bamName, " && ", " samtools view -bf 163 ", current.bam, " > ", A2.bamName, " && ", " samtools merge ", F.bamName, " ", A1.bamName, " ", A2.bamName, " && ", "samtools index ", F.bamName)
bamR.command <- paste("samtools view -bf 99 ", current.bam, " > ", S1.bamName, " && ", " samtools view -bf 147 ", current.bam, " > ", S2.bamName, " && ", " samtools merge ", R.bamName, " ", S1.bamName, " ", S2.bamName, " && ", "samtools index ", R.bamName)
				}                              
	if (!(paired.end)) {
bamF.command <- paste("samtools view -bf 16 ", current.bam, " > ", F.bamName, " && ", "samtools index ", F.bamName)
bamR.command <- paste("samtools view -bF 16 ", current.bam, " > ", R.bamName, " && ", "samtools index ", R.bamName)
                                   }                              
} # if libstrand R
# 
# for libraries representing F strand of the transcript
if (libstrand =="F") {
F.bamName <-   paste(substr(S.bamName, 1, nchar(S.bamName)-6),"F.bam", sep=".") 
R.bamName <-   paste(substr(A.bamName, 1, nchar(A.bamName)-6),"R.bam", sep=".") 
F.baiName <- paste(F.bamName, "bai", sep=".")
R.baiName <- paste(R.bamName, "bai", sep=".")
#
	if (paired.end) {	
bamR.command <- paste("samtools view -bf 83 ", current.bam, " > ", A1.bamName, " && ", " samtools view -bf 163 ", current.bam, " > ", A2.bamName, " && ", " samtools merge ", R.bamName, " ", A1.bamName, " ", A2.bamName, " && ", "samtools index ", R.bamName)
bamF.command <- paste("samtools view -bf 99 ", current.bam, " > ", S1.bamName, " && ", " samtools view -bf 147 ", current.bam, " > ", S2.bamName, " && ", " samtools merge ", F.bamName, " ", S1.bamName, " ", S2.bamName, " && ", "samtools index ", F.bamName)
				}
	if (!(paired.end)) {
bamR.command <- paste("samtools view -bf 16 ", current.bam, " > ", R.bamName, " && ", "samtools index ", R.bamName)
bamF.command <- paste("samtools view -bF 16 ", current.bam, " > ", F.bamName, " && ", "samtools index ", F.bamName)
                                   }    
} # if libstrand F
# 
# names of the wig files to generate: 
F.wig <-  paste(substr(F.bamName, 1, nchar(F.bamName)-4),"wig", sep=".")    
R.wig <-  paste(substr(R.bamName, 1, nchar(R.bamName)-4),"wig", sep=".")
# 
# base name (before ".genome.sorted.X.bam"): 
current.base <- substr(current.bam, 1, nchar(current.bam)-20)
F.base <- paste(current.base, "F", sep=".")
R.base <- paste(current.base, "R", sep=".")
# generating of wig files differs depending on quantification algorithm / resulting bam file peculiarities:
if (qAlgor == "RSEM") {
wigF.command <- paste("rsem-bam2wig ", F.bamName, F.wig, F.base, sep=" ")
wigR.command <- paste("rsem-bam2wig ", R.bamName, R.wig, R.base, sep=" ") }
if (qAlgor == "bwaHTSeq") {
wigF.command <- paste(bam2wigPath, F.bamName, F.wig)
wigR.command <- paste(bam2wigPath, R.bamName, R.wig) }
#
# file indicating completeness of a library processing: 
libReadyFile <- paste(current.base, "ready.txt", sep=".")
# pool of all "Ready Files" (presence of all of them in the dest.dir indicates completeness of visualization file preparation for all the libraries in the processing batch: 
allReadyFiles <- c(allReadyFiles, libReadyFile)
# 
# pack of commands for each library (will be sent to background execution, taking 1 core per library):
libraryCommands <- paste(bamF.command, " && ", bamR.command, " && ", wigF.command, " && ", wigR.command, " && ", "echo  >  ", libReadyFile, " &")
wig.files.thisRun.genome <- c(wig.files.thisRun.genome, F.wig, R.wig)
bam.files.thisRun.genome.split <-  c(bam.files.thisRun.genome.split, F.bamName, R.bamName)
bai.files.thisRun.genome.split <- c(bai.files.thisRun.genome.split, F.baiName, R.baiName)
system(libraryCommands) 
} # for bam (starndExtract) 
} # if strandExtract
#
# names of bam and bai files ready for visualization: 
if(!(strandExtract)) vis.files.2move <- c(bam.files.thisRun.genome, bai.files.thisRun.genome, wig.files.thisRun.genome)
if(strandExtract) vis.files.2move <- c(bam.files.thisRun.genome, bai.files.thisRun.genome, wig.files.thisRun.genome, bam.files.thisRun.genome.split, bai.files.thisRun.genome.split)
# write.table(vis.files.2move, file="vis.files.2move.txt", sep="\n") # for debugging
#
# if (genobam) try(system("rm *.transcript* & rm *.isoforms*  &"))
# (leaving transcript-level visualization in place for now) 
# settling down before moving of the bam files (used as sources in the previous step): 
Sys.sleep(20) 
#
# Checking if the visualization files are ready: 
vizReady.boo <- FALSE
while(!(vizReady.boo)) {
Sys.sleep(120) 
# number of libraries processed at visualization step at the moment:
nProcessed.viz <- sum(allReadyFiles %in% dir(dest.dir))
vizReady.boo <-  nProcessed.viz  == length(allReadyFiles) 
if (!(vizReady.boo)) cat(timestamp(), " Visualization file preparation is not completed yet", " \n", nProcessed.viz, " libraries are finished so far", "\n")
if (vizReady.boo) cat(timestamp(), "  Visualization file preparation is completed!", " \n",  nProcessed.viz, " libraries are processed", "\n") }
#
# moving the visualization files to the desired directory: 
if (vizReady.boo) {
for (visfile in 1:length(vis.files.2move)) {
Sys.sleep(1) 
current.mvfile <- vis.files.2move[visfile]
mv.command <- paste("mv ", current.mvfile, " ", destDirBam, " & ")
try(system(mv.command)) } 
#
# final housekeeping:
rm1 <- "rm *.A1.bam &"
rm2 <- "rm *.A2.bam &"
rm3 <- "rm *.S1.bam &"
rm4 <- "rm *.S2.bam &"
rm5 <- "rm *.ready.txt &"
rm6 <- "rm *.genome.bam &"
rm6b <- "rm *.transcript.bam &"
if (dataPrepare == "nodataprep") rm7 <- "rm *.1.fq && rm *.2.fq &" # the files were ready and kept in a separate place BEFORE this run; an extra copy of them is not needed 
if (dataPrepare == "dataprep") rm7 <- "/n" # leave final fastq files (potentially expensive after pre-treatment and QC) alone for keeping and further use in alternative computations
rm8 <- "rm *.fastq.ready &"
rm9 <- "rm *.DataReady &"
rm10 <- "rm *.completeExpression.* &"
# wildcarded stat directories names for all libraries: 
statdirs.common <- paste("*", text.add, "stat", sep=".")
statmove <- paste("mv",  statdirs.common, destDirLog, "&", "mv *.genes.results", destDirLog, "&& ") # "&&" ensures that DB update will proceed AFTER completeness of the visualization / housekeeping step
# writing analysis results back to the database:
writebacklog <- paste(destDirLog, text.add, ".writebacdir()kLog.txt", sep="")
writeback <- paste("Rscript GLSeq.writeback.R ", text.add, " 1>", writebacklog, " 2>", writebacklog, " &", sep="")
if (updateFromDb == "noupdate") writeback <- "date"
housekeeping.comm <- paste(rm1, rm2, rm3, rm4, rm5, rm6, rm6b, rm7, rm8, rm9, rm10, statmove, writeback)
Sys.sleep(10) 
try(system(housekeeping.comm)) 
}
#@@@@@@@@@@@@@@@@@@@@@@@
# BAM collection / wig generation / moving end
#@@@@@@@@@@@@@@@@@@@@@@@
#
# BigWig generation (assumes presence of 'chrom.sizes' file in the reference directory):
exePath=paste(base.dir, "wigToBigWig", sep="")
chr.sizes="../chrom.sizes"
setwd(destDirBam)
wig2BigWig.DIR(exePath, chr.sizes)
} # if completeExpression.boo 
#############################

