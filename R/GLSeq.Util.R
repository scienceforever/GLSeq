####################################################
# Utilities for GLSeq package
####################################################
#
#
####################
####################
# checking if the directory 
# names have trailing slashes
# and adding them if needed
####################
trailDirCheck <- function(dir.nm) {
dir.nm.endchar <- substr(dir.nm, nchar(dir.nm), nchar(dir.nm))
if (dir.nm.endchar != '/') dir.nm <- paste(dir.nm, '/', sep="")
dir.nm
}
#@@@@@@@@@@@@@@
#
###################
###################
# Collecting pathes to all the files that match a name pattern,
# relative to the current working directory (or any supplied as "wd"),
# and returning this as an R vector object
fileNamesCollect <- function(fNamePattern, wd="./") {
# single-quoted name pattern (to avoid "find: paths must precede expression" shell error):
initial.dir <- trailDirCheck(getwd())
wd <- trailDirCheck(wd)
setwd(wd) # this still may be a relative path like "./"!!!
wd <- trailDirCheck(getwd()) # now we have an absolute path
fNamePattern.squot <- paste("'",  fNamePattern, "'", sep="")
syscomm <- paste("find . -name", fNamePattern.squot)
out <- system(syscomm, intern=TRUE)
# removing of the "./" in the beginning of the output before adding absolute path in a general way:
out <- substr(out, 3,nchar(out))
# adding absolute path: 
out <- paste(wd, out, sep="")
# back to starting directory:
setwd(initial.dir)
out
}
######
# usage example: 
# source("~/run/GLSeq.Util.R")
# fileNamesCollect("Y128*")
######
#@@@@@@@@@@@@@@
#
###################
###################
# pooling reads from several libraries for transcriptome assembly, paired reads
readPoolPE <- function(poolName) {
left <- dir()[grep("1.fq", dir())]
right <- dir()[grep("2.fq", dir())]
left<- paste(left, collapse=" ")
right <- paste(right, collapse=" ")
leftPoolName <- paste(poolName, "_leftpool.fq", sep="")
rightPoolName <- paste(poolName, "_rightpool.fq", sep="")
bigLibrary <- paste("cat ", left, " >", leftPoolName, "& cat ", right, " >", rightPoolName,"&")
system(bigLibrary)
}
#@@@@@@@@@@@@@@
#
###################
###################
# Extracting file name from the full path
path2file <- function(fullPath) {
bricks <-  unlist(strsplit(fullPath, split="/"))
fName <- bricks[length(bricks)]
fName	
}
#@@@@@@@@@@@@@@
#
###################
###################
# 
# Converting all bam files within a directory to wig
bam2wig.batch <- function(path2util) {
# path2util is the path to the shell script that takes input and output file names
allfiles <- dir()
bamfiles <- allfiles[grep("bam$", allfiles, perl=TRUE)]
for (i in 1:length(bamfiles)) {
wigfile.i <- paste(substr(bamfiles[i], 1, nchar(bamfiles[i])-3), "wig", sep="")	
system(paste(path2util, bamfiles[i], wigfile.i, "&"))
}
}
#@@@@@@@@@@@@@
#
##################
##################
# updating destination and run directories
# for result cllections and results summarization
# with earlier-generated results
dirUpdate <- function(newRunDir, newDestDir, env = parent.frame()) {
env$base.dir <- trailDirCheck(newRunDir)
env$dest.dir <- trailDirCheck(newDestDir)
env$destDirBam <- paste(dest.dir, text.add, ".viz/", sep="" )
env$destDirCount <-  paste(dest.dir, text.add, ".counts/", sep="")
env$destDirLog <-  paste(dest.dir, text.add, ".stat/", sep="") 
} 
#@@@@@@@@@@@@@
#
##################
##################
# generating BigWIG files
wig2BigWig.DIR <- function(exePath="~/run/wigToBigWig", chr.sizes="./chrom.sizes") {
wigFiles <- dir()
wigFiles <- wigFiles[grep("\\.wig", wigFiles)]
for (i in 1:length(wigFiles)) {
currentWig <- wigFiles[i]
bwFile <- paste(substr(currentWig, 1, nchar(currentWig)-4), "bw", sep=".")
convertComm <- paste(exePath, currentWig, chr.sizes, bwFile, "&")
system(convertComm)
} }
#@@@@@@@@@@@@@
#
##################
##################
# Reporting the vectors of length 1 in the environment and their values
# (practical purpose: fast extracting of parameter information for the old runs in GLSeq from the rda files)
# maxNchar is the maximum number of characters we want to see in the value of an individual object 
# (objects with longer character values will be ignored, to avoid flooding the results with system command stacks etc)
show1Vec <- function(maxNchar=40) {
fullList <- ls(envir = .GlobalEnv)
select.IND <- rep(FALSE, length(fullList))
for (i in 1:length(fullList)) {
if (class(get(fullList[i])) != "function" & length(get(fullList[i])) ==1) { 
select.IND[i] <- TRUE } }
fullList <- fullList[select.IND]
fullList <- fullList[as.logical(nchar(sapply(fullList, get)) < maxNchar)]
droplist <- c("fullList", "i", "j", "out", "showEnv", "droplist") # handling second and subsequen runs
fullList <- fullList[!(fullList %in% droplist)]
for (j in 1:length(fullList)) {
	currentPiece <- c(fullList[j], get(fullList[j]))
        if (j == 1) out <- currentPiece
        if (j != 1) { out <- rbind(out, currentPiece) 
              if (i == 2) rownames(out)[1] <- 1
              rownames(out)[j] <- j } } 
out <- as.data.frame(out)
rownames(out) <- out[,1]
out <- out[,-1,drop=FALSE]
out }
#@@@@@@@@@@@@@


















