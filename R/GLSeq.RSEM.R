#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# RSEM quantification of the expression values
#########################################################
#
# Number of computation streams is set to 2 for the time being;
# with expected arrival of new hardware, a more general solution will be encoded 
# to take advantage of more than 16 cores available 
#
source("GLSeq.Util.R")
#
###################
# Generating RSEM option set:
###################
#
# Number of cores for every library:
RSEM.cores <- paste(' -p ', nCores,sep="") 
#
# Quality scores format?
# used to be (until revision 79):
# RSEM.qscores <- paste(' --', qScores, '-quals ', sep="")
# then, RSEM was updated:
### "Sept 25, 2013   RSEM v1.2.7 has a minor update. One line is added to the 'WHAT_IS_NEW' file to reflect a change made but forgotten to put in to 'WHAT_IS_NEW'. 
### The line added is "Renamed '--phred33-quals', '--phred64-quals', and '--solexa-quals' in 'rsem-calculate-expression' to '--bowtie-phred33-quals', '--bowtie-phred64-quals'..."
### 
# RSEM.qscores <- paste(' --bowtie-', qScores, '-quals ', sep="")
RSEM.qscores <- paste(' --', qScores, '-quals ', sep="")
# (reverted back in later versions)
#
# Paired-end sequencing? 
RSEM.pairedend <- NULL
if (paired.end) RSEM.pairedend <- ' --paired-end ' 
#
# output genomic BAM files?
RSEM.genobam <- ""
if (genobam) RSEM.genobam <-  ' --output-genome-bam ' 
#
# compute confidence intervals?
RSEM.confint <- ""
if (compConf) RSEM.confint <- " --calc-ci "
#
# strand-specifc library? which reference strand to align to? # --forward-prob 0.5 is the RSEM default (non-strand-specific libraries)
RSEM.libstrand <- ""
if (!(is.null(libstrand))) {
if (libstrand == "F") RSEM.libstrand <- " --forward-prob 1 "
if (libstrand == "R") RSEM.libstrand <- " --forward-prob 0 " }
# 
# advanced options:
RSEM.fragMaxLength <- ""
if (paired.end) RSEM.fragMaxLength <- paste(" --fragment-length-max", fragMaxLength)
RSEM.ciMem <- ""
if (compConf) RSEM.ciMem <- paste(" --ci-memory", ciMem)
options.RSEM <- paste(RSEM.cores, RSEM.qscores, RSEM.pairedend, RSEM.genobam, RSEM.confint, RSEM.libstrand, RSEM.fragMaxLength, RSEM.ciMem)
#@@@@@@@@@@@
# RSEM Options End
#@@@@@@@@@@@
#
# creating stacks of RSEM commands for each range of data files: 
# copy genome indices to the destimation dir: 
ref.dir <- paste(base.dir, rGenome, sep="")
indCopy <- paste("cd ", ref.dir, " && cp * ", dest.dir, sep="")
system(indCopy)
setwd(dest.dir)
# 
# stacks of shell commands:
#
comm.stack.pool <- NULL
for (zz in 1:nStreams) {
# shell commands for zz'th stream:
for (j in rangelist[[zz]]) {
# set of files (or file for SE library) for the current library with paths: 
if (ncol(fqfiles.table) == 1) fqfiles.paths <- paste(dest.dir, fqfiles.table[j,1], sep="") # single-end library
if (ncol(fqfiles.table) == 2) fqfiles.paths <- paste(paste(dest.dir, fqfiles.table[j,1], sep=""), paste(dest.dir, fqfiles.table[j,2], sep="")) # paired-end library
this.library <- substr(fqfiles.table[j,1], 1,libNchar)
this.resName <- paste(this.library, text.add, sep=".")
this.command <- paste(" date   >>  runLogFile1 ", " && ", "rsem-calculate-expression ", options.RSEM,  fqfiles.paths, rGenome, this.resName,  " 1>> ", runLogFile1, " 2>>", runErrFile1)
currentStackName <- paste("comm.stack", zz, sep="")
if (j==rangelist[[zz]][1]) assign(currentStackName,"\n")
assign(currentStackName, paste(get(currentStackName),  this.command,  " && ", sep="")) 
} # for j
# adding stack run completeness indicators: 
if (zz ==1) fileCompletenessID <- paste(text.add, ".completeExpression", sep="")
# adding a timestamp of finishing all the jobs in a job stack and creating the "jobs completed" marker file:
assign(currentStackName,paste(get(currentStackName), " date ", " && ", paste("echo  >  ", fileCompletenessID, ".", zz, sep=""), " & "))
comm.stack.pool <- paste(comm.stack.pool, get(currentStackName))
#
} # for zz
#
# collection of the results: 
collLog <- paste(destDirLog, text.add, ".ResultsCollectLog.txt", sep="")
collerr <- paste(destDirLog, text.add, ".ResultsCollectErrors.txt", sep="")
collResults <- paste("cd ", base.dir, " && ", "Rscript GLSeqResultsCollect.R ", text.add, base.dir, dest.dir, " 1>> ", collLog, " 2>> ", collerr, " &", sep="") 
# 
# pool of all the system commands (versions for +/- compute expression and +/- collect results): 
#
comm.stack.pool <- paste(comm.stack.pool, collResults, "\n", sep=" ") 
if (resCollect == "nocollect") collResults <- "\n"
if (exprRun == "noexprcalc") comm.stack.pool <- paste(collResults, "\n", sep=" ") 
#
# Making a record of the actual commands sent to the system by the script: 
commandFile <- paste(destDirLog, text.add, ".commandpool.txt", sep="")
write.table(comm.stack.pool, file=commandFile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
################################################################################
