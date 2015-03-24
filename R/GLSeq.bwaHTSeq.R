#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# Trimmomatic-BWA-HTSeq quantification of the expression values
# (paired end libraries); based on "Dana's script" + Irene input + Oleg's observations on the data
# (results of the latter: default headcrop is increased to 12 and Illumina artifacts removal is introduced)
# September 2013
#########################################################
#
#
source("GLSeq.Util.R")
# if (!(paired.end)) stop('BWA pipeline is currently supported for paired-end libraries only \n')
setwd(dest.dir)
# copy genome indices to the destimation dir: 
ref.dir <- paste(base.dir, rGenome, sep="")
indCopy <- paste("cd ", ref.dir, " && cp * ", dest.dir, sep="")
system(indCopy)
#
comm.stack.pool <- NULL # 
#
for (zz in 1:nStreams) {
# assembly and runing the system command, one library at a time:
for (i in rangelist[[zz]]) {
###################
# Alignment with SAM output
###################
# names of current fastq files:
fq.left <- fqfiles.table[i,1]
if (paired.end) fq.right <- fqfiles.table[i,2]
this.library <- substr(fqfiles.table[i,1], 1,libNchar)
this.resName <- paste(this.library, text.add, sep=".")
unsorted.sam <- paste(this.resName, "unsorted.sam", sep=".")
#
# names of the expected sai files:
sainame.left <- paste(fq.left,"sai",sep=".")
if (paired.end) sainame.right <- paste(fq.right,"sai",sep=".")
# alignment commands:
aln.left <- paste(bwaPath, "aln", refFASTAname, fq.left, ">", sainame.left) # System command #1
if (paired.end) aln.right <- paste(bwaPath, "aln", refFASTAname, fq.right, ">", sainame.right) # System command #2
# creating SAM file: 
if (paired.end) sam.create <- paste(bwaPath, "sampe", refFASTAname, sainame.left, sainame.right, fq.left, fq.right, ">>", unsorted.sam) # System command #3
if (!(paired.end)) sam.create <- paste(bwaPath, "samse", refFASTAname, sainame.left, fq.left, ">>", unsorted.sam)
#
###################
# SAM file cleanup
###################
# the executable:
cleanSAM <- paste("java -Xmx2g -jar ",picardToolsPath, "CleanSam.jar", sep="")
# name of the cleaned SAM file:
cleaned.sam <- paste(this.resName, "cleaned.sam", sep=".")
# 
# SAM cleanup system command:
cleansam.comm <- paste(cleanSAM, " I=", unsorted.sam, " O=", cleaned.sam, sep="") # System command #4
#
###################
# Adding RG Header + sorting
###################
# the executable: 
headersortSAM <- paste("java -Xmx2g -jar ",picardToolsPath, "AddOrReplaceReadGroups.jar", sep="")
# name of the processed (final) SAM file:
final.sam <- paste(this.resName, "final.sam", sep=".")
finalsam.comm <- paste(headersortSAM, " I=", cleaned.sam, " O=", final.sam, " SO=coordinate LB=", refFASTAname, " PL=ILLUMINA PU=unknown SM=", this.resName, " VALIDATION_STRINGENCY=LENIENT", sep="") # System command #5
#
###################
# SAM => BAM file with index 
###################
sorted.arg <- paste(this.resName, "sorted", sep=".")
ref.index <- paste(refFASTAname, "fai", sep=".") 
sorted.bam <- paste(this.resName, "sorted.bam", sep=".") 
bam.create <- paste("samtools view -uS -t ", ref.index, final.sam, " | samtools sort - ", sorted.arg) # System command #6
bam.index <- paste("samtools index", sorted.bam) # System command #7
#
###################
# Converting the final bam (coordinate-sorted)
# to name-sorted  sam:  
###################
#
countable.sam <- paste(this.resName, "countable.sam", sep=".")
paired.arg <- paste(this.resName, "paired", sep=".")
paired.bam <- paste(this.resName, "paired.bam", sep=".")
# let's avoid pipe here to save RAM
countable.comm <- paste("samtools sort -n", sorted.bam, paired.arg, "&&", "samtools view -h", paired.bam, ">", countable.sam) # System command #8
#
###################
# Expresion counting (HTSeq) 
###################
# name-sorted ("cleaned") sam is used
# to support mate pair information 
countfile <- paste(this.resName, "counts", sep=".") 
# initializing strandness option with a non-strand-specific value:
countOpt <-  "--stranded=no"
# overwriting countOpt with information on strand-specificity:
if (libstrand == "R") countOpt <-  "--stranded=reverse"
if (libstrand == "F") countOpt <-  "--stranded=yes"
# selecting ID attribute for counts reporting: 
countOpt <- paste(countOpt, " --idattr=", idAttr, sep="")
count.comm <- paste("python -m HTSeq.scripts.count", countOpt, countable.sam, refGFFname, " > ", countfile) # System command #9
#

###################
# Expression counting (FeatureCount) 
###################




###################
# Housekeeping 
###################
if (paired.end) spaceCleanup <- paste("rm", sainame.left, "&& rm", sainame.right, "&& rm", unsorted.sam, "&& rm", cleaned.sam, "&& rm", final.sam, "&& rm", paired.bam) # System command #10
if (!(paired.end)) spaceCleanup <- paste("rm", sainame.left, "&& rm", unsorted.sam, "&& rm", cleaned.sam, "&& rm", final.sam, "&& rm", paired.bam) # System command #10
#
# current command: 
if (paired.end) comm.i <- paste(aln.left, "&&", aln.right, "&&", sam.create, "&&", cleansam.comm, "&&", finalsam.comm, "&&", bam.create, "&&", bam.index, "&&", countable.comm, "&&", count.comm, "&&", spaceCleanup)
if (!(paired.end)) comm.i <- paste(aln.left, "&&", sam.create, "&&", cleansam.comm, "&&", finalsam.comm, "&&", bam.create, "&&", bam.index, "&&", countable.comm, "&&", count.comm, "&&", spaceCleanup)
# for the very first assembly in the stack: 
if (i == rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool, " date && ", comm.i)
# for subsequent assemblies of every stack: 
if (i != rangelist[[zz]][1])  comm.stack.pool <- paste(comm.stack.pool, " && date && ", comm.i)
# system(comm.i)
} # for i 
if (zz ==1) fileCompletenessID <- paste(text.add, ".completeExpression", sep="")
comm.stack.pool <- paste(comm.stack.pool,  " && echo  >  ", fileCompletenessID, ".", zz, " & ", sep="")
} # for zz
#
# collection of the results: 
collLog <- paste(destDirLog, text.add, ".ResultsCollectLog.txt", sep="")
collerr <- paste(destDirLog, text.add, ".ResultsCollectErrors.txt", sep="")
collResults <- paste("cd ", base.dir, " && ", "Rscript GLSeqResultsCollect.R ", text.add, base.dir, dest.dir, " 1>> ", collLog, " 2>> ", collerr, " &", sep="") 
if (resCollect == "nocollect") collResults <- "\n"
# 
# pool of all the system commands (versions for +/- compute expression and +/- collect results): 
#
comm.stack.pool <- paste(comm.stack.pool, collResults, "\n", sep=" ") 
if (exprRun == "noexprcalc") comm.stack.pool <- paste(collResults, "\n", sep=" ") 
#
###########################################################################




