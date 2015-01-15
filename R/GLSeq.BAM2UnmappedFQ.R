###################################################
# Extracts unmapped reads from a BAM file and converts them to FASTQ format
###################################################
#
##############
# objects required: 
#     dest.dir
#     picardToolsPath
#     refFASTAname
##############
#
source("GLSeq.Util.R")
#
setwd(dest.dir)
bamfiles2 <- dir(pattern="sorted.bam.bai") # names of the indices
bamfiles <- substr(bamfiles2, 1, nchar(bamfiles2)-4)
#
for (bamfile in 1:length(bamfiles)) {
currentBAM <- bamfiles[bamfile]
current.filebase <- substr(currentBAM, 1, nchar(currentBAM)-4) # introduced later on and used in system command #2 only
currentUnmappedSAM <- paste(substr(currentBAM, 1, nchar(currentBAM)-3), "unmapped.sam", sep="")
currentUnmappedSAM.comm <- paste("samtools view -f4 ", currentBAM, " > ", currentUnmappedSAM, sep="") # system command #1 
#
# adding headers to satisfy subsequent picard processing:
betterUnmappedSAM <- paste(substr(currentUnmappedSAM, 1, nchar(currentUnmappedSAM) - 4), "HEAD.sam", sep="")
headerAdd.comm <- paste("java -Xmx2g -jar ",picardToolsPath, "AddOrReplaceReadGroups.jar", " I=", currentUnmappedSAM, " O=", betterUnmappedSAM, " SO=coordinate LB=", refFASTAname, " PL=ILLUMINA PU=unknown SM=", current.filebase, " VALIDATION_STRINGENCY=LENIENT", sep="") # System command #2
#
currentUnmappedFQ <- paste(substr(currentBAM, 1, nchar(currentBAM)-3), "unmapped.fq", sep="")
unmappedFQ.comm <- paste("java -Xmx2g -jar ",picardToolsPath, "SamToFastq.jar INPUT=", betterUnmappedSAM, " FASTQ=", currentUnmappedFQ, sep="") # System command #3
#
housekeeping.comm <- paste("rm", currentUnmappedSAM, "&& mv", betterUnmappedSAM, currentUnmappedSAM) # System command #4
comm.bamfile <- paste(currentUnmappedSAM.comm, "&&", headerAdd.comm, "&&",  unmappedFQ.comm, "&&", housekeeping.comm, "&")
system(comm.bamfile)
}
##############
	
