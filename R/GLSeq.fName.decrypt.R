#######################################################
# Adding meaningful experimental condition IDs 
# (those used by GLSeq.remap.R)
# to the beginning of the library-centric file names 
# in a given directory
# June 18, 2016; highly simplified GLSeq.remap.R
#######################################################
#
args <- commandArgs(trailingOnly = TRUE)
#
# directory with the data: 
theDir <- as.character(args[1])
#
# remapping CSV file (library names to condition IDs)
mapCSV <- as.character(args[2])
# 
setwd(theDir)
avail.files <- dir()
theMap <- read.table(mapCSV,sep="\t",header=FALSE,row.names=1,as.is=TRUE)
#
# determining the number of unique chracters in library IDs:
nUnique <- nchar(rownames(theMap)[1])
#
for (i in 1:length(avail.files)) {
currentFname <- avail.files[i]
libUniqueString <- substr(currentFname,1,nUnique)
conditionID <- theMap[libUniqueString,1]
newFname <- paste(conditionID, currentFname, sep=".")
if (!(is.na(newFname))) system(paste("mv", currentFname, newFname))
}
#########################################################
# Example usage: 
# Rscript GLSeq.fName.decrypt.R  /data/viz/  /maps/yeast_remap.csv 



