#######################################################
# Remapping counts data (with library IDs) to experimental condition IDs 
# July 16, 2013; last updated July 27, 2016
#######################################################
#
args <- commandArgs(trailingOnly = TRUE)
source("GLSeq.Util.R")
#
# directory with count data: 
countDir <- as.character(args[1])
#
# remapping CSV file (library names to condition IDs)
mapCSV <- as.character(args[2])
# 
# choice of column reordering: 
reorderMe <- as.logical(as.numeric(as.character(args[3])))
# 
setwd(countDir)
try(system("mkdir remapped")) 
avail.files <- dir()[grep("csv", dir())]
theMap <- read.table(mapCSV,sep="\t",header=FALSE,check.names=FALSE, colClasses=c("character", "character"), row.names=1)
#
# determining the number of unique chracters in library IDs:
nUnique <- nchar(rownames(theMap)[1])
#
for (i in 1:length(avail.files)) {
currentFname <- avail.files[i]
current.obj <- read.csv(currentFname, row.names=1,  check.names = FALSE)
used.libs <- substr(colnames(current.obj), 1, nUnique)
if (sum(used.libs %in% rownames(theMap)) != length(used.libs)) stop('Not all the actually used libraries are found in the supplied map file!')
new.colnames <- theMap[used.libs,1]
# cbind(colnames(current.obj), new.colnames) # sanity check => passed 
colnames(current.obj) <- new.colnames
#
# reordering the columns for downstream convenience, if selected: 
if (reorderMe) {
cols.ordered <-  colnames(current.obj)[order(colnames(current.obj))]
current.obj <- current.obj[,cols.ordered] }
baseFname <- substr(currentFname, 1, nchar(currentFname)-4)
newFname <- paste("remapped/m", baseFname, "csv", sep=".")
write.csv(current.obj, file=newFname)
}
#########################################################
# Example usage: 
# Rscript GLSeq.remap.R  /mnt/bigdata/processed_data/a3_rna-seq/Tsato/GLBRCY73.EB4.01/GLBRCY73.EB4.01.counts/  /mnt/bigdata/processed_data/a3_rna-seq/Tsato/yeast_remap.csv 1



