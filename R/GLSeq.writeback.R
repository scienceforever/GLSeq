#########################################################
# Great Lakes Seq package for low-level processing of RNA-Seq data
# After-run updating of the datavbase
# this script updates the rda file with the timestamp of run completeness, run status
# and adds all the missing objects introduced later (for earlier runs)
#########################################################
#
args <- commandArgs(trailingOnly = TRUE)
text.add <- as.character(args[1])
source("GLSeq.Util.R")
#
# later-introduced objects that may not be present in every 
# variable-containing file of the earlier runs:
libList <- NULL
runDate <-  NULL
dest.dir.base<-  NULL
expID <-  NULL
medium <- "someMedium"
timepoint <- "someTimepoint"
replicate <- "someReplicate" 
strain <- "MT203pPET"
# (those will be overwritten with the actual values below, if available): 
# loading variables for the current run (from the current directory) 
currentRun.dataFile <- paste("GLSeq.vars.", text.add, ".rda", sep="")
load(currentRun.dataFile) # 
#
# preparing the objects that may be missing in the 'currentRun.dataFile'
# but derived form other objects already present there: 
# in particular, "libfile" may be derived from "fqfiles" by removing the trailing".1.fq" or ".2.fq"
# since every library is mentioned twice in the split file list, only unique records should be kept:
if(is.null(libList)) {
fqfiles.base <- unique(substr(fqfiles, 1,nchar(fqfiles) - 5))
# inferring file names in the raw directory in the case of pre-processed (uncompressed and split) data:
libList <- paste(fqfiles.base, ".gz", sep="") }
#
# capturing the date/time of the processing completeness (i.e. calling of this script): 
runDateEnd <- paste(strsplit(date()," ")[[1]], collapse="_")
runDateEnd <- gsub(":", "_", runDateEnd) 
#
# option for the MySQL update:
if (updateFromDb == "update") {
library(RMySQL) 
dbcon <- dbConnect(MySQL(), group="glseq")
statusUpdate <- paste("update glseqresults set status='Completed', runend=", "'", runDateEnd, "'", " where textadd=", "'", text.add, "'", " and  status= 'Submitted'", sep="")
dbSendQuery(dbcon, statusUpdate)
Sys.sleep(2)
dbDisconnect(dbcon)
} # if update
#
#################
# Option for GLOW update
#################
# 
if (updateFromDb == "GLOW") {
a <- 1 # other development rout was taken, for the time being
} # if GLOW
# saving the updated attributes to the new run data file
# (including attributes introduced after the de facto data processing - will be deprecated later on)
currentRun.dataFile.updated <- paste("GLSeq.vars.", text.add, ".updated.rda", sep="")
save.image(file=currentRun.dataFile.updated)
# 



