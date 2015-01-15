#################################
# Oleg Moskvin; info@scienceforever.com
# Post-run summarizing of the results / file locations 
# for upload to a database 
# Dec. 04, 2013
# Run: Rscript GLSeq.runSummary.R dest.dir.base text.add 
# dest.dir.base may be diferent from the base destination directory at the runtime! 
# this makes the script useful for "wayback" summaries, after moving the results folder to a permanent location
#
# running: Rscript GLSeq.runSummary.R text.add newDestDir outName workflow_template_id experiment_id goGLOW new_run_dir raw_files_table workflow_id refGenID
# i.e. Rscript GLSeq.runSummary.R GLBRCY128.EF4.09 /mnt/bigdata/processed_data/a3_rna-seq/Tsato/GLBRCY128.EF4.09 GLBRC128.ACSH.QC.RSEM (tempID) (experID) 0 /home/GLBRCORG/omoskvin/run/ GLBRCY128.raw.txt 3 2
# Rscript GLSeq.runSummary.R GLBRCY128.EF4.09 /mnt/bigdata/processed_data/a3_rna-seq/Tsato/GLBRCY128.EF4.09 Y128.ACSH.QC.RSEM 1 4 0 /home/GLBRCORG/omoskvin/run/ Y128.cellulose.raw.txt 3 2
#################################
################## 'file_sample_table.txt' is expected to reside in the working directory
#################################
args <- commandArgs(trailingOnly = TRUE)
text.add <- as.character(args[1])
#
# the final destinaton directory (for long-term storage) where the results should already reside at the moment of the script invocation:
newDestDir <- as.character(args[2])
#
# name of the summary table (may be chosen to be more meaningful than the text.add)
reportName <- as.character(args[3])
#
# workflow template ID:
workflowTemplateID <- as.character(args[4])
#
# experiment ID:
experimentID <- as.numeric(args[5])
# (if experimentID is set to 0, then it should be extracted from either local ID table ("expSampleIDs.txt") or taken from GLOW - see below) 
#
# stop at generating and saving of the report table (0) or go ahead to generate XML file and send an update request to GLOW (1)? 
# (will be converted to boolean):
goGLOW <- as.logical(as.numeric((args[6])))
#
# setting the new run directory for cases 
# when the rda file with run parameters contains no-longer-existent temporary directory for the run assigned to 'base.dir'
newRunDir <- as.character(args[7])
#
# file with raw file names actually used for the run (i.e. concateneated with the library names,
# according to former JGI standards):
rawFnames <- as.character(args[8])
#
# workflow ID (0 if a new workflow):
wfID <- as.character(args[9])
#
# reference genome ID:
refGenID <- as.character(args[10])
#
# storing the raw file names' logfiles under newRunDir/rawfnames/
rawFnames <- paste("rawfnames/", rawFnames, sep="")
#
source("GLSeq.Util.R")
currentRun.dataFile <- paste("GLSeq.vars.", text.add, ".rda", sep="")
load(currentRun.dataFile) 
require(XML)
#
currentRun.dataFile.updated <- paste("GLSeq.vars.", text.add, ".updated.rda", sep="")
try(load(currentRun.dataFile.updated), silent=TRUE)
# updating directories:
dirUpdate(newRunDir, newDestDir)
#
# sample ID retrieval:
fileSampleTable <- read.delim(paste(base.dir,"expSampleIDs.txt", sep=""), as.is=TRUE)
# extracting fastq file names from the full path: 
# fastQ.all <- as.character(sapply(fileSampleTable$file.path, path2file))
#
fastQ.all <- as.character(fileSampleTable$File.Name)
sampleID.key <- as.data.frame(cbind(fastQ.all, Sample.ID=fileSampleTable$Sample.ID, Experiment.Set.ID=fileSampleTable$Experiment.Set.ID))
sampleID.key[,1] <- as.character(sampleID.key[,1])
sampleID.key[,2] <- as.numeric(as.character(sampleID.key[,2]))
sampleID.key[,3] <- as.numeric(as.character(sampleID.key[,3]))
# 
# overwriting the relevant range of the raw file within sampleID.key
# with older (actually used and recorded in the rda file of the run) raw file names: 
# (if this procedure is not applicable, the 8-th argument of the script should be set to NULL)
if (!(is.null(rawFnames))) {
rawnames <- as.character(read.table(rawFnames)$V9)
for (i in 1:length(rawnames)) {
sampleID.key[grep(substr(rawnames[i], 6, nchar(rawnames[i])), sampleID.key[,1]),1] <- rawnames[i]
}
} # if rawFnames
#
completionStatus <- "Complete"
#
#####################
#taking exp/samp IDs from GLOW:
#####################
if (goGLOW) {
# dummy assignment (for now):
	a <- 2
}
#
#####################
#generating workflow parameter file:
#####################
attrlist <- c("space1", "space2", "space1","updateFromDb", "dataPrepare", "exprRun", "resCollect", "expID", "protID", "space1", "space3", "space1", "raw.dir", "readyData.dir", "base.dir", "raw.fNames", "strain", "rGenome", "refFASTAname", "refGFFname", "gtfFeatureColumn", "idAttr", "paired.end", "qScores", "nCores", "nStreams", "nStreamsDataPrep", "compConf", "rawGen", "qAlgor", "genobam", "libstrand", "strandBase", "strandExtract", "dest.dir.base", "runAttempt", "libNchar", "libList", "runDate", "trimPath", "picardToolsPath", "fastqcPath", "bwaPath", "readTrim", "trimHead", "artificial.fq", "revisionNumber")
#
space1 <- "#################################################################################"
space2 <- "Attributes used as command line arguments to run GLSeq.top.R script:"
space3 <- "Attributes recorded in GLSeq.attr.R file or in the database (see updateFromDb value to distinguish between the two) at the moment of GLSeq.top.R invocation"
#
# restricting to the objects available for that run because of continuous evolution of the attribute list:
attrlist <- attrlist[attrlist %in% ls()]
#
# output table of the attributes:
attrTable <- matrix(NA, length(attrlist), 2)
rownames(attrTable) <- attrlist
attrTable <- attrTable[,-2,drop=FALSE]
for (ii in 1:length(attrlist)) { 
	if(is.null(get(attrlist[ii]))) attrTable[ii,1] <- "NULL" # straightforward assignment of the get() output doesn't work in this case
#	if(!(is.null(get(attrlist[ii]))) & class(get(attrlist[ii])) != "try-error") attrTable[ii,1] <- get(attrlist[ii]) 
	if(!(is.null(get(attrlist[ii])))) attrTable[ii,1] <- get(attrlist[ii]) }
#
attrTable.fName <- paste(destDirLog, text.add, ".runParam.txt", sep="")
write.table(attrTable, file=attrTable.fName, row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")
#
# final housekeeping (those file may or may not be left there at this moment): 
setwd(newDestDir)
try(system("rm *.ready.txt"), silent = TRUE)
try(system("rm *.fastq.ready"), silent = TRUE)
try(system("rm runLogFile1"), silent = TRUE)
try(system("rm *.S1.bam"), silent = TRUE)
try(system("rm *.S2.bam"), silent = TRUE)
try(system("rm *.A1.bam"), silent = TRUE)
try(system("rm *.A2.bam"), silent = TRUE)
try(system("rm *.transcript.bam"), silent = TRUE)
try(system("rm *.genome.bam"), silent = TRUE)
try(system("rm *.DataReady"), silent = TRUE)
try(system("rm *completeExpression*"), silent = TRUE)
#
# protocolReport <- matrix(NA, 1, 100) 
#
####################
# Collecting library-centric files:
####################
# 
# charchter strings to match library-related files:
libMatch <- substr(fqfiles.table[,1],1,libNchar)
# adding wildcard for the rest of the file name:
libMatch.wild <- paste(libMatch, "*", sep="")
fCollect <- list()
#
##############
# Extracting experiment ID from correspondence table
##############
if (experimentID == 0) {
expIDs <- sampleID.key[substr(sampleID.key[,1],1,libNchar) %in% libMatch,3]
if (length(unique(expIDs)) > 1) stop('simultaneous reporting of samples that belong to multiple experiment sets are not supported at this time \n')
experimentID <- expIDs[1]
}
# 
####################
# Defining output XML files:
####################
#
experimentSamplesXML <-  xmlNode("experiment", attrs=c(experiment_id=experimentID))
workflow.name <- paste(reportName, text.add, sep=".")
experimentWorkflowXML <- xmlNode("experiment", attrs=c(experiment_id=experimentID))
# workflowXML <- xmlNode("workflow", attrs=c(workflow_id=workflow.name), xmlNode("name", workflow.name), xmlNode("workflow_template_id", workflowTemplateID), xmlNode("completion_status", completionStatus))
if (wfID ==0) workflowXML <- xmlNode("workflow", xmlNode("name", workflow.name), xmlNode("workflow_template_id", workflowTemplateID), xmlNode("completion_status", completionStatus), xmlNode("parameter_list", attrTable.fName), xmlNode("experiment_id", experimentID), xmlNode("reference_genome_id", refGenID))
if (wfID !=0) workflowXML <- xmlNode("workflow", attrs=c(workflow_id=wfID), xmlNode("name", workflow.name), xmlNode("workflow_template_id", workflowTemplateID), xmlNode("completion_status", completionStatus), xmlNode("parameter_list", attrTable.fName), xmlNode("experiment_id", experimentID), xmlNode("reference_genome_id", refGenID))
#
for (libs in 1:nrow(fqfiles.table)) {
lib.fNames <- fileNamesCollect(libMatch.wild[libs], newDestDir)
fCollect[[libs]] <- lib.fNames
# XML generation:
	fastq.key.IND <- grep(fqfiles.base[libs], sampleID.key[,1])
	samp.id <- sampleID.key[fastq.key.IND,2][1] # the latter is added for the time being
	sampleXML <- xmlNode("sample", attrs=c(sample_id=samp.id))
	for (thisfile in 1:length(lib.fNames)) {
		fType <- NULL
		fDescr <- "no description provided"
		fName <- path2file(lib.fNames[thisfile])
		if (sum(grep(".bam$",lib.fNames[thisfile], perl=TRUE)) > 0) fType <- "BAM" 
		if (sum(grep(".bam.bai",lib.fNames[thisfile])) > 0) fType <- "BAMI" 
		if (sum(grep(".genes.results",lib.fNames[thisfile])) > 0) fType <- "Counts"
		if (sum(grep(".counts.txt",lib.fNames[thisfile])) > 0) fType <- "Counts"
			if (sum(grep(".counts",lib.fNames[thisfile])) > 0) fType <- "Counts" # BWA-HTSeq version
		if (sum(grep(".wig$",lib.fNames[thisfile], perl=TRUE)) > 0) fType <- "WIG"	
		if (sum(grep(".bw$",lib.fNames[thisfile], perl=TRUE)) > 0) fType <- "BigWig"
		if (sum(grep(".F.bw$",lib.fNames[thisfile], perl=TRUE)) > 0) fType <- "ForwardBigWig"
		if (sum(grep(".R.bw$",lib.fNames[thisfile], perl=TRUE)) > 0) fType <- "ReverseBigWig"
		# file descriptions:
		if (sum(grep(".sorted.bam$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Binary alignment map sorted by coordinate"
		if (sum(grep(".sorted.bam.bai$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Index of the binary alignment map that was sorted by coordinate"
		if (sum(grep(".sorted.F.bam$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Binary alignment map restricted to the Forward strand of the genome, sorted by coordinate"
		if (sum(grep(".sorted.F.bam.bai$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Index of the binary alignment map restricted to the Forward strand of the genome, sorted by coordinate"
		if (sum(grep(".sorted.R.bam$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Binary alignment map restricted to the Reverse strand of the genome, sorted by coordinate"
		if (sum(grep(".sorted.R.bam.bai$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Index of the binary alignment map restricted to the Reverse strand of the genome, sorted by coordinate"
		if (sum(grep(".sorted.wig$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "wiggle file of genome coverage on both strands"	
		if (sum(grep(".sorted.F.wig$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "wiggle file of genome coverage restricted to the Forward strand of the genome"
		if (sum(grep(".sorted.R.wig$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "wiggle file of genome coverage restricted to the Reverse strand of the genome"	
		if (sum(grep(".sorted.bw$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "BigWig file of genome coverage on both strands"
		if (sum(grep(".sorted.F.bw$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "BigWig file of genome coverage restricted to the Forward strand of the genome"
		if (sum(grep(".sorted.R.bw$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "BigWig file of genome coverage restricted to the Reverse strand of the genome"
		if (sum(grep(".genome.sorted.F.bam$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Annotated transcriptome based binary alignment map sorted by coordinate and restricted to the forward strand of the genome"
		if (sum(grep(".genome.sorted.F.bam.bai$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Index of the annotated transcriptome based binary alignment map sorted by coordinate and restricted to the forward strand of the genome"
		if (sum(grep(".genome.sorted.R.bam$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Annotated transcriptome based binary alignment map sorted by coordinate and restricted to the reverse strand of the genome"	
		if (sum(grep(".genome.sorted.R.bam.bai$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Index of the annotated transcriptome based binary alignment map sorted by coordinate and restricted to the reverse strand of the genome"	
		if (sum(grep(".genome.sorted.bam",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Annotated transcriptome based binary alignment map sorted by coordinate"
		if (sum(grep(".genome.sorted.bam.bai$",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Index of the annotated transcriptome based binary alignment map sorted by coordinate"
		if (sum(grep(".transcript.sorted.bam",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Annotated transcriptome based binary alignment map based on the transcripts inferred by RSEM and sorted by coordinate"
		if (sum(grep(".transcript.sorted.bam.bai",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Index of the annotated transcriptome based binary alignment map based on the transcripts inferred by RSEM and sorted by coordinate"
		if (sum(grep(".genome.sorted.F.wig",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Wiggle file of the annotated transcriptome coverage restricted to the Forward strand of the genome"
		if (sum(grep(".genome.sorted.R.wig",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Wiggle file of the annotated transcriptome coverage restricted to the Reverse strand of the genome"
		if (sum(grep(".genes.results",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "Library centric expression matrix with different measures in separate columns"
		if (sum(grep(".genome.sorted.R.bw",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "BigWig file of the annotated transcriptome coverage restricted to the Reverse strand of the genome"
		if (sum(grep(".genome.sorted.F.bw",lib.fNames[thisfile], perl=TRUE)) > 0) fDescr <- "BigWig file of the annotated transcriptome coverage restricted to the Forward strand of the genome"
		# assembling XML
		if (!(is.null(fType))) {
		outfileXML <- xmlNode("datafile", xmlNode("name", fName),  xmlNode("description", fDescr), xmlNode("file_path", lib.fNames[thisfile]), xmlNode("file_type", fType), xmlNode("owner", "omoskvin"), xmlNode("input_output", "output"))
		# adding sample >> files tree to XML for all the samples (sampleXML):
		sampleXML <- addChildren(sampleXML, outfileXML)
		# adding datafiles part to the workflow XML: 
		workflowXML <- addChildren(workflowXML, outfileXML)	
		} # avoiding records for the file types outside of the controlled vocabulary
	} # for thisfile    
experimentSamplesXML <- addChildren(experimentSamplesXML, sampleXML)	
} # for libs
#
#####################
# adding workflow-specific files 
# to workflow.XML:
#####################
wfMatch.wild <- paste(text.add, "*", sep="")
wf.fNames <- fileNamesCollect(wfMatch.wild, newDestDir)
all.fNames <-  fileNamesCollect("*", newDestDir) # revison 78
wf.fNames <- all.fNames # just like that, for the time being
for (wfFile in 1:length(wf.fNames)) {
fType <- NULL
fDescr <- "no description provided"
fName <- path2file(wf.fNames[wfFile])
if (sum(grep(".csv$",wf.fNames[wfFile], perl=TRUE)) > 0 & qAlgor == "RSEM" ) fType <- "Merged Counts" # catch-all
	# fro RSEM:
if (sum(grep(".FPKM.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts" 
	fDescr <- "Normalized counts FPKM" }
if (sum(grep(".FPKM_lower.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts"
	fDescr <- "Lower estimate of FPKM at 95 percent confidence level" }
if (sum(grep(".FPKM_upper.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts"
	fDescr <- "Upper estimate of FPKM at 95 percent confidence level" }
if (sum(grep(".FPKM_pme.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts"
	fDescr <- "Posterior Mean Estimate of FPKM" }
if (sum(grep(".counts.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts" 
	fDescr <- "Counts" }
if (sum(grep(".counts_pme.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts" 
	fDescr <- "Posterior Mean Estimate of counts" }
if (sum(grep(".TPM.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts"
	fDescr <- "Transcripts Per Million" }
if (sum(grep(".TPM_pme.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts"
	fDescr <- "Posterior Mean Estimate of Transcripts Per Million" }
if (sum(grep(".TPM_upper.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts" 
	fDescr <- "Upper estimate of Transcripts Per Million at 95 percent confidence level" }
if (sum(grep(".TPM_lower.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts" 
	fDescr <- "Lower estimate of Transcripts Per Million at 95 percent confidence level" }
	# for BWA-HTSeq:
if (sum(grep(".RPKM.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "Merged Counts" 
	fDescr <- "Normalized counts RPKM" }
if (sum(grep(".exceptionReport.csv$",wf.fNames[wfFile], perl=TRUE)) > 0) { fType <- "QC" 
	fDescr <- "Read alignment statistics for every library in the workflow" }
if (sum(grep("^m.",wf.fNames[wfFile], perl=TRUE)) > 0) fDescr <- paste(fDescr, "remapped to meaningful sample identifiers") # doesn't work when called via Rscript!
if (!(is.null(fType))) {
	wfFileXML <- xmlNode("datafile", xmlNode("name", fName), xmlNode("description", fDescr), xmlNode("file_path", wf.fNames[wfFile]), xmlNode("file_type", fType), xmlNode("owner", "omoskvin"), xmlNode("input_output", "output"))
	workflowXML <- addChildren(workflowXML, wfFileXML)	
	}
}
#
#####################
# adding to final XML files: 
#####################
experimentWorkflowXML <- addChildren(experimentWorkflowXML, workflowXML)
#
####################
# exporting XML's to files:
####################
setwd(base.dir)
expSamples.xmlName <- paste("xml/", workflow.name, ".expSamples.xml", sep="")
expWorkflow.xmlName <- paste("xml/", workflow.name, ".expWorkflow.xml", sep="")
saveXML(experimentSamplesXML, file=expSamples.xmlName, prefix=NULL)
saveXML(experimentWorkflowXML, file=expWorkflow.xmlName, prefix=NULL)
#
#####################
# updating the GLOW:
#####################
expSamples.upsertName <- paste("xml/", workflow.name, ".expSamples.upsert", sep="")
expWorkflow.upsertName <- paste("xml/", workflow.name, ".expWorkflow.upsert", sep="")
#
# generating upsert command files:
system(paste("cat xml/expxml.start", expSamples.xmlName, " >  xml/expSamples.upsertTMP2"))	
system(paste("perl -pe 's/>\n/>/g'  xml/expSamples.upsertTMP2 >  xml/expSamples.upsertTMP"))
system(paste("perl -pe 's/>\\s+/>/g' xml/expSamples.upsertTMP >", expSamples.upsertName)) 
system(paste("rm  xml/expSamples.upsertTMP && mv", expSamples.upsertName, "xml/expSamples.upsertTMP"))
system(paste("perl -pe 's/\n//g' xml/expSamples.upsertTMP >", expSamples.upsertName))
system("rm xml/expSamples.upsertTMP")
system("rm xml/expSamples.upsertTMP2")
#
system(paste("cat xml/expxml.start", expWorkflow.xmlName, " >  xml/expWorkflow.upsertTMP2"))	
system(paste("perl -pe 's/>\n/>/g'  xml/expWorkflow.upsertTMP2 >  xml/expWorkflow.upsertTMP"))
system(paste("perl -pe 's/>\\s+/>/g' xml/expWorkflow.upsertTMP >", expWorkflow.upsertName)) 
system(paste("rm  xml/expWorkflow.upsertTMP && mv", expWorkflow.upsertName, "xml/expWorkflow.upsertTMP"))
system(paste("perl -pe 's/\n//g' xml/expWorkflow.upsertTMP >", expWorkflow.upsertName))
system("rm xml/expWorkflow.upsertTMP")
system("rm xml/expWorkflow.upsertTMP2")
#
if (goGLOW) {
sa.comm <- paste("curl --cookie ~/GLOWtest/cjar --data @", expSamples.upsertName, " https://glow-ng-test.glbrc.org/upsert_experiment", sep="")
wf.comm <- paste("curl --cookie ~/GLOWtest/cjar --data @", expWorkflow.upsertName, " https://glow-ng-test.glbrc.org/upsert_experiment", sep="")
system(paste(sa.comm, "&&", wf.comm)) 
}
#####################################









