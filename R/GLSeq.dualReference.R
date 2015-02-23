######################################
# preparing "dual reference" gtf file
# Oleg Moskvin - 2013
######################################
#
args <- commandArgs(trailingOnly = TRUE)
single.gtf.name <- as.character(args[1])
double.gtf.name <- as.character(args[2])
single.gtf = read.table(single.gtf.name, sep="\t", as.is=TRUE)
#
exon.boo <- single.gtf$V3 == "exon" 
exon.ind <- (1:length(exon.boo))[exon.boo]
double.gtf <- single.gtf
# adding reverse strands as transcripts: 
for(i in 1:length(exon.ind)) {
current.lineNr <- exon.ind[i]
myLine <- as.character(double.gtf[current.lineNr,])
myLineNew <- myLine
if(is.na(myLine[7])) i
if(myLine[7] == "-") myLineNew[7] <- "+"
if(myLine[7] == "+") myLineNew[7] <- "-"
new.attr <- as.character(myLineNew[9])
attr.pairs <- strsplit(new.attr, split=";")[[1]]
attr.pairs.new <- paste(attr.pairs,"ANTISENSE", sep=".")
attr.line.new <- paste(attr.pairs.new, collapse=";")
myLineNew[9] <- attr.line.new
# adding the second line for annotation purpose: 
myLineNew.second <- myLineNew
# assigning "ncRNA" attribute to the antisense transcript:
myLineNew.second[3] <- "ncRNA"
# handlig the cases when the initial record has or has not its duplicate annotation record
# (we may have slightly different starting positions for the 2 records!!!
# so, easy matching by position won't work here
current.geneName <- strsplit(attr.pairs[1], split=" ")[[1]][2]
nextLine <- as.character(double.gtf[current.lineNr+1,])
# doing the attribute split with the next line: 
next.attr <- as.character(nextLine[9])
next.pairs <- strsplit(next.attr, split=";")[[1]]
next.geneName <- strsplit(next.pairs[1], split=" ")[[1]][2]
# if there is no duplicate record (exception rather than rule): 
if (current.geneName != next.geneName) {
topBrick <- head(double.gtf, current.lineNr)
bottomBrick <- tail(double.gtf, -current.lineNr)
double.gtf <- rbind(topBrick, myLineNew, myLineNew.second, bottomBrick) }
# if there is a duplicate record (typical case): 
if (current.geneName == next.geneName) {
topBrick <- head(double.gtf, current.lineNr+1)
bottomBrick <- tail(double.gtf, -(current.lineNr+1))
double.gtf <- rbind(topBrick, myLineNew, myLineNew.second, bottomBrick) }
# we have just inserted 2 rows; updating indexes for the next cycle:
exon.ind <- exon.ind+2 
} 
#
######
# putting single quotes around textual attributes:
attr.quote.add <- function(attribute.vec) {
for (i in 1: length(attribute.vec)) {
curr.attributes <- attribute.vec[i]
attr.pairs <- strsplit(curr.attributes, split=";")[[1]]
trim.leading <- function (x)  sub("^\\s+", "", x)
attr.pairs <- trim.leading(attr.pairs)
for (j in 1:length(attr.pairs)) {
# current pair of ID and textual attribute:
this <- attr.pairs[j]
# handling IDs and their textual values separately:
this.ID <- strsplit(this, split=" ")[[1]][1]
textual <- strsplit(this, split=" ")[[1]][2]
textual <- paste("\'", textual, "\'", sep="")
this.updated <- paste(this.ID, textual, collapse=" ")
attr.pairs[j] <- this.updated
}
attribute.vec[i] <- paste(attr.pairs, collapse=";")
}
attribute.vec
}
double.gtf[,9] <- attr.quote.add(double.gtf[,9])
write.table(double.gtf, file=double.gtf.name, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")




