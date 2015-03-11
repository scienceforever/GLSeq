#'  A function to rename fastq files to satisfy GLSEQ conventon
#'
#' \code{rename.fastq.files} renames files in a given directory to satisfy GLSEQ naming conventon
#' 
#' @param fastq.dir char path to fastq folder
#' @param left.read.file.extension char original file extention of left read files, to be replaced by .1.fq
#' @param right.read.file.extension char original file extention of right read files, to be replaced by .2.fq
#' @return 1
#' @export
rename.fastq.files = function(fastq.dir,left.read.file.extension,right.read.file.extension) {
  #  remember current directory
  current.dir = getwd()
  #  go to the fastq directory
  setwd(fastq.dir)
  #  get the original file names
  original.file.names = dir()
  #  strip the extensions
  new.file.names = sub(left.read.file.extension,"",original.file.names,fixed=T)
  new.file.names = sub(right.read.file.extension,"",new.file.names,fixed=T)
  #  pad shorter file names to make them the same length as the longest one
  fn.lengths = nchar(new.file.names)
  pad.length = max(fn.lengths) - fn.lengths
  pad.strings = sapply(pad.length,function(x){paste(rep("_",x),collapse="")})
  new.file.names = paste(new.file.names,pad.strings,sep="")
  #  add appropriate extensions
  left.file.ind = grep(left.read.file.extension,original.file.names)
  new.file.names[left.file.ind] = paste(new.file.names[left.file.ind],"1.fq",sep=".")
  right.file.ind = grep(right.read.file.extension,original.file.names)
  new.file.names[right.file.ind] = paste(new.file.names[right.file.ind],"2.fq",sep=".")
  #  rename the files
  names = data.frame(old=original.file.names,new=new.file.names)
  apply(names,1,function(x){file.rename(x["old"],x["new"])})
  #  go back to the original directory
  setwd(current.dir)
  1
}
