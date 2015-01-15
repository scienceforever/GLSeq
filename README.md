# GLSeq: Pipeline for low-level RNA-Seq data processing 

![GLSeq scheme](GLSeq.gif)

**The overview of GLSeq**. _Supplementary figure, Moskvin et. al. (2014) Preprint doi: [10.1101/010488](http://dx.doi.org/10.1101/010488) 



----
## Introduction

**GLSeq** ("GL" stands for "Great Lakes") is a pipeline facilitating low-level processing of the RNA-Seq data. This includes (optional) trimming of the sequencing reads, genome (or transcriptome) alignment and feature coverage counting, using either hard-threshold or probabilistic approaches.  

The pipeline was used to compare low-level data processing options in our recent study (see the reference below). At the moment, it is focused on providing gene-level information (except for annotation-agnostic read coverage visualizations when genome alignemnt option is chodsen). Soon, it will be extended to generate experiment-wide summaries of transcript-level results useful for alternative splicing studies. 

The pipeline is written by biologist for biologists, using R language. Unlike usual R scripts, the GLSeq scripts are not intended to be loaded into an R session; they must be used as command-line scripts with options, on top of "Rscript" command (see operation manual). 

## Features
- support for single- or paired-end libraries
- support for strand-specific libraries 
- easy setup of multithreading by the end-user via providing the number of computation streams and CPUs per stream in the run-specific attributes file
- unattended processing of large sets of libraires regardless of the available computational resources
- option to run splitting by strand coverage when using concatenated FASTQ files for paired-end libraires (at the "data preparation" step)
- generating both library-specific count files and experiment-centric tables of counts, FPKM etc. 
- providing 3 types of visualization files (bam, WIG, BigWIG), with strand-specific versions, where applicable

## Prerequisites

- Linux workstation 
- R [R-Project](http://r-project.org)
- Ruby [Ruby Language](https://www.ruby-lang.org/en/)
- Perl [Perl Language](https://www.perl.org/)
- Java [Oracle Java](https://www.oracle.com/java/index.html)
- Python [Python Language](https://www.python.org/)

## Citing GLSeq
If you use the GLSeq package, please cite:

- Moskvin O.V., McIlwain S., Ong I.M. Making sense of RNA-Seq data: from low-level processing to functional analysis. doi: [10.1101/010488](http://dx.doi.org/10.1101/010488)

If you use Bowtie-RSEM as part of your customized pipeline, please also cite:

- Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome biology 2009; 10:R25
- Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC bioinformatics 2011; 12:323

If you use BWA-HTSeq as part of the pipeline, please also cite:

- Li H, Durbin R. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics 2009; 25:1754-1760
- Anders SP, P.T.; Huber, W. HTSeq - A Python framework to work with high-throughput sequencing data. doi: [10.1101/002824](http://dx.doi.org/10.1101/002824)
- Picard project [Picard Tools](http://broadinstitute.github.io/picard/)

If you use BWA-FeatureCounts option (coming soon) please also cite:

- Liao Y., Smyth G.K., Shi W. (2014) featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30(7): 923-930

If you use Trimmomatic, please cite: 

- Bolger AM, Lohse M, Usadel B. (2014) Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30:2114-2120



## Usage 

## Under the hood


