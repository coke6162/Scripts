# This script calls the ATACseqQC fragSizeDist function 
# It takes two inputs:
# 1st - The name of a bam file.
# 2nd - The output location where to store the pdf with the plot.

# Load required packages
library(ATACseqQC)
library(BSgenome.Btaurus.UCSC.bosTau8)
library(TxDb.Btaurus.UCSC.bosTau8.refGene)
                       
# Read name of input bam and remove last suffix.
args = commandArgs(trailingOnly=TRUE)
prefix <- gsub(".sorted.bam", "", args[1])

# Read in bam file living in the ATACseqQC directory
bamFile <- system.file("extdata", args[1], package="ATACseqQC", mustWork=TRUE)
bamFileName <- gsub(".sorted.bam", "", basename(bamFile))

# Plot estimated library complexity
pdfName_libComplex <- paste0(args[2], prefix, "_libComplex.pdf")
pdf(pdfName_libComplex, width = 6.66, height = 5, useDingbats = FALSE)

libComplex <- estimateLibComplexity(readsDupFreq(bamFile))

dev.off()

# Plot the size distribution from the loaded bam file.
pdfName_fragSizeDist <- paste0(args[2], prefix, "_fragSizeDist.pdf")
pdf(pdfName_fragSizeDist, width = 6.66, height = 5, useDingbats = FALSE)

fragSize <- fragSizeDist(bamFile, prefix)

dev.off()

sessionInfo()