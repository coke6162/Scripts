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
bamFileName <- gsub(".bam", "", basename(bamFile))

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

# Adjust the read start sites
tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
which <- as(seqinfo(Btaurus), "GRanges")
gal <- readBamFile(bamFile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
gal1 <- shiftGAlignmentsList(gal)

# Generate promoter/transcript body (PT) score
pdfName_promTran <- paste0(args[2], prefix, "_PT.pdf")
pdf(pdfName_promTran, width = 6.66, height = 5, useDingbats = FALSE)

txs <- transcripts(TxDb.Btaurus.UCSC.bosTau8.refGene)
pt <- PTscore(gal1, txs)
plot(pt$log2meanCoverage, pt$PT_score, 
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")

dev.off()

# Generate nucleosome free regions (NFR) score
pdfName_NFR <- paste0(args[2], prefix, "_NFR.pdf")
pdf(pdfName_NFR, width = 6.66, height = 5, units = in, useDingbats = FALSE)

nfr <- NFRscore(gal1, txs)
plot(nfr$log2meanCoverage, nfr$NFR_score, 
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))

dev.off()

# Generate transcription start site (TSS) enrichment score
tsse <- TSSEscore(gal1, txs)
tsse_summary <- summary(tsse$TSS.enrichment.score)
write.table(tsse, paste0(args[2], prefix, "_tssEnrich.txt"), sep = "\t")
write.table(tsse, paste0(args[2], prefix, "_tssEnrich_summary.txt"), sep = "\t")

sessionInfo()