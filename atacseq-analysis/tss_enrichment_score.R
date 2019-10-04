# This script calculates the TSS enrichment score from the TSS regions given an ATAC-seq signal.
# It takes three inputs.
# 1st: The path that leads to the input file locations.
# 2nd: A matrix obtained by deeptools computeMatrix, which contains the cumulative TSS signal coverage.
# 3rd: The path that leads to the output file locations.

# Read in computeMatrix output matrix, remove the header section, convert it onto a matrix R object, and change "NA" to zeroes.
args = commandArgs(trailingOnly=TRUE)
inputFile <- paste0(args[1], args[2])
prefix <- gsub(".matrix", "", args[2])
TSSmatrix <- read.delim(inputFile, skip=3, header=F)
TSSmatrix <- as.matrix(TSSmatrix)
TSSmatrix[is.na(TSSmatrix)] <- 0

# Get the mean of the columns of the matrix.
TSSmatrixColMeans <- colMeans(TSSmatrix)

# Get the mean of the first 200 columns means.
backgroundMean <- mean(TSSmatrixColMeans[1:200])

# Normalize the means for the entire 4000 bp by dividing the values by the background mean.
TSSmatrixNorm <- TSSmatrixColMeans/backgroundMean

# Find the highest point in the plot. This will be considered the TSS enrichment score
maxPoint <- max(TSSmatrixNorm)
maxPoint <- format(round(maxPoint, 2), nsmall = 2)
maxPointPos <- which.max(TSSmatrixNorm)

# Save plots onto a pdf file.
pdfName <- paste0(args[3], prefix, ".pdf")
pdf(pdfName, width = 6.66, height = 5, useDingbats = FALSE)

# Plot the normalized means for the aggregate TSS +/- 2000 bp window.
par(mfrow = c(1,2),oma = c(0, 0, 6, 0))
plot(TSSmatrixColMeans, 
     type = "l", 
     xlab = "bp position relative to TSS", 
     ylab = "raw mean coverage", 
     xaxt = "n", 
     lwd = 2)
abline(v = maxPointPos, 
       col = "red", 
       lty = 5, 
       lwd = 3)
axis(1, at = c(0,1000,2000,3000,4000), labels=c(-2000,-1000,"TSS",1000,2000))
plot(TSSmatrixNorm, 
     type = "l", 
     xlab = "bp position relative to TSS", 
     ylab = "coverage normalized relative to background", 
     xaxt = "n", 
     lwd = 2)
abline(v = maxPointPos, 
       col = "red", 
       lty = 5, 
       lwd = 3)
axis(1, at = c(0,1000,2000,3000,4000), labels=c(-2000,-1000,"TSS",1000,2000))
mtext(paste0("meta-TSS coverage plots\n", prefix, "\nTSS enrichment score: ", toString(maxPoint)), 
      outer = TRUE, 
      cex = 2)
dev.off()

# Save the TSS enrichment score as a text file.
outFile <- paste0(args[3], prefix, ".score")
writeLines(paste0(args[2], "\n", toString(maxPoint)), outFile)