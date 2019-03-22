# SpTPM.R
#' Scatterplot of TPM of inputs vs. average TPM of sample.
#'
#' \code{SpTPM} Scatterplot of TPM of inputs vs. average TPM of sample. Count reads per trancript and get the TPM; plot the scatterplot for average TPM samples and TPM inputs.
#'
#' @param inputFile The input of alignment files as a list in bam form
#' @param sampleFiles The samples of alignment files as a list in bam form
#' @param version default is "hg19". Could be "mm10", "mm9", "hg38" and "hg19"
#'
#' @examples
#' \dontrun{
#' inputFile <- "../BCBProjectData/SRR4393137.Aligned.sortedByCoord.out.bam"
#' sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam",
#' "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam",
#' "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")
#' SpTPM(inputFile, sampleFiles)
#' }
#' @export
SpTPM <- function(inputFile, sampleFiles, version = "hg19") {
  print("Plot CLIP average TPM")

  # Initial a data frame to store the sequence count and gene length
  TPMDF <- data.frame(Input1 = double(), stringsAsFactors = FALSE)
  bamFiles <- c(inputFile, sampleFiles)
  for (i in seq_along(bamFiles)) {
    # 1. get the counts from Rsubread::featureCounts
    fc <- Rsubread::featureCounts(files=bamFiles[i],annot.inbuilt="hg19")
    x <- edgeR::DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

    # 2. store the data to data frame
    if (i == 1) {
      TPMDF <- data.frame(geneLen = x$genes$Length, stringsAsFactors = FALSE)
    }

    b <- data.frame(input = x$counts,stringsAsFactors = FALSE)
    colnames(b) <- paste0("Input", i)
    TPMDF <- cbind(TPMDF, b)
  }

  # Get the TPM for each gene
  for (i in seq(2,length(bamFiles)+1, 1)) {
    # 1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
    TPMDF[,i] <- TPMDF[,i] / (TPMDF[,1]/1000.0)
    # 2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
    sum <- sum(TPMDF[,i])/1000000.0
    # 3. Divide the RPK values by the “per million” scaling factor. This gives you TPM.
    TPMDF[,i] <- TPMDF[,i]/sum
  }

  # save TPMDF
  # save(TPMDF, file = "./inst/extdata/TPMDF.Rdata")
  # load(system.file("extdata", "TPMDF.Rdata", package = "BCBProject"))

  # Get the average of the samples
  TPMDF$sum <- 0
  for (i in seq(3,length(sampleFiles)+2, 1)) {
    TPMDF$sum <- TPMDF$sum + TPMDF[,i]
  }

  TPMDF$average <- TPMDF$sum/length(sampleFiles)

  # case 1: plot all data, add 1 for log value
  x <- TPMDF[,2] + 1
  y <- TPMDF$average + 1

  plot(x, y,
      xlab="Input average TPM",
      ylab="CLIP average TPM",
      log="xy")

  reg = lm(log(x) ~ log(y))
  abline(reg, col=2)
  legend("topright", bty="n", legend=paste("R =", format(summary(reg)$adj, digits=4)))

  # case2: plot non-zero data
  x <- TPMDF[TPMDF[,2]!=0 & TPMDF$average!=0,2]
  y <- TPMDF$average[TPMDF[,2]!=0 & TPMDF$average!=0]

  plot(x, y,
    xlab="Input non-zero average TPM",
    ylab="CLIP non-zero average TPM",
    log="xy")

  reg = lm(log(x) ~ log(y))
  abline(reg, col=2)
  legend("topright", bty="n", legend=paste("R =", format(summary(reg)$adj, digits=4)))
}

# [END]
