# SpRPKM.R
#' Scatterplot of RPKM of inputs vs. average RPKM of sample.
#'
#' \code{SpRPKM} Scatterplot of RPKM of inputs vs. average RPKM of sample. Count reads per trancript and get the RPKM; plot the scatterplot for average RPKM samples and RPKM inputs.
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
#' SpRPKM(inputFile, sampleFiles)
#' }
#' @export
SpRPKM <- function(inputFile, sampleFiles, version = "hg19") {
  print("Plot CLIP average RPKM")

  # Initial a data frame to store the sequence RPKM
  RPKMDF <- data.frame(Input1 = double(), stringsAsFactors = FALSE)
  bamFiles <- c(inputFile, sampleFiles)
  for (i in seq_along(bamFiles)) {
    # 1. get the counts from Rsubread::featureCounts
    # invisible(utils::capture.output(Rsubread::funciton()))
    if (version == "hg19") {
      fc <- Rsubread::featureCounts(files=bamFiles[i],annot.inbuilt="hg19")
    } else {
      fc <- Rsubread::featureCounts(files=bamFiles[i],annot.inbuilt=version)
    }
    x <- edgeR::DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

    # 2. get rpkm for each bam file
    x_rpkm <- edgeR::rpkm(x,x$genes$Length)

    # 3. store the data to data frame
    if (i == 1) {
      colnames(x_rpkm) <- paste0("Input",i)
      RPKMDF <- rbind(RPKMDF, x_rpkm)
    } else {
      RPKMDF[,paste0("Input",i)] <- x_rpkm
    }
  }

  # save RPKMDF
  # save(RPKMDF, file = "./inst/extdata/RPKMDF.Rdata")
  # load(system.file("extdata", "RPKMDF.Rdata", package = "BCBProject"))

  # Get the average of the samples
  RPKMDF$sum <- 0
  for (i in seq(2,length(sampleFiles)+1, 1)) {
    RPKMDF$sum <- RPKMDF$sum + RPKMDF[,i]
  }

  RPKMDF$average <- RPKMDF$sum/length(sampleFiles)

  # case 1: plot all data, add 1 for log value
  x <- RPKMDF[,1] + 1
  y <- RPKMDF$average + 1

  plot(x, y,
    xlab="Input average RPKM",
    ylab="CLIP average RPKM",
    log="xy")
  reg = lm(log(x) ~ log(y))
  abline(reg, col=2)
  legend("topright", bty="n", legend=paste("R =", format(summary(reg)$adj, digits=4)))

  # case 2: plot data with non-zero value in input and samples
  x <- RPKMDF[RPKMDF[,1]!=0 & RPKMDF$average!=0,1]
  y <- RPKMDF$average[RPKMDF[,1]!=0 & RPKMDF$average!=0]

  plot(x, y,
    xlab="Input non-zero average RPKM",
    ylab="CLIP non-zero average RPKM",
    log="xy")
  reg = lm(log(x) ~ log(y))
  abline(reg, col=2)
  legend("topright", bty="n", legend=paste("R =", format(summary(reg)$adj, digits=4)))
}

# [END]
