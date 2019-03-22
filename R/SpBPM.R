# SpBPM.R
#' Scatterplot of BPM of inputs vs. average BPM of sample.
#'
#' \code{SpBPM} Scatterplot of BPM of inputs vs. average BPM of sample. Get the BPM from deeptools; plot the scatterplot for average BPM samples and BPM inputs.
#'
#' @param inputFile The input of alignment files as a list in bam form
#' @param sampleFiles The samples of alignment files as a list in bam form
#' @param version default is "hg19". Could be "mm10", "mm9", "hg38" and "hg19"
#' @param binSize default is 30. The bin size for each bin
#' @param mergeFile default "merged.bam". Merged file of all samples to do bamCompare function
#'
#' @examples
#' \dontrun{
#' inputFile <- "../BCBProjectData/SRR4393137.Aligned.sortedByCoord.out.bam"
#' sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam",
#' "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam",
#' "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")
#' SpBPM(inputFile, sampleFiles)
#' }
#' @export
SpBPM <- function(inputFile, sampleFiles, version = "hg19", binSize = 30, mergeFile = "merged.bam") {
  print("Plot CLIP average BPM")
  # Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/Library/Frameworks/Python.framework/Versions/3.6/bin", sep=":"))

  # Get the BPM for each bin from multiBamSummary
  print("multiBamSummary - Note: took very long time and will output a file results.bed that store the counts for each bin")
  bamFiles <- c(inputFile, sampleFiles)
  command <- "multiBamSummary bins -b "
  for (file in bamFiles) {
    command <- paste0(command, file, " ")
  }
  command <- paste0(command, "-o results.npz -bs ", binSize, " --outRawCounts results.bed")
  system(command)

  # Read the multiRead file from multiBamSummary
  multiReads <- utils::read.table("results.bed",
                            header = FALSE,
                            sep="\t",
                            stringsAsFactors=FALSE)

  # save the data and load if necessary
  # save(multiReads, file = "./inst/extdata/multiBamSummary.Rdata")
  # load(system.file("extdata", "multiBamSummary.Rdata", package = "BCBProject"))

  # Get the BPM for all inputs
  # BPM (per bin) = number of reads per bin / sum of all reads per bin (in millions)
  for (i in seq(4,length(bamFiles)+3, 1)) {
    sum <- sum(multiReads[,i])
    multiReads[,i] <- multiReads[,i] / (sum/1000000.0)
  }

  # Get the average
  multiReads$sum <- 0
  for (i in seq(5,length(sampleFiles)+4, 1)) {
    multiReads$sum <- multiReads$sum + multiReads[,i]
  }

  multiReads$average <- multiReads$sum/length(sampleFiles)

  # Plot data with non-zero value in input and samples
  plotData <- multiReads[(multiReads$average != 0) & (multiReads$V4 != 0),]

  x <- plotData$V4
  y <- plotData$average

  plot(x, y,
    xlab="Input non-zero average BPM",
    ylab="CLIP non-zero average BPM",
    log="xy")

  reg = lm(log(x) ~ log(y))
  abline(reg, col=2)
  legend("topright", bty="n", legend=paste("R =", format(summary(reg)$adj, digits=4)))

  # Get log2ratio of BPM from bamcompare
  if(!file.exists(mergeFile)){
    Rsamtools::mergeBam(sampleFiles, mergeFile,overwrite=TRUE)
    Rsamtools::indexBam(mergeFile)
  }

  print("bamCompare - Note: took very long time and will output a file log2ratio.bw that store the score of log of ratio of two files")
  command <- paste0("bamCompare -b1 ", mergeFile, " -b2 ", inputFile, " -o log2ratio.bw -bs ", binSize, " --scaleFactorsMethod None --normalizeUsing BPM")
  system(command)
  t <- rtracklayer::import("log2ratio.bw", format="BigWig")
  plot(sort(t$score), xlab = "index", ylab = "CLIP log2ratio BPM")
}

# [END]
