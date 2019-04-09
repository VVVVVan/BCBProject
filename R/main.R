# main.R
#' Generate plots to analysis the bam files.
#'
#' \code{main} Generate several plots to analysis the bam files. 1 a. Scatterplot of RPKM of inputs vs. average RPKM of samples; b. Scatterplot of TPM of inputs vs. average TPM of samples; c. Scatterplot of BPM of inputs vs. average BPM of samples with binsize 30. 2 a. gene coverage of samples over annotated regions of genome (default is hg19); b. gene coverage calculated per kilobase of annotated features; c. gene coverage of different length of 5' UTR and 3' UTR. 3 a. Average profile from samples over enriched genes and inenriched genes body by ngs.plot. b. Average profile from samples over 5' UTR; c. Average profile from samples over 3' UTR
#'
#' @param inputFile The input of alignment files as a list in bam form
#' @param sampleFiles The samples of alignment files as a list in bam form
#' @param annotaedGenomes The list of annotated genomes for "5' UTR", "3' UTR", "exons", "introns", "promoters", "lncRNA", "tRNA", "snc / miRNA". The order of the annotated genome is important.
#' @param UTRDiffL UTR feature file with different length in bed format, default if NULL
#' @param version default is "hg19". Could be "mm10", "mm9", "hg38" and "hg19"
#' @param binSize default is 30. The bin size for each bin
#' @param mergeFile default "BPMmerged.bam". Merged file of all samples to do bamCompare function
#'
#' @examples
#' \dontrun{
#' inputFile <- "../BCBProjectData/SRR4393137.Aligned.sortedByCoord.out.bam"
#' sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")
#'
#' annotaedGenomes <- c("../BCBProjectData/5UTR2000.bed", "../BCBProjectData/3UTR2000.bed", "../BCBProjectData/GENCODE.Exon.bed", "../BCBProjectData/GENCODE.Intron.bed", "../BCBProjectData/GENCODE.Promoter.bed", "../BCBProjectData/lincRNA.bed", "../BCBProjectData/tRNA.bed", "../BCBProjectData/miRNA.bed")
#' UTRDiffL <- c("../BCBProjectData/3UTR2000.bed","../BCBProjectData/3UTR1000.bed","../BCBProjectData/3UTR500.bed","../BCBProjectData/3UTR200.bed","../BCBProjectData/5UTR2000.bed","../BCBProjectData/5UTR1500.bed","../BCBProjectData/5UTR1000.bed","../BCBProjectData/5UTR500.bed")
#' main(inputFile, sampleFiles, annotaedGenomes, UTRDiffL)
#' }
#' @export
#'
main <- function(inputFile, sampleFiles, annotaedGenomes, UTRDiffL=NULL, version="hg19", binSize = 30, mergeFile = "merged.bam") {
  # 1. Install required package
  # a. For ngs.plot
  for (pkg in c("doMC", "caTools", "utils")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      utils::install.packages(pkg, dep=T)
    }
  }

  # b. For functions in this package generally
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    utils::install.packages("BiocManager")
  }

  for (pkg in c("edgeR", "limma", "Rsubread", "BSgenome", "Rsamtools", "ShortRead")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, version = "3.8")
    }
  }

  # c. biomaRt
  if (! requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
  }

  # 2. Par the plot
  graphics::par(mfrow=c(3,3))

  # 3. Scatterplots (6 plots)
  # a. RPKM (2 plots)
  SpRPKM(inputFile, sampleFiles, version=version)

  # b. TPM (2 plots)
  SpTPM(inputFile, sampleFiles, version=version)

  # c. BPM (2 plots)
  SpTPM(inputFile, sampleFiles, version=version, binSize=binSize, mergeFile=mergeFile)

  # 4. Gene coverage (2 or 3 plots)
  # a. gene coverage of samples over annotated regions of genome
  # b. gene coverage calculated per kilobase of annotated features
  # c. gene coverage of different length of 5' UTR and 3' UTR if UTRDiffL exists
  genomeCoverage(sampleFiles, annotaedGenomes, UTRDiffL=UTRDiffL, mergeFile=mergeFile)

  # 5. Average Profile
  print(averageProfile(sampleFiles))
  print("# Please copy above code to terminal to run the ngs.plot")
}


# [END]
