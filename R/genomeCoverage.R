# genomeCoverage.R
#' Genome coverage of bam files
#'
#' \code{genomeCoverage} Plot the genome coverage plot for 5'UTR, 3'UTR, exons, introns, promoters. lncRNA, sno/miRNA, ambiguous(overlapping annotation), intergenic(not map to annotated regions).
#
#'
#' @param bamFiles the name of the file a gtf file
#' @param annotaedGenomes The list of annotated genomes for "5' UTR", "3' UTR", "exons", "introns", "promoters", "lncRNA", "tRNA", "snc / miRNA". The order of the annotated genome is important.
#'
#'
#' @examples
#' \dontrun{
#' bamFiles <- list("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")
#' annotaedGenomes <- list("../BCBProjectData/GENCODE.5UTR.bed", "../BCBProjectData/GENCODE.3UTR.bed", "../BCBProjectData/GENCODE.Exon.bed", "../BCBProjectData/GENCODE.Intron.bed", "../BCBProjectData/GENCODE.Promoter.bed", "../BCBProjectData/lincRNA.bed", "../BCBProjectData/tRNA.bed", "../BCBProjectData/miRNA.bed")
#' genomeCoverage(bamfiles,annotaedGenomes)
#' }
#' @export
genomeCoverage <- function(bamFiles, annotaedGenomes) {
  # Combine three CLIP result to one file
  command <- "samtools merge merged.bam "
  for (file in bamFiles) {
    command <- paste0(command, file, " ")
  }
  system(command)

  # Get the genome coverage of the bam file
  system("bedtools genomecov -ibam merged.bam -bga -split | awk '$4!=0' > coverage.bed")

  i <- 1
  coverage <- c()
  coverageKB <- c()
  name <- c("5' UTR", "3' UTR", "exons", "introns", "promoters", "lncRNA", "tRNA", "snc / miRNA")
  notCoverage <- character()
  overlap <- character()
  union <- character()
  file <- "../BCBProjectData/GENCODE.3UTR.bed"
  for (file in annotaedGenomes) {
    system(paste0("bedtools coverage -a coverage.bed -b ", file, " > output.bed"))
    a <- read.table("output.bed",
                    header = FALSE,
                    sep="\t",
                    stringsAsFactors=FALSE)
    fileDF <- read.table(file,
      header = FALSE,
      sep="\t",
      stringsAsFactors=FALSE)

    numUTR <- sum(abs((fileDF[,3])-(fileDF[,2]))) # There are duplicates, TODO: need to get rid of repeats
    fileDF[7] <- paste0(fileDF[,3], fileDF[,2])
    fileDF <- fileDF[(! duplicated(fileDF[,7])),]
    numFileDF <- sum(abs((fileDF[,3])-(fileDF[,2])))

    coverage <- c(coverage, sum(a[,ncol(a)-2] * a[,4])/sum(a[,ncol(a)-1] * a[,4]))
    coverageKB <- c(coverageKB, (sum(a[(a[ncol(a)] != 0),4])) / (numFileDF/1000))

    if (i == 1) {
      notCoverage <- rownames(a[(a[ncol(a)] == 0),])
      union <- rownames(a[(a[ncol(a)] != 0),])
    } else {
      notCoverage <- intersect(notCoverage, rownames(a[(a[ncol(a)] == 0),]))
      inter <- intersect(union, rownames(a[(a[ncol(a)] != 0),]))
      overlap <- union(overlap, inter)
      union <- union(union, rownames(a[(a[ncol(a)] != 0),]))
    }
    i <- i+1
  }
  # Add ambigious and intergenic coverage
  coverage <- c(coverage, sum(a[as.integer(overlap),][,ncol(a)-1] * a[as.integer(overlap),][,4])/sum(a[,ncol(a)-1] * a[,4]))
  coverage <- c(coverage, sum(a[as.integer(notCoverage),][,ncol(a)-1]* a[as.integer(notCoverage),][,4])/sum(a[,ncol(a)-1]* a[,4]))

  # Plot
  p = barplot(coverage*100,names.arg = "", ylab = "% Coverage", c(0, max(coverage*100)))
  text(p[,1], -3.7, srt=60, adj=1, xpd=TRUE, labels = c("5' UTR", "3' UTR", "exons", "introns", "promoters", "lncRNA", "tRNA", "snc / miRNA", "ambiguous", "intergenic"))

  P2 = barplot(coverageKB[1:4], names.arg = "",ylab = "% Coverage per kb") # TODO: not look like the paper
  text(p[,1], -3.7, srt=60, adj=1, xpd=TRUE, labels = c("5' UTR", "3' UTR", "exons", "introns"))

  # Remove the result file if necessary
  # file.remove(list.files("./", pattern = "merged.bam"))
  # file.remove(list.files("./", pattern = "coverage.bed"))
  # file.remove(list.files("./", pattern = "output.bed"))
}


#[END]
