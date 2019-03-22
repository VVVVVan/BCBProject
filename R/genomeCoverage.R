# genomeCoverage.R
#' Gene coverage of samples over annotated regions of genome.
#'
#' \code{genomeCoverage} Plot the genome coverage plot for 5'UTR, 3'UTR, exons, introns, promoters. lncRNA, sno/miRNA.
#
#'
#' @param sampleFiles The samples of alignment files as a list in bam form
#' @param annotaedGenomes The list of annotated genomes for "5' UTR", "3' UTR", "exons", "introns", "promoters", "lncRNA", "tRNA", "snc / miRNA". The order of the annotated genome is important.
#' @param mergeFile default "merged.bam". Merged file of all samples to do gene coverage
#' @param UTRDiffL UTR feature file with different length in bed format, default if NULL
#'
#' @examples
#' \dontrun{
#' sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam",
#' "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam",
#' "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")
#' annotaedGenomes <- c("../BCBProjectData/5UTR2000.bed", "../BCBProjectData/3UTR2000.bed",
#' "../BCBProjectData/GENCODE.Exon.bed", "../BCBProjectData/GENCODE.Intron.bed",
#' "../BCBProjectData/GENCODE.Promoter.bed", "../BCBProjectData/lincRNA.bed",
#' "../BCBProjectData/tRNA.bed", "../BCBProjectData/miRNA.bed")
#' UTRDiffL <- c("../BCBProjectData/3UTR2000.bed","../BCBProjectData/3UTR1000.bed",
#' "../BCBProjectData/3UTR500.bed","../BCBProjectData/3UTR200.bed",
#' "../BCBProjectData/5UTR2000.bed","../BCBProjectData/5UTR1500.bed",
#' "../BCBProjectData/5UTR1000.bed","../BCBProjectData/5UTR500.bed")
#' genomeCoverage(sampleFiles, annotaedGenomes, UTRDiffL)
#' }
#' @export
genomeCoverage <- function(sampleFiles, annotaedGenomes, UTRDiffL=NULL, mergeFile = "merged.bam") {
  print("Plot % coverage and % coverage per kb")

  # 1. Combine sample CLIP result to one merged file
  if(!file.exists(mergeFile)){
    Rsamtools::mergeBam(sampleFiles, mergeFile,overwrite=TRUE)
    Rsamtools::indexBam(mergeFile)
  }

  # 2. Get the genome coverage of the bam file
  command <- paste0("bedtools genomecov -ibam ", mergeFile, " -bga -split | awk '$4!=0' > coverage.bed")
  system(command)

  # 3. Compare the coverage file to the features file
  i <- 1
  coverage <- coverageKB <- ratio <- c()
  notCoverage <- overlap <- union <- character() # For ambiguous and intergenic later

  for (file in annotaedGenomes) {
    print(paste0("Coverage bam file to", file))
    # 3.1 form coverage file of merged file to a feature file
    system(paste0("bedtools coverage -a coverage.bed -b ", file, " > output.bed"), ignore.stderr=TRUE)

    # 3.2 read the features count
    a <- utils::read.table("output.bed",
                    header = FALSE,
                    sep="\t",
                    stringsAsFactors=FALSE)

    # 3.3 read the features file
    fileDF <- utils::read.table(file,
                         header = FALSE,
                         sep="\t",
                         stringsAsFactors=FALSE)

    # There are duplicates in feature files, get rid of repeats
    fileDF[ncol(fileDF)+1] <- paste0(fileDF[,3], fileDF[,2])
    fileDF <- fileDF[(! duplicated(fileDF[,ncol(fileDF)])),]
    numFileDF <- sum(abs((fileDF[,3])-(fileDF[,2])))

    # 3.4 Store the coverage and coverge per KB value
    # col in a-> name, start, end, block count, reads, coverage of B, length of A, ratio
    coverage <- c(coverage, sum(a[,ncol(a)-2])/sum(a[,ncol(a)-1]))
    coverageKB <- c(coverageKB, (sum(a[(a[ncol(a)] != 0),5])) / (numFileDF/1000))

    # Store data for ambiguous and intergenic
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

  # 5. Add ambigious and intergenic coverage
  coverage <- c(coverage, sum(a[as.integer(overlap),][,ncol(a)-1] * a[as.integer(overlap),][,4])/sum(a[,ncol(a)-1] * a[,4]))
  coverage <- c(coverage, sum(a[as.integer(notCoverage),][,ncol(a)-1]* a[as.integer(notCoverage),][,4])/sum(a[,ncol(a)-1]* a[,4]))

  # 6. Plot
  P1 <- graphics::barplot(coverage*100,ylab = "% Coverage")
  graphics::text(P1[,1], -3.7, srt=60, adj=1, xpd=TRUE, labels = c("5' UTR", "3' UTR", "exons", "introns", "promoters", "lncRNA", "tRNA", "snc / miRNA", "ambiguous", "intergenic"))

  P2 <- graphics::barplot(coverageKB[1:4]*100, names.arg = "",ylab = "% Coverage per kb")
  graphics::text(P2[,1], -3.7, srt=60, adj=1, xpd=TRUE, labels = c("5' UTR", "3' UTR", "exons", "introns"))

  # 6. Check for different length of 5'UTR & 3'UTR if there exist UTRDiffL
  if (length(UTRDiffL)>1) {
    diffUTR <- coverageKB[3:4]
    for (file in UTRDiffL) {
      print(paste0("Coverage bam file to", file))
      system(paste0("bedtools coverage -a coverage.bed -b ", file, " > output.bed"), ignore.stderr=TRUE)
      # 6.1 read the features count
      a <- read.table("output.bed",
        header = FALSE,
        sep="\t",
        stringsAsFactors=FALSE)
      # 6.2 read the features file
      fileDF <- read.table(file,
        header = FALSE,
        sep="\t",
        stringsAsFactors=FALSE)

      # There are duplicates in feature files, get rid of repeats
      fileDF[ncol(fileDF)+1] <- paste0(fileDF[,3], fileDF[,2])
      fileDF <- fileDF[(! duplicated(fileDF[,ncol(fileDF)])),]
      numFileDF <- sum(abs((fileDF[,3])-(fileDF[,2])))

      # 6.3 Store the coverage and coverge per KB value
      diffUTR <- c((sum(a[(a[ncol(a)] != 0),5])) / (numFileDF/1000), diffUTR)
    }

    P3 <- graphics::barplot(diffUTR*100, names.arg = "",ylab = "% Coverage per kb", ylim=c(0,max(diffUTR*100)+10))
    graphics::text(P3[,1], -3.7, srt=60, adj=1, xpd=TRUE, labels = c("5UTR500","5UTR1000","5UTR1500", "5UTR2000","3UTR200","3UTR500","3UTR1000","3UTR2000", "exons", "introns"))
  }

  # Remove the result file if necessary
  # file.remove(list.files("./", pattern = "merged.bam"))
  # file.remove(list.files("./", pattern = "coverage.bed"))
  # file.remove(list.files("./", pattern = "output.bed"))
}

#[END]
