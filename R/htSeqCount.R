# htSeqCount.R
#' Count reads per trancript with HTSeq package
#'
#' \code{htSeqCount} Count reads per transcript
#'
#' @param bamFiles The alignment files as a list in bam form
#' @param gffFile the gff file or gtf file
#' @param outputFile the output file name, default is count.csv
#'
#' @examples
#' \dontrun{
#' files <- list("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")
#' htSeqCount(files,"../BCBProjectData/gencode.v19.annotation.gtf")
#' }
#' @export
htSeqCount <- function(bamFiles, gffFile, outputFile = "count.csv") {
  command <- "python3 -m HTSeq.scripts.count -f bam "
  for (file in bamFiles) {
    command <- paste0(command, file, " ")
  }
  command <- paste0(command, gffFile, " > ", outputFile)

  system(command)
}

# [END]
