# getTSS.R
#' Get TSS (transcript start site)
#'
#' \code{getTSS} Find TSS sites and store in bed file (tab seperated)
# TSS: Find the 5' UTR and get the next one
#      chrom - chromosome name
#      start - start site of a TSS
#      end - end site of a TSS + 1
#      name - name of the gene
#      score - use "period" as name holder
#      strand - defined as + (forward) or - (reverse)
#'
#' @param fileName the name of the file a gtf file
#' @param outputFile the name of output file, default is TSS.bed
#'
#' @seealso \code{\link[rtracklayer]{import}} is used to import the gtf file.
#'
#' @examples
#' \dontrun{
#'  getTSS("../BCBProjectData/gencode.v19.annotation.gtf")
#' }
#' @export
getTSS <- function(fileName, outputFile = "TSS.bed") {
  # fileName <- gencode.v19.annotation.gtf"
  gtfFile <- system.file("extdata", fileName, package = "BCBProject")
  gtf <- rtracklayer::import(gtfFile)
  gtf_df <- as.data.frame(gtf)

  # save(list="gtf_df", file="./inst/extdata/gtfDataFrame.Rdata")
  # load(system.file("extdata", "gtfDataFrame.Rdata", package = "BCBProject"))

  UTRsequences <- gtf_df[which(gtf_df["type"] == "UTR"), ]
  # Initial a int for numebr of 5' UTRs and output data frame
  accumutaiveUTRs <- as.integer(1)
  output <- data.frame(chrom = character(),
                       start = integer(),
                       end = integer(),
                       name = character(),
                       score = character(),
                       strand = character(),
                       stringsAsFactors = FALSE)

  for (i in 1:nrow(UTRsequences)) {
    if (UTRsequences[i, "strand"] == "+") {
      # Check if the sequence after the end of UTR is CDS start site
      item <- UTRsequences[i,]
      startCDS <- as.integer(item["end"] + 1)
      existCDSs <- which(gtf_df["start"] == startCDS)
      # print(i) # To check where the loop is
      for (CDS in existCDSs) {
        if (gtf_df[CDS, "type"] == "CDS") {
          # add the requried things
          x <- data.frame(chrom  = item["seqnames"],
                          start  = as.integer(item["start"]),
                          end    = as.integer(item["start"] + 1),
                          name   = item["gene_name"],
                          score  = "period",
                          strand = item["strand"],
                          stringsAsFactors = FALSE)
          #if (! (nrow(merge(x, output)) > 0)) {
            row.names(x) <- accumutaiveUTRs
            output <- rbind(output, x)
            rm(x)
            accumutaiveUTRs <- accumutaiveUTRs + 1
          #}
        }
      }
    } else if (UTRsequences[i, "strand"] == "-") {
      # Check if the sequecne before the start of UTR is CDS end site
      item <- UTRsequences[i,]
      startCDS <- as.integer(item["start"] - 1)
      existCDSs <- which(gtf_df["end"] == startCDS)
      for (CDS in existCDSs) {
        if (gtf_df[CDS, "type"] == "CDS") {
          # add the requried things
          y <- data.frame(chrom  = item["seqnames"],
                          start  = as.integer(item["end"] - 1),
                          end    = as.integer(item["end"]),
                          name   = item["gene_name"],
                          score  = "period",
                          strand = item["strand"],
                          stringsAsFactors = FALSE)
          #if (! (nrow(merge(y, output)) > 0)) {
            row.names(y) <- accumutaiveUTRs
            output <- rbind(output, y)
            rm(y)
            accumutaiveUTRs <- accumutaiveUTRs + 1
          #}
        }
      }
    }
  }

  # save(list="output", file="./inst/extdata/TSSDataFrame.Rdata")
  # load(system.file("extdata", "TSSDataFrame.Rdata", package = "BCBProject"))

  TSSexport <- unique(output) # get the unique rows
  TSSexport[,"score"] <- 0.0
  colnames(TSSexport) <- c("chrom", "start", "end", "name", "score", "strand")
  rtracklayer::export(TSSexport, "TSSexport.bed")

  # Modifiy the output bed file by change the start -1 to left
  system(paste0("bedtools slop -i TSSexport.bed -g ../BCBProjectData/human.hg19.genome -l -1 -r 0 > ", outputFile))
  file.remove(list.files("./", pattern = "TSSexport.bed"))
}

#' Expand TSS (transcript start site)
#'
#' \code{expandTSS} Expand the TSS to 500
#'
#' @param outputFile the name of output file, default is TSSextend.bed
#'
#' @examples
#' \dontrun{
#'  expandTSS()
#' }
#' @export
expandTSS <- function(outputFile) {
  system(paste0("bedtools slop -i TSS.bed -g ../BCBProjectData/human.hg19.genome -l 0 -r 499 -s > ", outputFile))
}

# [end]
