# averageProfile.R
#' averageProfile to get profile from RNMT-enriched targets and not enriched targets.
#' Pre-request: download ngs.plot from https://drive.google.com/drive/folders/0B1PVLadG_dCKN1liNFY0MVM1Ulk
#'
#' \code{averageProfile} gets profile from RNMT-enriched target (50 counts per million), the aligned region could be tss, tes, genebody, exon, cgi, enhancer, dhs or bed.
#'
#' @param sampleFiles The samples of alignment files as a list in bam form
#' @param region genome region to plot, should be tss, tes, genebody, exon, cgi, enhancer, dhs or bed
#' @param mergeFile File name of merged file
#' @param mergeNGSFile File name of merged file with all start with chr, only bam file with same pattern could use ngs.plot
#' @param enrichedFile File name of genes with enriched RNMT
#' @param unenrichedFile File name of genes with unenriched RNMT
#'
#' @examples
#' \dontrun{
#' sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam",
#' "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam",
#' "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")
#' averageProfile(sampleFiles)
#' }
#' @export
averageProfile <- function(sampleFiles, region = "genebody", mergeFile = "merged.bam", mergeNGSFile = "mergedNGS.bam", enrichedFile = "enriched.txt", unenrichedFile = "unenriched.txt") {
  # Combine sample CLIP result to one merged file
  if(!file.exists(mergeFile)){
    Rsamtools::mergeBam(sampleFiles, mergeFile,overwrite=TRUE)
    Rsamtools::indexBam(mergeFile)
  }

  # 1. Get cpm for merged bam file
  fc <- Rsubread::featureCounts(files=mergeFile,annot.inbuilt="hg19")
  x <- edgeR::DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])
  x_cpm <- edgeR::cpm(x)

  rich_cpm <- as.integer(row.names(x_cpm)[which(x_cpm>50)])
  unrich_cpm <- as.integer(row.names(x_cpm)[which(x_cpm<=50)])

  # 2. Filter the rich & unrich gene
  myMart <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")

  filters <- biomaRt::listFilters(myMart) #"entrezgene"
  attributes <- biomaRt::listAttributes(myMart) # "ensembl_gene_id"

  tmp <- biomaRt::getBM(filters = "entrezgene",
    attributes = c("entrezgene","ensembl_gene_id"),
    values = rich_cpm,
    mart = myMart)

  utils::write.table(tmp$ensembl_gene_id, enrichedFile, append = FALSE, sep = "\n", dec = ".",
    row.names = FALSE, col.names = FALSE,quote=FALSE)

  tmp <- biomaRt::getBM(filters = "entrezgene",
    attributes = c("entrezgene","ensembl_gene_id"),
    values = unrich_cpm,
    mart = myMart)

  utils::write.table(tmp$ensembl_gene_id, unenrichedFile, append = FALSE, sep = "\n", dec = ".",
    row.names = FALSE, col.names = FALSE,quote=FALSE)

  # 3. Filter the bam file with 'chr'
  # samtools view -h -o out.sam merged.bam
  system(paste0("samtools view -h -o out.sam ", mergeFile))

  out <- utils::read.delim("out.sam",
                    header=FALSE,
                    sep="\n",
                    quote="",
                    stringsAsFactors = FALSE)
  b <- grepl("chr|@HD|@PG|@CO", out[,1])
  newOut <-out[b,]

  utils::write.table(newOut, "test.sam", append = FALSE, sep = "\n", dec = ".",
    row.names = FALSE, col.names = FALSE,quote=FALSE)

  # samtools view -Sb  test.sam  > mergedNGS.bam
  system(paste0("samtools view -Sb  test.sam  > ", mergeNGSFile))

  # 4. plot
  if (region %in% c("tss", "tes", "genebody", "exon", "cgi", "enhancer", "dhs", "bed")) {
    print(paste0("ngs.plot.r -G hg19 -R genebody -C ", mergeNGSFile, " -E ", enrichedFile, " -O t1 -D ensembl"))
    print(paste0("ngs.plot.r -G hg19 -R genebody -C ", mergeNGSFile, " -E ", unenrichedFile, " -O t1 -D ensembl"))
  } else {
    print("The region should be one of tss, tes, genebody, exon, cgi, enhancer, dhs or bed!")
  }

  utils::capture.output(file.remove(list.files("./", pattern = "out.sam")))
  utils::capture.output(file.remove(list.files("./", pattern = "test.sam")))
}

# [END]
