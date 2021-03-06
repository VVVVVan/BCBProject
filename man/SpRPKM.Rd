% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpRPKM.R
\name{SpRPKM}
\alias{SpRPKM}
\title{Scatterplot of RPKM of inputs vs. average RPKM of sample.}
\usage{
SpRPKM(inputFile, sampleFiles, version = "hg19")
}
\arguments{
\item{inputFile}{The input of alignment files as a list in bam form}

\item{sampleFiles}{The samples of alignment files as a list in bam form}

\item{version}{default is "hg19". Could be "mm10", "mm9", "hg38" and "hg19"}
}
\description{
\code{SpRPKM} Scatterplot of RPKM of inputs vs. average RPKM of sample. Count reads per trancript and get the RPKM; plot the scatterplot for average RPKM samples and RPKM inputs.
}
\examples{
\dontrun{
inputFile <- "../BCBProjectData/SRR4393137.Aligned.sortedByCoord.out.bam"
sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam",
"../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam",
"../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")
SpRPKM(inputFile, sampleFiles)
}
}
