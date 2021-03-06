% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpBPM.R
\name{SpBPM}
\alias{SpBPM}
\title{Scatterplot of BPM of inputs vs. average BPM of sample.}
\usage{
SpBPM(inputFile, sampleFiles, version = "hg19", binSize = 30,
  mergeFile = "merged.bam")
}
\arguments{
\item{inputFile}{The input of alignment files as a list in bam form}

\item{sampleFiles}{The samples of alignment files as a list in bam form}

\item{version}{default is "hg19". Could be "mm10", "mm9", "hg38" and "hg19"}

\item{binSize}{default is 30. The bin size for each bin}

\item{mergeFile}{default "merged.bam". Merged file of all samples to do bamCompare function}
}
\description{
\code{SpBPM} Scatterplot of BPM of inputs vs. average BPM of sample. Get the BPM from deeptools; plot the scatterplot for average BPM samples and BPM inputs.
}
\examples{
\dontrun{
inputFile <- "../BCBProjectData/SRR4393137.Aligned.sortedByCoord.out.bam"
sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam",
"../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam",
"../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")
SpBPM(inputFile, sampleFiles)
}
}
