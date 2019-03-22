# `BCBProject`

&nbsp;

###### [Fan Shen](https://orcid.org/0000-0001-8720-5874), Undergraduate of Bioinformatics and Computaional Biology, University of Toronto, Canada. &lt;van.shen@mail.utoronto.ca&gt;

----

## 1. Description
BCB430 2018 Fall - 2019 Winter.
This pacakge is going to analysis the bam file after alignment.

Supervisor: Zhaolei Zhang, Hyunmin Lee

&nbsp;

# 2. Prerequests
## Download database
a. Download the gtf file and genome file of human gene v19, i.e gencode.v19.annotation.gtf, human.hg19.genome.

b. Download the 5' UTR, 3'UTR, exons, introns, promoters (1.5kb upstream), lncRNA, tRNA, sno/miRNA from (table brower of UCSC)[http://genome.ucsc.edu/cgi-bin/hgTables], with GENCODE gene v19 for first five.

* Note: Please selcet proper version for your data.

&nbsp;

## Download tools
a. Install (deeptools)[https://deeptools.readthedocs.io/en/develop/content/installation.html].
b. Install (ngs.plot)[https://github.com/shenlab-sinai/ngsplot]

* Note: If deeptools cannot call from R s.t meet error like "sh: bamCompare: command not found", call this:
```R
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/the/bin/folder/of/deeptools", sep=":"))
```

&nbsp;

# 3. Generate several plots to analysis the bam files. 
## 1 Scatterplot 
Scatterplot of RPKM, TPM and BPM of inputs vs. average RPKM of samples.

&nbsp;

### a. Scatterplot of RPKM of inputs vs. average RPKM of samples; 
Example code:
```R
inputFile <- "../BCBProjectData/SRR4393137.Aligned.sortedByCoord.out.bam"
sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")

SpRPKM(inputFile, sampleFiles)
```
Example outcomes:

![](./inst/img/RPKM.jpeg?sanitize=true "Scatterplot of RPKM")

&nbsp;

### b. Scatterplot of TPM of inputs vs. average TPM of samples; 
Example code:
```R
inputFile <- "../BCBProjectData/SRR4393137.Aligned.sortedByCoord.out.bam"
sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")

SpTPM(inputFile, sampleFiles)
```
Example outcomes:

![](./inst/img/TPM.jpeg?sanitize=true "Scatterplot of TPM")

&nbsp;

### c. Scatterplot of BPM of inputs vs. average BPM of samples with binsize 30. 
Example code:
```R
inputFile <- "../BCBProjectData/SRR4393137.Aligned.sortedByCoord.out.bam"
sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")

# This would took a while
SpBPM(inputFile, sampleFiles)
```
Example outcomes:

![](./inst/img/BPM_nonzero.jpeg?sanitize=true "Scatterplot of BPM from multiBamSummary")
![](./inst/img/bamCompare.jpeg?sanitize=true "Scatterplot of BPM from bamCompare")

&nbsp;

## 2 Gene coverage
Gene coverage of samples over 
a. annotated regions of genome (default is hg19); 
b. calculated per kilobase of annotated features; 
c. different length of 5' UTR and 3' UTR. 

&nbsp;

Example code:
```R
sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")

annotaedGenomes <- c("../BCBProjectData/5UTR2000.bed", "../BCBProjectData/3UTR2000.bed", "../BCBProjectData/GENCODE.Exon.bed", "../BCBProjectData/GENCODE.Intron.bed", "../BCBProjectData/GENCODE.Promoter.bed", "../BCBProjectData/lincRNA.bed", "../BCBProjectData/tRNA.bed", "../BCBProjectData/miRNA.bed")

UTRDiffL <- c("../BCBProjectData/3UTR2000.bed","../BCBProjectData/3UTR1000.bed","../BCBProjectData/3UTR500.bed","../BCBProjectData/3UTR200.bed","../BCBProjectData/5UTR2000.bed","../BCBProjectData/5UTR1500.bed","../BCBProjectData/5UTR1000.bed","../BCBProjectData/5UTR500.bed")

genomeCoverage(sampleFiles, annotaedGenomes, UTRDiffL=UTRDiffL)
```
Example outcomes:

![](./inst/img/genomeCoverage.jpeg?sanitize=true "Genome Coverage")

&nbsp;

## 3 Average profile
Average profile of samples over enriched and not-enriched genes by ng.plot.

&nbsp;

Example code:
```R
sampleFiles <- c("../BCBProjectData/SRR4393138.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393139.Aligned.sortedByCoord.out.bam", "../BCBProjectData/SRR4393140.Aligned.sortedByCoord.out.bam")

averageProfile(sampleFiles)
# Would output some command like this. Please copy below code to terminal to run the ngs.plot.
# [1] "ngs.plot.r -G hg19 -R genebody -C mergedNGS.bam -E enriched.txt -O t1 -D ensembl"             
# [1] "ngs.plot.r -G hg19 -R genebody -C mergedNGS.bam -E unenriched.txt -O t1 -D ensembl"

```
Example outcomes:

![](./inst/img/t1.avgprof.pdf?sanitize=true "ngs.plot enriched")
![](./inst/img/t2.avgprof.pdf?sanitize=true "ngs.plot unenriched")

&nbsp;

# 4. Get TSS and extend

This function `getTSS()` is going to get the transcript start site from human gene and result a bed file i.e. TSS.bed. Also, `expandTSS` is going to expand the TSS to 500 and result a bed file i.e. TSSextend.bed.

```R
getTSS("../BCBProjectData/gencode.v19.annotation.gtf")
expandTSS()
```

&nbsp;

# 5. Reference
Ramírez, Fidel, Devon P. Ryan, Björn Grüning, Vivek Bhardwaj, Fabian Kilpert, Andreas S. Richter, Steffen Heyne, Friederike Dündar, and Thomas Manke. deepTools2: A next Generation Web Server for Deep-Sequencing Data Analysis. Nucleic Acids Research (2016). doi:10.1093/nar/gkw257.

Shen, L.*, Shao, N., Liu, X. and Nestler, E. (2014) ngs.plot: Quick mining and visualization of next-generation sequencing data by integrating genomic databases, BMC Genomics, 15, 284.

Varshney, D., Lombardi, O., Schweikert, G., Dunn, S., Suska, O., & Cowling, V. H. (2018). mRNA Cap Methyltransferase, RNMT-RAM, Promotes RNA Pol II-Dependent Transcription. Cell reports, 23(5), 1530-1542.

