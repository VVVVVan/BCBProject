# BCBProject


# Prerequests
## Download database
Download the gtf file and genome file of human gene v19, i.e gencode.v19.annotation.gtf, human.hg19.genome.
Download the 5' UTR, 3'UTR, exons, introns, promoters (1.5kb upstream), lncRNA, tRNA, sno/miRNA from (table brower of UCSC)[http://genome.ucsc.edu/cgi-bin/hgTables], (with GENCODE gene v19 for first five).

## Install pakcages
Install bedtools
Install samtools
Install (HTSeq)[https://htseq.readthedocs.io/en/release_0.11.1/install.html].
  * Only count is used, no need of install matplotlib.

# 1. Get TSS and extend

This function `getTSS()` is going to get the transcript start site from human gene and result a bed file i.e. TSS.bed. Also, `expandTSS` is going to expand the TSS to 500 and result a bed file i.e. TSSextend.bed.

```R
getTSS("../BCBProjectData/gencode.v19.annotation.gtf")
expandTSS()
```

# 2. 


# Reference
S Anders, T P Pyl, W Huber: HTSeq — A Python framework to work with high-throughput sequencing data. bioRxiv 2014. (doi: 10.1101/002824)[https://www.biorxiv.org/content/10.1101/002824v2].
