[![Build Status](https://travis-ci.org/rmflight/fmdatabreastcaparp1.svg?branch=master)](https://travis-ci.org/rmflight/fmdatabreastcaparp1)

# Fondufe-Mittendorf Breast Cancer PARP1 Data Package

This package contains the data used for the correlation analysis in the publication:

**Genome-wide profiling of PARP1 reveals an interplay with gene regulatory regions and DNA methylation**, Nalabothula Narasimharao, Taha Al-jumaily, Robert M. Flight, Shao Xiaorong, Hunter N. B. Moseley, F. Lisa Barcellos, Yvonne Fondufe-Mittendorf, Submitted

This data is used by the [correlation analysis package](https://github.com/rmflight/fmcorrelationbreastcaparp1) to generate the vignettes that report the actual correlations between PARP1 reads and other data. Unless otherwise noted, all data was downloaded for the MCF-7 cell line.

Please see the vignette for details on sources of data and processing. In most cases, what is done is simply to take data from UCSC and generate a `GenomicRanges` object that can be worked with in the analysis.

This package contains the following `RData` objects:

* `parp1_ln4_reads_all`: all PARP1 reads from MCF-7
* `parp1_ln5_reads_all`: all PARP1 reads from MDA-MB231
* `parp1_ln4_unique`: summarized reads, with counts at locations capped at 6
* `parp1_ln5_unique`: summarized reads, with counts at locations capped at 6
* `ctcf_rep1`: CTCF ChIP-Seq peaks, replicate 1
* `ctcf_rep2`: CTCF ChIP-Seq peaks, replicate 2
* `methyl_rep1`: methylation reads, replicate 1
* `methyl_rep2`: methylation reads, replicate 2
* `histone_marks`: `GRangesList` of ChIP-Seq peaks of various histone marks
* `expr_data`: Affymetrix HGU-133 2.0 expression data
* `tss_windows`: 2KB windows around Ensembl transcrip start sites
* `nuc_signal`: summed MNase nucleosome signal over 500 base tiles of the genome in GM12878

## Citation

If this package is used in other analyses, it should be cited as:

**fmdatabreastcaparp1: Fondufe-Mittendorfe lab data for PARP1 correlation in MCF-7 and MDA-MB231**, R. M. Flight, T. Al-jumaily, Y. Fondufe-Mittendorf, H. N. B. Moseley  doi:10.6084/m9.figshare.1267545

## Installation

To install this package, you will need to first download the data files from `figshare`, and then clone this package and install:

```
git clone https://github.com/rmflight/fmdatabreastcaparp1.git
mkdir fmdatabreastcaparp1/data
wget https://dl.dropbox.com/s/89xgvje4rfl0lml/fmdatabreastcaparp1.zip
unzip fmdatabreastcaparp1.zip -d fmdatabreastcaparp1/data

R
library(devtools)
install("fmdatabreastcaparp1")
```

## Notes on Installation and Running Vignette

This package requires a fair amount of memory for the vignette to be run. As far as I can tell, it will need ~4GB to run the part of the vignette that determines what to cap the reads at.

Additionally, I have not been able to do more than a `build` and `install` on Travis-ci. However, the package does pass `R CMD check` on my own machine.
