#' get chromosome
#' 
#' from a filename, get the chromosome part of the filename
#' 
#' @param filename the filename to parse
#' @param split_chr the character separating the chr* from the rest
#' @return string of the chromosome
#' @export
get_chr <- function(filename, split_chr){
  file_split <- strsplit(filename, split_chr)[[1]]
  n_split <- length(file_split)
  chr_part <- file_split[n_split]
  substr(chr_part, 1, nchar(chr_part)-4)
}


#' All the PARP1 reads from LN4
#'
#' A dataset containing all the PARP1 reads from LN4
#' @name parp1_ln4_reads_all
#' @format A GRanges object.
#'
#' @source Processing reads from YFM
NULL

#' All the PARP1 reads from LN5
#'
#' A dataset containing all the PARP1 reads from LN4
#'
#' @format A GRanges object.
#' @name parp1_ln5_reads_all
#' @source Processing reads from YFM
NULL


#' The unique location PARP1 reads from LN4
#'
#' A dataset containing only the unique location PARP1 reads from LN4, with the number of reads
#' at each location capped at 6.
#'
#' @format A GRanges object with an mcols DataFrame with a column named n_count
#'
#' @source Processing reads from YFM
#' @name parp1_ln4_unique
NULL

#' The unique location PARP1 reads from LN5
#'
#' A dataset containing only the unique location PARP1 reads from LN5, with the number of reads
#' at each location capped at 6.
#'
#' @format A GRanges object with an mcols DataFrame with a column named n_count
#'
#' @source Processing reads from YFM
#' @name parp1_ln5_unique
NULL


#' CTCF rep1
#'
#' CTCF ChIP-Seq peak intensity.
#'
#' @format A GRanges object with an mcols DataFrame with columns:
#' \describe{
#'   \item{mcols.signal}{signal}
#'   \item{mcols.pValue}{the pValue of the peak}
#' }
#'
#' @source \url{http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeUwTfbs}
#' @name ctcf_rep1
NULL


#' CTCF rep2
#'
#' CTCF ChIP-Seq peak intensity.
#'
#' @format A GRanges object with an mcols DataFrame with columns:
#' \describe{
#'   \item{mcols.signal}{signal}
#'   \item{mcols.pValue}{the pValue of the peak}
#' }
#'
#' @source \url{http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeUwTfbs}
#' @name ctcf_rep2
NULL


#' Methylation rep1
#'
#' Methylation reads
#'
#' @format A GRanges object with an mcols DataFrame with columns:
#' \describe{
#'   \item{mcols.readCount}{the number of reads}
#'   \item{mcols.percentMeth}{the percentage of reads that were methylated}
#' }
#'
#' @source \url{http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=396587911_8YHquTEUSQjmIfAHtJJT7vWY7N8U&clade=mammal&org=Human&db=hg19&hgta_group=allTracks&hgta_track=wgEncodeHaibMethylRrbs&hgta_table=wgEncodeHaibMethylRrbsMcf7DukeSitesRep1}
#' @name methyl_rep1
NULL

#' Methylation rep2
#'
#' Methylation reads
#'
#' @format A GRanges object with an mcols DataFrame with columns:
#' \describe{
#'   \item{mcols.readCount}{the number of reads}
#'   \item{mcols.percentMeth}{the percentage of reads that were methylated}
#' }
#'
#' @source \url{http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=396587911_8YHquTEUSQjmIfAHtJJT7vWY7N8U&clade=mammal&org=Human&db=hg19&hgta_group=allTracks&hgta_track=wgEncodeHaibMethylRrbs&hgta_table=wgEncodeHaibMethylRrbsMcf7DukeSitesRep2}
#' @name methyl_rep2
NULL


#' Transcript start sites
#'
#' 2KB windows around Ensembl transcript start sites (TSS)
#'
#' @format A named GRanges object, where the names are the Ensembl transcript IDs
#'
#' @source ensembl transcripts
#' @name tss_windows
NULL

#' Histone Marks
#' 
#' ChIP-Seq peaks of various histone marks
#' 
#' @format GRangesList of length 6, with H3k09me3, H3k27ac, H3k27me3, H3k36me3,
#' H3k4me3_r1 and H3k4me3_r2. Each is a GRanges object with mcols of:
#' \describe{
#'   \item{mcols.score}{the score}
#'   \item{mcols.signal}{the peak intensity}
#'   \item{mcols.pvalue}{-1 * log10 of the pvalue}
#' }
#' 
#' @name histone_marks
#' @source \url{http://genome.ucsc.edu/cgi-bin/hgTracks?hgsid=339462435&hgt_=1371663888&db=hg19&tsCurTab=advancedTab&hgt_tsDelRow=&hgt_tsAddRow=&hgt_tsPage=&tsSimple=&tsName=histone&tsDescr=&tsGroup=Any&tsType=Any&hgt_mdbVar1=cell&hgt_mdbVal1=MCF-7&hgt_mdbVar2=antibody&hgt_mdbVal2=Any&hgt_tSearch=search}
NULL

#' Expression Data
#' 
#' DNA Microarray Expression Data from GEO, Affymetrix hgu133plus2, GSM307014
#' measured on MCF-7 cells
#' 
#' @format data.frame with:
#' \describe{
#'   \item{ID_REF}{the affymetrix probe set id}
#'   \item{VALUE}{the normalized signale}
#'   \item{ABS_CALL}{absence presence call}
#'   \item{P-VALUE}{p-value of being present}
#' }
#' @name expr_data
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM307014}
NULL
#'   