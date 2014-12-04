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
#'
#' @format A GRanges object.
#'
#' @source Processing reads from YFM
"parp1_ln4_reads_all"

#' All the PARP1 reads from LN5
#'
#' A dataset containing all the PARP1 reads from LN4
#'
#' @format A GRanges object.
#'
#' @source Processing reads from YFM
"parp1_ln5_reads_all"


#' The unique location PARP1 reads from LN4
#'
#' A dataset containing only the unique location PARP1 reads from LN4, with the number of reads
#' at each location capped at 6.
#'
#' @format A GRanges object with an mcols DataFrame with a column named n_count
#'
#' @source Processing reads from YFM
"parp1_ln4_unique"

#' The unique location PARP1 reads from LN5
#'
#' A dataset containing only the unique location PARP1 reads from LN5, with the number of reads
#' at each location capped at 6.
#'
#' @format A GRanges object with an mcols DataFrame with a column named n_count
#'
#' @source Processing reads from YFM
"parp1_ln5_unique"


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
#' @source http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeUwTfbs
"ctcf_rep1"


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
#' @source http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeUwTfbs
"ctcf_rep2"


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
#' @source http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=396587911_8YHquTEUSQjmIfAHtJJT7vWY7N8U&clade=mammal&org=Human&db=hg19&hgta_group=allTracks&hgta_track=wgEncodeHaibMethylRrbs&hgta_table=wgEncodeHaibMethylRrbsMcf7DukeSitesRep1
"methyl_rep1"

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
#' @source http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=396587911_8YHquTEUSQjmIfAHtJJT7vWY7N8U&clade=mammal&org=Human&db=hg19&hgta_group=allTracks&hgta_track=wgEncodeHaibMethylRrbs&hgta_table=wgEncodeHaibMethylRrbsMcf7DukeSitesRep2
"methyl_rep2"