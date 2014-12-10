## ----setup---------------------------------------------------------------
library(fmdatabreastcaparp1)
library(GenomicRanges)
library(GEOquery)
raw_data_dir <- "/mlab/data/rmflight/Documents/projects/work/fondufe-mittendorf_lab/parp1_data/"
save_dir <- "/mlab/data/rmflight/Documents/projects/work/fondufe-mittendorf_lab/fmdatamcf7parp1/data"

## ----ln4_reads, eval=FALSE-----------------------------------------------
#  ln4_files <- file.path(raw_data_dir, dir(raw_data_dir, pattern = "YFM_LN4"))
#  
#  parp1_ln4_reads_all <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1, end = 1), strand = "*")
#  
#  for (i_file in ln4_files){
#    tmp_reads <- read.table(i_file, sep = ",", header = TRUE)
#  	use_chr <- get_chr(i_file, "_")
#  	parp1_ln4_reads_all <- c(parp1_ln4_reads_all, GRanges(seqnames = use_chr,
#  					  ranges = IRanges(start = tmp_reads[, "startx"], width = 1),
#  					  strand = "*"))
#  }
#  parp1_ln4_reads_all <- parp1_ln4_reads_all[2:(length(parp1_ln4_reads_all))]
#  save(parp1_ln4_reads_all, file = file.path(save_dir, "parp1_ln4_reads_all.RData"))
#  rm(parp1_ln4_reads_all)

## ----ln5_reads, eval=FALSE-----------------------------------------------
#  ln5_files <- file.path(raw_data_dir, dir(raw_data_dir, pattern = "YFM_LN5"))
#  
#  parp1_ln5_reads_all <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1, end = 1), strand = "*")
#  
#  for (i_file in ln5_files){
#    tmp_reads <- read.table(i_file, sep = ",", header = TRUE)
#    use_chr <- get_chr(i_file, "_")
#  	parp1_ln5_reads_all <- c(parp1_ln5_reads_all, GRanges(seqnames = use_chr,
#  					  ranges = IRanges(start = tmp_reads[, "startx"], width = 1),
#  					  strand = "*"))
#  }
#  parp1_ln5_reads_all <- parp1_ln5_reads_all[2:(length(parp1_ln5_reads_all))]
#  save(parp1_ln5_reads_all, file = file.path(save_dir, "parp1_ln5_reads_all.RData"))
#  rm(parp1_ln5_reads_all)

## ----load_chrom_3--------------------------------------------------------
data(parp1_ln4_reads_all)
chr3_reads <- parp1_ln4_reads_all[seqnames(parp1_ln4_reads_all) == "chr3"]
chr3_unique <- unique(chr3_reads)
chr3_overlap <- countOverlaps(chr3_unique, chr3_reads)
mcols(chr3_unique)$overlap <- chr3_overlap
chr3_unique <- sort(chr3_unique)

## ----check_chr3_counts---------------------------------------------------
chr3_overlap_rle <- rle(sort(chr3_overlap))
chr3_counts <- chr3_overlap_rle$values
names(chr3_counts) <- chr3_overlap_rle$lengths
chr3_counts

## ----chr3_99-------------------------------------------------------------
chr3_totals <- data.frame(n_loc = chr3_overlap_rle$lengths, count_at_loc = chr3_overlap_rle$values)
chr3_totals$tot_reads <- chr3_totals$n_loc * chr3_totals$count_at_loc
chr3_totals$cum_reads <- cumsum(chr3_totals$tot_reads)
chr3_totals$perc_reads <- chr3_totals$cum_reads / sum(chr3_totals$tot_reads) * 100
head(chr3_totals)

## ----ln4_all_unique------------------------------------------------------
ln4_unique <- unique(parp1_ln4_reads_all)
ln4_overlap <- countOverlaps(ln4_unique, parp1_ln4_reads_all)
ln4_overlap_rle <- rle(sort(ln4_overlap))
ln4_counts <- data.frame(n_loc = ln4_overlap_rle$lengths, count_at_loc = ln4_overlap_rle$values)

ln4_counts$tot_reads <- ln4_counts$n_loc * ln4_counts$count_at_loc
ln4_counts$cum_reads <- cumsum(ln4_counts$tot_reads)
ln4_counts$perc_reads <- ln4_counts$cum_reads / sum(ln4_counts$tot_reads) * 100
head(ln4_counts)

## ----parp1_ln4_unique, eval=FALSE----------------------------------------
#  data(parp1_ln4_reads_all)
#  max_count <- 6
#  parp1_ln4_unique <- unique(parp1_ln4_reads_all)
#  parp1_ln4_counts <- countOverlaps(parp1_ln4_unique, parp1_ln4_reads_all)
#  parp1_ln4_counts[parp1_ln4_counts > max_count] <- max_count
#  mcols(parp1_ln4_unique)$n_count <- parp1_ln4_counts
#  save(parp1_ln4_unique, file = file.path(save_dir, "parp1_ln4_unique.RData"))

## ----parp1_ln5_unique, eval=FALSE----------------------------------------
#  data(parp1_ln5_reads_all)
#  parp1_ln5_unique <- unique(parp1_ln5_reads_all)
#  parp1_ln5_counts <- countOverlaps(parp1_ln5_unique, parp1_ln5_reads_all)
#  parp1_ln5_counts[parp1_ln5_counts > max_count] <- max_count
#  mcols(parp1_ln5_unique)$n_count <- parp1_ln5_counts
#  save(parp1_ln5_unique, file = file.path(save_dir, "parp1_ln5_unique.RData"))

## ----ctcf_data, eval=FALSE-----------------------------------------------
#  ctcf_names <- c("chrom", "start", "end", "name", "score", "strand", "signal", "pValue", "qValue")
#  ctcf_path <- file.path(raw_data_dir, "ctcf_peaks")
#  ctcf_files <- file.path(ctcf_path, dir(ctcf_path, pattern = "narrowPeak.gz"))
#  rep1 <- read.table(ctcf_files[1], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#  
#  names(rep1) <- ctcf_names
#  ctcf_rep1 <- GRanges(seqnames = rep1$chrom,
#                       ranges = IRanges(start = rep1$start, end = rep1$end),
#                       mcols = DataFrame(rep1[, c("signal", "pValue")]))
#  
#  
#  rep2 <- read.table(ctcf_files[2], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#  names(rep2) <- ctcf_names
#  ctcf_rep2 <- GRanges(seqnames = rep2$chrom,
#                       ranges = IRanges(start = rep2$start, end = rep2$end),
#                       mcols = DataFrame(rep2[, c("signal", "pValue")]))
#  
#  save(ctcf_rep1, file = file.path(save_dir, "ctcf_rep1.RData"))
#  save(ctcf_rep2, file = file.path(save_dir, "ctcf_rep2.RData"))
#  rm(ctcf_rep1, ctcf_rep2)

## ----save_methylation----------------------------------------------------
methyl_names <- c("bin", "chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRGB", "readCount", "percentMeth")
methyl_path <- file.path(raw_data_dir, "methylation_data")
methyl_files <- file.path(methyl_path, dir(methyl_path, pattern = "mcf7_methyl"))

rep1 <- read.table(methyl_files[1], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(rep1) <- methyl_names
methyl_rep1 <- GRanges(seqnames = rep1$chrom,
                       strand = rep1$strand,
                      ranges = IRanges(start = rep1$start, width = 1),
                      mcols = DataFrame(rep1[, c("readCount", "percentMeth")]))

rep2 <- read.table(methyl_files[2], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
names(rep2) <- methyl_names
methyl_rep2 <- GRanges(seqnames = rep2$chrom,
                       strand = rep2$strand,
                      ranges = IRanges(start = rep2$start, width = 1),
                      mcols = DataFrame(rep2[, c("readCount", "percentMeth")]))

save(methyl_rep1, file = file.path(save_dir, "methyl_rep1.RData"))
save(methyl_rep2, file = file.path(save_dir, "methyl_rep2.RData"))
rm(methyl_rep1, methyl_rep2)

## ----read_histones, eval=FALSE-------------------------------------------
#  histone_patterns <- list(H3k09me3 = "*H3k09me3",
#                          H3k27ac = "*H3k27ac",
#                          H3k27me3 = "*H3k27me3",
#                          H3k36me3 = "*H3k36me3",
#                          H3k4me3_r1 = "*H3k4me3",
#                          H3k4me3_r2 = "*H3k04me3")
#  
#  histone_dir <- file.path(raw_data_dir, "histone_marks")
#  histone_names <- c("chrom", "start", "end", "name", "score", "strand", "signal", "pvalue", "qvalue", "other")
#  histone_marks <- lapply(histone_patterns, function(in_pattern){
#    use_file <- dir(histone_dir, pattern = in_pattern)
#    tmp <- read.table(file.path(histone_dir, use_file), header = FALSE, sep = "\t")
#    names(tmp) <- histone_names
#    GRanges(seqnames = tmp$chrom,
#    			  ranges = IRanges(start = tmp$start, end = tmp$end),
#  					mcols = DataFrame(tmp[, c("score", "signal", "pvalue")]))
#  })
#  histone_marks <- GRangesList(histone_marks)
#  save(histone_marks, file = file.path(save_dir, "histone_marks.RData"))
#  rm(histone_marks)

## ----get_expression------------------------------------------------------
expr_data <- dataTable(getGEO(GEO = "GSM307014"))@table
expr_data$VALUE <- as.numeric(expr_data$VALUE)
expr_data[, "DETECTION P-VALUE"] <- as.numeric(expr_data[, "DETECTION P-VALUE"])
rownames(expr_data) <- as.character(expr_data$ID_REF)
expr_data$ID_REF <- rownames(expr_data)
save(expr_data, file = file.path(save_dir, "expr_data.RData"))

## ----get_tss_windows, eval = FALSE---------------------------------------
#  tss_file <- file.path(raw_data_dir, "ensGene_TTSS.csv")
#  tss_data <- read.table(tss_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
#  tss_regions <- GRanges(seqnames = tss_data$chrom,
#                         strand = tss_data$strand,
#                         ranges = IRanges(start = tss_data$txStart,
#                                          end = tss_data$txEnd),
#                         names = tss_data$name)
#  rm(tss_data)
#  
#  tss_windows <- GRanges(seqnames = seqnames(tss_regions),
#                         strand = strand(tss_regions),
#                         ranges = IRanges(start = start(tss_regions) - 1000,
#                                          end = start(tss_regions) + 1000))
#  names(tss_windows) <- mcols(tss_regions)$names
#  save(tss_windows, file = file.path(save_dir, "tss_windows.RData"))
#  rm(tss_windows, tss_regions)

## ------------------------------------------------------------------------
Sys.time()
sessionInfo()

