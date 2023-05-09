# > END analysis pfu < #

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# functions ----
## Gff file ====
get_gff <- function(my_gff, my_Strand = c("+", "-")){
  
  my_gff %>%
    dplyr::mutate(seqnames = "CP023154.1",
                  type = "CDS") %>%
    dplyr::rename(gene_id = locus_tag) %>%
    dplyr::select(seqnames, gene_id, start, end, strand, type) %>%
    dplyr::filter(strand == my_Strand) %>%
    rowwise() %>%
    dplyr::mutate(start = ifelse(my_Strand == "+",start, start - 300),
                  end = ifelse(my_Strand == "+",end + 300, end)) %>%
    ungroup() %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  
}

# > min normalised coverage: 5
# > min detected exact position in 3 of 4 replicates
read_term <- function(input_F, min_coverage = 5){
  
  vroom(input_F,  col_names = F, show_col_types = F, num_threads = 8) %>%
    dplyr::select(X1, X2, X3, X5, X6, X10, X11, X13, X14, X15, X16) %>%
    dplyr::rename(chrom = X1, peak_start = X2, peak_end = X3, prominence = X5, strand = X6, 
                  width = X10, replicate = X11, start_cov  = X13, end_cov = X14, cov_peak = X15, cov_width = X16) %>%
    group_by(chrom, peak_start, peak_end, replicate) %>%
    dplyr::filter(cov_peak == max(cov_peak)) %>%
    dplyr::filter(cov_peak >= min_coverage) %>%
    group_by(chrom, start_cov, end_cov) %>%
    dplyr::filter(n() >= 3) %>%
    group_by(chrom, start_cov, end_cov, replicate) %>%
    dplyr::filter(n_distinct(end_cov) == 1) %>%
    group_by(chrom, start_cov, end_cov) %>%
    dplyr::mutate(decision_v = if_else(strand == "+", 
                                       max(end_cov), min(end_cov))) %>% 
    dplyr::filter(end_cov == decision_v) %>%                             
    ungroup() %>%
    arrange(end_cov) %>%
    dplyr::rename(seqnames = chrom, start = start_cov, end = end_cov) %>%
    dplyr::group_by(seqnames, start, end, strand) %>%
    summarise(cov_peak = mean(cov_peak)) %>%
    dplyr::select(seqnames, start, end, strand, cov_peak) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
}


## TTS annotation ====
get_peaks <- function(in_File_F, in_File_R, my_Set, min_coverage){
  
  input_peaks_F <- read_term(in_File_F, min_coverage) 
  input_peaks_R <- read_term(in_File_R, min_coverage) 
  
  annot_forward <- annotatr::annotate_regions(input_peaks_F, pfu_gff_forward) %>%
    as_tibble() %>%
    dplyr::mutate(end_feature = annot.end - 300) %>%
    dplyr::mutate(UTR3 = end-end_feature) %>%
    dplyr::rename(start_feature = annot.start, TTS = end, gene_id = annot.gene_id, biotype = annot.type) %>%
    dplyr::select(seqnames, start_feature, end_feature, TTS, cov_peak, UTR3, strand, gene_id, biotype) %>%
    group_by(gene_id) 
  
  annot_forward_h <- annot_forward %>%
    filter(UTR3 >= 0) %>%
    dplyr::filter(TTS == TTS[which.max(cov_peak)]) %>%
    dplyr::mutate(is_utr3 = T)
  
  annot_forward_f <- left_join(annot_forward, annot_forward_h) %>%
    dplyr::mutate(TTS_class = if_else(TTS >= end_feature & cov_peak == max(cov_peak), "primary", "rest")) %>%
    dplyr::mutate(TTS_class = if_else(TTS_class == "primary", "primary", 
                                      if_else(sum(TTS_class == "primary") == 0 & !is.na(is_utr3), "primary",
                                              if_else(TTS >= end_feature & is.na(is_utr3), "secondary",
                                                      if_else(TTS > start_feature & TTS < end_feature, "intragenic", "other")))))  
  
  
  annot_reverse <- annotatr::annotate_regions(input_peaks_R, pfu_gff_reverse) %>%
    as_tibble() %>%
    dplyr::mutate(start_feature = annot.start + 300) %>%
    dplyr::mutate(UTR3 = start_feature-end) %>%
    dplyr::rename(end_feature = annot.end, TTS = end, gene_id = annot.gene_id, biotype = annot.type) %>%
    dplyr::select(seqnames, start_feature, end_feature, TTS, cov_peak, UTR3, strand, gene_id, biotype) %>%
    group_by(gene_id)
  
  annot_reverse_h <- annot_reverse %>%
    filter(UTR3 >= 0) %>%
    dplyr::filter(TTS == TTS[which.max(cov_peak)]) %>%
    dplyr::mutate(is_utr3 = T)
  
  annot_reverse_f <- left_join(annot_reverse, annot_reverse_h) %>%
    dplyr::mutate(TTS_class = if_else(TTS <= start_feature & cov_peak == max(cov_peak), "primary", "rest")) %>%
    dplyr::mutate(TTS_class = if_else(TTS_class == "primary", "primary", 
                                      if_else(sum(TTS_class == "primary") == 0 & !is.na(is_utr3), "primary",
                                              if_else(TTS <= start_feature & is.na(is_utr3), "secondary",
                                                      if_else(TTS > start_feature & TTS < end_feature, "intragenic", "other")))))  
  
  return(bind_rows(annot_forward_f, annot_reverse_f) %>%
           dplyr::select(-is_utr3) %>%
           dplyr::mutate(set = my_Set))
}


# data ----
ext_dir <- "/path/to/working/directory/"

## genome data ====
### gff file ####
pfu_gff <- read.gff(here("data/genome/GCF_008245085.1_ASM824508v1_genomic.gff")) %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(locus_tag = str_split_fixed(str_split_fixed(attributes, ";old_",2)[,1],"locus_tag=",2)[,2],
                old_locus_tag = str_split_fixed(attributes, "old_locus_tag=",2)[,2],
                width = abs(start-end)) %>%
  dplyr::select(locus_tag,old_locus_tag,start, end, strand,width)

### fasta file ####
pfu_fasta <- readDNAStringSet(here("data/genome/GCF_008245085.1_ASM824508v1_genomic.fna"))
names(pfu_fasta) <- "chr"

### Modify gff file for forward and reverse strand ####
pfu_gff_forward <- get_gff(pfu_gff, "+")
pfu_gff_reverse <- get_gff(pfu_gff, "-")

## Peak calling results (3'end annotation) - output from termseq_peaks ====
files_peak <- list.files(paste0(ext_dir, "termseq_peaks23"), pattern = ".narrowPeak.counts", full.names = T, recursive = T)

tts_peaks_all <- pmap_dfr(list(files_peak[str_detect(files_peak, "plus")],
                               files_peak[str_detect(files_peak, "neg")], 
                               "T95",
                               6), 
                          get_peaks)

tts_peaks_all_pr <- tts_peaks_all %>%
  dplyr::filter(TTS_class == "primary") %>%
  distinct(gene_id, .keep_all = T) %>%
  dplyr::mutate(UTR3_sequence = if_else(strand == "+",str_replace_all(as.character(pfu_fasta$chr[(TTS-20):(TTS+20)]), "T", "U"),
                                        str_replace_all(as.character(reverseComplement(pfu_fasta$chr[(TTS-20):(TTS+20)])), "T", "U"))) 

vroom_write(tts_peaks_all_pr, here("data/Rdata/end3_pfu_short_read_95.tsv"))

# Supplementary Table 2 ----
## annotation ====
pfu_annotation_small <- vroom(here("data/Rdata/pfu_annotation_small.tsv"))

## end3 data ====
pfu_utr3_fine <- vroom(here("data/Rdata/end3_pfu_short_read_95.tsv")) %>%
  dplyr::mutate(locus_tag = gene_id) %>%
  left_join(pfu_annotation_small, by = c("locus_tag")) %>%
  dplyr::rename(pfu_tag = pf_name,
                annotation = description) %>%
  dplyr::select(locus_tag, old_locus_tag, pfu_tag, annotation, start_feature:strand)


wb <- createWorkbook()
addWorksheet(wb, "TERM_seq")
writeDataTable(wb, "TERM_seq", pfu_utr3_fine)
saveWorkbook(wb, here("supplemental_tables/Supplementary_Table_2.xlsx"))

