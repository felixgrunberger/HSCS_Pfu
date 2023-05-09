# > operon plotting < #

# libraries ----
library(here)
source(here("scripts/load_libraries.R"))

# functions ----
read_single_cov_short <- function(input){
  frw_cov <- vroom(cov_files[str_detect(cov_files, input) & str_detect(cov_files, "reverse")],
                   show_col_types = F,num_threads = 4, col_names = F) %>%
    dplyr::rename(chr = 1, pos = 2, cov = 3) %>%
    dplyr::mutate(strand = "+",
                  cpm = cov/sum(cov) * 1000000)
  rev_cov <- vroom(cov_files[str_detect(cov_files, input) & str_detect(cov_files, "forward")],
                   show_col_types = F,num_threads = 4, col_names = F) %>%
    dplyr::rename(chr = 1, pos = 2, cov = 3) %>%
    dplyr::mutate(strand = "-",
                  cpm = cov/sum(cov) * 1000000)
  return(bind_rows(frw_cov, rev_cov) %>% dplyr::mutate(sample_ids = input))
}

read_single_cov_long <- function(input){
  
  frw_cov <- vroom(cov_files_long[str_detect(cov_files_long, input) & str_detect(cov_files_long, "forward")],
                   show_col_types = F,num_threads = 4, col_names = F) %>%
    dplyr::rename(chr = 1, pos = 2, cov = 3) %>%
    dplyr::mutate(strand = "+",
                  cpm = cov/sum(cov) * 1000000)
  rev_cov <- vroom(cov_files_long[str_detect(cov_files_long, input) & str_detect(cov_files_long, "reverse")],
                   show_col_types = F,num_threads = 4, col_names = F) %>%
    dplyr::rename(chr = 1, pos = 2, cov = 3) %>%
    dplyr::mutate(strand = "-",
                  cpm = cov/sum(cov) * 1000000)
  return(bind_rows(frw_cov, rev_cov) %>% dplyr::mutate(sample_ids = input))
}

read_bam_files <- function(inputBAM, method){
  
  # read in files
  init <- readGAlignments(inputBAM, use.names = T, param = ScanBamParam(tag=c("NM"), what="mapq"))
  init_t <- GenomicAlignments::as.data.frame(init) %>%
    dplyr::mutate(minion_read_name = names(init),
                  mapped_gene = seqnames) 
  
  left  <- paste(str_split_fixed(string = init_t$cigar, pattern = "M", n = 2)[,1],"M", sep = "")
  right <- paste(str_split_fixed(string = init_t$cigar, pattern = "M", n = 2)[,2],"1M", sep = "")
  
  #................................calculate cigar tables / SOFT AND HARD CLIPPING!!!
  init_t$soft_l <- as_tibble(cigarOpTable(left))$S
  init_t$hard_l <- as_tibble(cigarOpTable(left))$H
  init_t$soft_r <- as_tibble(cigarOpTable(right))$S
  init_t$hard_r <- as_tibble(cigarOpTable(right))$H
  
  # calculate number of aligned reads based on CIGAR operations (M,I)
  init_t$aligned_reads <- unlist(lapply(explodeCigarOpLengths(init_t$cigar, ops = c("M", "I")), function(x) sum(x)))
  
  # calc read identity..
  init_t_final <- init_t %>%
    dplyr::mutate(identity = (1 - NM/aligned_reads)*100) %>%
    dplyr::group_by(minion_read_name) %>%
    dplyr::filter(identity == max(identity),
                  aligned_reads == max(aligned_reads)) %>%
    dplyr::distinct(minion_read_name, .keep_all = T) %>%
    dplyr::mutate(sample = method,
                  gene = str_split_fixed(mapped_gene,"-",2)[,2])
  
  # return table
  return(init_t_final)
}

annotate_bams <- function(input_mapped, input_remapped, dataset){
  mapped_t   <- read_bam_files(input_mapped, dataset)
  remapped_t <- read_bam_files(input_remapped, dataset)
  mapped_t_a <- mapped_t %>%
    dplyr::select(-mapped_gene, -gene) %>%
    left_join(remapped_t %>%
                dplyr::select(minion_read_name, mapped_gene, gene), 
              by = "minion_read_name") 
}


plot_full_operon <- function(selected_gene, selected_range){
  
  w_strand <- pfu_gff %>%
    dplyr::filter(locus_tag %in% selected_gene) %>%
    dplyr::select(strand) %>%
    deframe()
  
  ### ont single reads ####
  ont_w <- ont_data %>%
    dplyr::filter(sample %in% "RNA_95_2") %>%
    dplyr::filter(start %in% selected_range,
                  end %in% selected_range,
                  strand %in% w_strand) %>%
    arrange(start,aligned_reads) %>%
    dplyr::mutate(rown = 1:nrow(.))
  
  ### coverage  ####
  cov_w <- bind_rows(cov_illumina %>%
                       dplyr::filter(pos %in% selected_range,
                                     strand %in% w_strand) %>%
                       left_join(sample_ids_clean, by = "sample_ids") %>%
                       group_by(Condition, chr, pos, strand) %>%
                       summarise(mean_cpm = mean(cpm, na.rm = T)) %>%
                       dplyr::mutate(method = "short"),
                     cov_long %>%
                       dplyr::filter(pos %in% selected_range,
                                     strand %in% w_strand) %>%
                       dplyr::mutate(Condition = "Ctrl") %>%
                       group_by(Condition, chr, pos, strand) %>%
                       summarise(mean_cpm = mean(cpm, na.rm = T)) %>%
                       dplyr::mutate(method = "long"))
  
  ### selected sites ####
  tss_w <- tss %>%
    dplyr::filter(end5 %in% selected_range,
                  strand %in% w_strand)
  tts_w <- tts %>%
    dplyr::filter(end3 %in% selected_range,
                  strand %in% w_strand)
  
  ### gff ####
  pfu_gff_w <- pfu_gff %>%
    dplyr::filter(start %in% selected_range,
                  end %in% selected_range) 
  
  ### plot ###
  ggplot(data = cov_w,
         aes(x = pos, y = log10(mean_cpm+1), 
             fill = Condition, 
             color = Condition)) +
    facet_grid(rows = vars(method), scales = "free_y") +
    geom_area(position = position_identity(), alpha = 0.5) +
    geom_vline(xintercept = tss_w$end5, linetype = "solid") +
    geom_vline(xintercept = tts_w$end3, linetype = "dashed") +
    scale_fill_manual(values = c("#B3B3B3","#E9995E")) +
    scale_color_manual(values = c("black","#E9995E")) +
    theme_light() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"),
          panel.grid.minor = element_blank()) +
    ggplot(data = pfu_gff_w) +
    geom_segment(aes(x = start, xend = end, y = 1, yend = 1)) +
    ggplot() +
    geom_segment(data = ont_w, aes(y = rown, yend = rown, x = start, xend = end), size = .5, alpha = 0.5) +
    theme_light() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = "dashed"),
          panel.grid.minor = element_blank()) +
    patchwork::plot_layout(ncol = 1)
  
}

# data ----
## Sample info ====
sample_ids    <- read_xlsx(here("data/samples/sample_id_info.xlsx")) %>%
  dplyr::mutate(Condition = case_when(Condition == "CS_short" ~ "CS1",
                                      Condition == "CS_middle" ~ "CS2",
                                      Condition == "CS_long" ~ "CS3",
                                      Condition == "CS_Rec_short" ~ "CSR",
                                      Condition == "Ctrl" ~ "Ctrl",
                                      Condition == "HS1" ~ "HS1",
                                      Condition == "HS2" ~ "HS2",
                                      Condition == "HS3" ~ "HSR")) 

## Bad samples info ====
remove_rna_samples <- data.table(Condition = c("HS1", "HS2", "HSR", "CSR"),
                                 bio_rep = c(1, 4, 1, 2),
                                 remove = T) %>%
  left_join(sample_ids) %>%
  dplyr::select(sample_ids) %>%
  deframe()

sample_ids_clean <- sample_ids %>%
  dplyr::filter(!sample_ids %in% remove_rna_samples)

## Illumina coverage files (calculated using samtools depth) ====
cov_files <- list.files(paste0(ex_dir, "data/coverage"), recursive = T, full.names = T, pattern = ".tsv")
cov_samples <- str_split_fixed(str_split_fixed(cov_files, "\\/", 8)[,7], "_", 6)[,3]
ctrl <- unique(cov_samples[cov_samples %in% sample_ids_clean$sample_ids[sample_ids_clean$Condition == "Ctrl"]])
hs1 <- unique(cov_samples[cov_samples %in% sample_ids_clean$sample_ids[sample_ids_clean$Condition == "HS1"]])
cov_illumina <- pmap_dfr(list(c(ctrl, hs1)), read_single_cov_short) %>% 
  dplyr::mutate(method = "short")

## Nanopore coverage files (calculated as described in the documentation) ====
cov_files_long <- list.files(paste0(ex_dir_long, "coverage/ont"), recursive = T, full.names = T, pattern = ".tsv")
cov_samples_long <- str_split_fixed(str_split_fixed(cov_files_long, "\\/", 8)[,7], "_", 2)[,2]
normal_cond <- c("95_1", "95_2")
cov_long <- pmap_dfr(list(normal_cond), read_single_cov_long) %>%
  dplyr::mutate(method = "long")

## Nanopore single read files ====
### mapped files | simply from unfiltered fastq to minimap2/remap ####
files          <- list.files(paste0(dir,"mapped/new2023"), recursive = T, full.names = T,pattern = ".sorted.bam$")
files95 <- files[str_detect(files, "_95")]
mapped_frame   <- pmap_dfr(list(files95[!str_detect(files95, "remapped")],
                                files95[str_detect(files95, "remapped")],
                                unique(str_split_fixed(str_split_fixed(files95, "\\/", n = 8)[,7],"_fu",2)[,1])),
                           annotate_bams)

mapped_frame_annotation <- mapped_frame %>%
  left_join(pfu_annotation %>%
              distinct(locus_tag, .keep_all = T), 
            by = c("gene" = "locus_tag"))

### save data frame ####
fwrite(mapped_frame_annotation, here("data/Rdata/mapped_data_trimmed.tsv"), sep = "\t",col.names = T, nThread = 8)
mapped_frame <- vroom( here("data/Rdata/mapped_data_trimmed.tsv"))

ont_data <- vroom(here("data/Rdata/mapped_data_trimmed.tsv")) %>%
  ungroup() %>%
  dplyr::filter(njunc == 0) %>%
  dplyr::select(strand.x, start, end, sample, aligned_reads) %>%
  dplyr::rename(strand = strand.x)

## TSS info (from Grünberger et al. 2019 - Frontiers in Microbiology) ====
tss <- vroom(here("data/genome/TSS_MasterTable.tsv"), show_col_types = F) %>%
  dplyr::filter(Primary == 1) %>%
  dplyr::filter(Locus_tag %in% pfu_gff$old_locus_tag) %>%
  dplyr::rename(old_locus_tag = Locus_tag,
                end5 = Pos,
                strand = SuperStrand) %>%
  left_join(pfu_gff %>% dplyr::select(locus_tag, old_locus_tag)) %>%
  dplyr::filter(locus_tag != "") %>%
  dplyr::select(locus_tag, end5, strand)

## TTS info (from 3' end detection) ====
tts <- vroom(here("data/Rdata/end3_pfu_short_read_95.tsv")) %>%
  dplyr::mutate(locus_tag = gene_id) %>%
  dplyr::rename(end3 = TTS) %>%
  dplyr::filter(locus_tag != "") %>%
  dplyr::select(locus_tag, end3, strand)

## Gene selection ====

### hsp20 <- "PFDSM3638_RS09505"
hsp20 <- "PFDSM3638_RS09505"
plot_range_hsp20 <- 1698815:1703444
plot_full_operon(hsp20,plot_range_hsp20)

### thermo <- "PFDSM3638_RS09970"
thermo <- "PFDSM3638_RS09970"
plot_range_thermo <- 1790478:1797743
plot_full_operon(thermo,plot_range_thermo)