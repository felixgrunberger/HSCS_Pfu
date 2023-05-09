# > downstream RNA-seq analysis - DESeq2 < # 

# load libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# read in data ----
## make custom gtf ====

## detailed genome annotation ====
pfu_fasta <- readDNAStringSet(here("data/genome/GCF_008245085.1_ASM824508v1_genomic.fna"))

pfu_gff <- read.gff(here("data/genome/GCF_008245085.1_ASM824508v1_genomic.gff")) %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(locus_tag = str_split_fixed(str_split_fixed(attributes, ";old_",2)[,1],"locus_tag=",2)[,2],
                old_locus_tag = str_split_fixed(attributes, "old_locus_tag=",2)[,2],
                width = abs(start-end),
                biotype = str_split_fixed(str_split_fixed(attributes, "gene_biotype=",2)[,2],";locus_tag=",2)[,1]) %>%
  dplyr::filter(biotype == "protein_coding") %>%
  dplyr::select(locus_tag,old_locus_tag,start, end, strand,width, biotype)

## make custom gtf ====
custom_gtf <- data.table("1" = "NZ_CP023154.1", 
                         "2" = "NCBI",
                         "3" = "all",
                         "4" = pfu_gff$start,
                         "5" = pfu_gff$end,
                         "6" = NA,
                         "7" = pfu_gff$strand,
                         "8" = NA,
                         "9" = paste0("transcript_id ",pfu_gff$locus_tag, "; ",
                                      "gene_id ",pfu_gff$locus_tag, "; ",
                                      "gene_name ",pfu_gff$locus_tag, ";"))

write.table(custom_gtf, file = here("data/genome/pfu_custom_CDS.gtf"),
            col.names = F, row.names = F, quote = F, sep = "\t")

gtf_file <- here("data/genome/pfu_custom_CDS.gtf")

# DESeq2 pipeline ----

## mapped files ====
bam_list    <- list.files(path = "/path/to/workspace/mapped/", 
                          recursive = T, full.names = T, pattern = ".sorted.bam$")

## calculate count matrix using featurecounts and label features ====
### calc ####
counts <- featureCounts(bam_list,verbose = F, annot.ext = gtf_file, strandSpecific = 2,
                        primaryOnly = T,isGTFAnnotationFile = T, nthreads = 8,
                        GTF.featureType = "all", allowMultiOverlap = F, 
                        isLongRead = F)$counts

colnames(counts) <- str_split_fixed(str_split_fixed(colnames(counts), "_S",2)[,1], "_",3)[,3]

### save ####
save(counts, file = here("data/Rdata/featurecounts_rnaseq.Rdata"))

### load pre-calculated featurecounts object ====
load(here("data/Rdata/featurecounts_rnaseq.Rdata"))

## Sample design  ====
### Read in info table ####
sample_ids    <- read_xlsx(here("data/samples/sample_id_info.xlsx")) %>%
  dplyr::mutate(Condition = case_when(Condition == "CS_short" ~ "CS1",
                                      Condition == "CS_middle" ~ "CS2",
                                      Condition == "CS_long" ~ "CS3",
                                      Condition == "CS_Rec_short" ~ "CSR",
                                      Condition == "Ctrl" ~ "Ctrl",
                                      Condition == "HS1" ~ "HS1",
                                      Condition == "HS2" ~ "HS2",
                                      Condition == "HS3" ~ "HSR")) 
sample_design <- colnames(counts)

### make design matrix ####
colDataRaw <- data.table(sample_ids = sample_design) %>%
  left_join(sample_ids) %>%
  dplyr::filter(!is.na(time))

## DESeq object ====
### calc DESeq2 ####
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts[,colnames(counts) %in% colDataRaw$sample_ids], 
                                      colData = colDataRaw,
                                      design = ~Condition)

dds <- dds[rowSums(counts(dds)) > 1,]
dds_adjusted <- estimateSizeFactors(dds)

### get normalized counts table ####
normalized_counts <- counts(dds_adjusted, normalized = T)
vroom_write(normalized_counts %>%
              as.data.frame() %>%
              rownames_to_column("gene"), 
            here("data/Rdata/deseq_normalized_counts.tsv"))
