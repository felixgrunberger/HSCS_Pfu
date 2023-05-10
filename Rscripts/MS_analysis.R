# libraries ----
library(here)
source(here("Rscripts/load_libraries.R"))

# function ----
get_pairwise_DE_table_MS <- function(input_test, input_counts = dat.log, input_design = design, input_ctrl = "Ctrl"){
  
  contrast_x <- c(paste0(input_test, "-", input_ctrl))
  contrast   <-   makeContrasts(contrasts=contrast_x, levels=input_design)
  fit1 <- lmFit(input_counts, input_design)
  fit2 <- contrasts.fit(fit1, contrasts = contrast)
  fit3 <- eBayes(fit2)
  
  output <- topTable(fit3, sort.by = "P", n = Inf) %>%
    rownames_to_column("rowid") %>%
    as_tibble() %>%
    left_join(protein_names, by = "rowid") %>%
    arrange(desc(logFC)) %>%
    mutate(typeReg = case_when(logFC <= 0 & adj.P.Val < 0.05 ~ "down",
                               logFC >=  0 & adj.P.Val < 0.05 ~ "up",
                               adj.P.Val >= 0.05 ~ "rest"),
           comp = paste0("Condition_Ctrl_vs_", input_test)) %>%
    dplyr::select(c(logFC, adj.P.Val, locus_tag, typeReg, comp)) %>%
    dplyr::mutate(set2 = input_test) %>%
    dplyr::rename(log2FC = 1, padj = 2)
  return(output)
}

# data ----
## load MS quantities saved in Table S3 ====
xl_file <- read_xlsx(here("supplemental_tables/Supplementary_Table_3.xlsx"), sheet = "QUANTITY")

## log2 transform count data ====
### counts ####
counts_MS <- xl_file[,-c(1:4)]

### extract Condition & Replicate from col name ####
ms_sample_info <- data.table(run_label = colnames(counts_MS)) %>%
  dplyr::mutate(Condition = str_remove_all(str_split_fixed(run_label, "_biorep", 2)[,1],"_"),
                bio_rep = str_split_fixed(str_split_fixed(run_label, "_biorep", 2)[,2], "_",2)[,1])

### log2 ####
ex <- log2(counts_MS)

### calc PCA ####
pca <- prcomp(t(ex), center = TRUE, scale. = FALSE)

pca_tibble <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column("run_label") %>%
  left_join(ms_sample_info, by = "run_label") %>%
  dplyr::filter(!is.na(Condition))

eigs <- pca$sdev^2
pc1_var <- round(eigs[1] / sum(eigs) *100,1)
pc2_var <- round(eigs[2] / sum(eigs) *100,1)

# plot ----
ggplot(pca_tibble,  
       aes(x = PC1, y = PC2, fill = Condition, group = Condition, 
           shape = as.factor(bio_rep))) +
  geom_mark_ellipse(aes(color = Condition, group = Condition), 
                    alpha = 0.25) +
  geom_point(size = 4, color = "black") +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Replicate") +
  guides(fill = guide_legend(override.aes=list(shape=c(21))),
         shape = guide_legend(override.aes=list(shape=c(21,22,23,24),
                                                fill = "grey50"))) +
  theme_linedraw() +
  scale_color_manual(values = batlow_custom) +
  scale_fill_manual(values = batlow_custom) +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(-20,20)) +
  scale_y_continuous(limits = c(-20,20)) +
  xlab(paste0("PC1: ", pc1_var, "% variance")) +
  ylab(paste0("PC2: ", pc2_var, "% variance")) 

# Remove outliers ----
rem_ms_pca <- data.table(bio_rep = c(1,4,2,3,1),
                         Condition = c("CS2", "CS3", "CSR", "HS2", "HSR"),
                         label = "remove") %>%
  left_join(ms_sample_info %>% 
              dplyr::mutate(bio_rep = as.numeric(bio_rep)), 
            by = c("bio_rep", "Condition")) %>%
  pull(run_label)

counts_MS_remBad <- counts_MS[!colnames(counts_MS) %in% rem_ms_pca]

ex2 <- log2(counts_MS_remBad)
pca2 <- prcomp(t(ex2), center = TRUE, scale. = FALSE)
eigs <- pca2$sdev^2
pc1_var <- round(eigs[1] / sum(eigs) *100,1)
pc2_var <- round(eigs[2] / sum(eigs) *100,1)

pca_tibble2 <- pca2$x %>%
  as.data.frame() %>%
  rownames_to_column("run_label") %>%
  left_join(ms_sample_info, by = "run_label") 

ggplot(pca_tibble2,  
       aes(x = PC1, y = PC2, fill = Condition, group = Condition, 
           shape = as.factor(bio_rep))) +
  geom_mark_ellipse(aes(color = Condition, group = Condition), 
                    alpha = 0.25) +
  geom_point(size = 4, color = "black") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  guides(fill = guide_legend(override.aes=list(shape=c(21))),
         shape = guide_legend(override.aes=list(shape=c(21,22,23,24),
                                                fill = "grey50"))) +
  theme_linedraw() +
  scale_color_manual(values = batlow_custom) +
  scale_fill_manual(values = batlow_custom) +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(expand = c(.2,0)) +
  scale_y_continuous(expand = c(.2,0)) +
  xlab(paste0("PC1: ", pc1_var, "% variance")) +
  ylab(paste0("PC2: ", pc2_var, "% variance")) 


## export normalized counts tables ====
### get protein names ####
protein_names <- xl_file[,1] %>%
  dplyr::mutate(rowid = as.character(1:nrow(.)))

normalized_counts_MS_remBad <- counts_MS_remBad %>%
  rownames_to_column("rowid") %>%
  as_tibble() %>%
  left_join(protein_names %>% dplyr::select(locus_tag, rowid), by = "rowid") %>%
  dplyr::select(-rowid) %>%
  pivot_longer(-locus_tag) %>%
  left_join(ms_sample_info, by = c("name" = "run_label")) %>%
  group_by(Condition, locus_tag) %>%
  summarise(mean_counts = mean(value, na.rm = T)) %>%
  dplyr::rename(gene = locus_tag) %>%
  dplyr::select(Condition, gene, mean_counts) %>%
  ungroup()

save(normalized_counts_MS_remBad,
     file = here("data/Rdata/normalized_counts_MS_remBad.Rdata"))

## export normalized counts tables with IBAQ Scores ====
ibaq_MS <- read_xlsx(here("supplemental_tables/Supplementary_Table_3.xlsx"), sheet = "IBAQ")
ibaq_MS_remBad <- ibaq_MS[!colnames(ibaq_MS) %in% rem_ms_pca]


normalized_ibaq_MS_remBad <- ibaq_MS_remBad %>%
  dplyr::select(-c(old_locus_tag, pfu_tag, annotation)) %>%
  pivot_longer(-locus_tag) %>%
  dplyr::mutate(value = as.numeric(as.character(value))) %>%
  left_join(ms_sample_info, by = c("name" = "run_label")) %>%
  dplyr::filter(value != 0) %>%
  group_by(Condition, locus_tag) %>%
  summarise(mean_counts = mean(value, na.rm = T))  %>%
  dplyr::rename(gene = locus_tag) %>%
  dplyr::select(Condition, gene, mean_counts) %>%
  ungroup()

save(normalized_ibaq_MS_remBad,
     file = here("data/Rdata/normalized_ibaq_MS_remBad.Rdata"))

# prepare counts table and statistical analysis using limma ----
## prepare MS data ====
dat.log = log2(counts_MS_remBad)

## design table ====
cond <- data.table(run_label = colnames(dat.log)) %>%
  left_join(ms_sample_info, by = "run_label") %>%
  dplyr::select(Condition) %>%
  deframe() %>%
  as.factor()

design = model.matrix(~0+cond) 
colnames(design) = gsub("cond","",colnames(design))

## make contrast table ====
all_conds <- levels(as.factor(cond))
test_conds <- all_conds[all_conds != "Ctrl"]
protein_names <- xl_file[,1] %>%
  dplyr::mutate(rowid = as.character(1:nrow(.)))

## perform limma ====
limma_all <- pmap_dfr(list(test_conds),get_pairwise_DE_table_MS)

## export FC tables ====
comp_table_MS <- limma_all %>%
  dplyr::rename(gene = locus_tag, Condition = set2) %>%
  dplyr::select(gene, log2FC, padj, typeReg, Condition) %>%
  arrange(gene) 

vroom::vroom_write(comp_table_MS, 
                   file = here("data/Rdata/comp_table_MS.tsv"))
