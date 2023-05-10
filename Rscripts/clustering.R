# libraries ----
library(here)
library(NbClust)
library(factoextra)
library(cluster)
source(here("Rscripts/load_libraries.R"))


# functions ----
plot_heatmap_arcog_simple <- function(input, set, sig_value){
  
  my_colors <- colorRampPalette(c("#F9F9F9", "#CEC9EC", "#ABA0DE", "#8B79D3", "#6B5DA5", "#4D4573", "#312D44"))
  
  
  ggplot(data = input %>%
           dplyr::mutate(plot_cat = paste0(arcog_cat_big, "_", plot_order,"_",category),
                         circle_color = case_when(-log10(over_represented_pvalue) > (max(-log10(over_represented_pvalue))/2) ~ "white",
                                                  .default = "black")),
         aes(y = (plot_cat), 
             x = cluster, 
             fill = -log10(over_represented_pvalue), 
             size = numDEInCat, shape = over_represented_pvalue < sig_value)) +
    geom_tile(color = "black", size = 0.2) +
    geom_point(aes(color = circle_color), stroke = 0.5) +
    scale_shape_manual(values = c(NA,21)) +
    theme_light() +
    theme_linedraw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    scale_color_manual(values = c("black", "white")) +
    ggtitle(set)  +
    coord_equal() +
    scale_fill_gradientn(colours = my_colors(100)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.ticks = element_blank()) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    guides(shape = "none", color = "none") +
    xlab("") +
    ylab("") 
}

fc_to_wide <- function(inP, type, cond){
  inP %>%
    dplyr::filter(Condition %in% cond) %>%
    dplyr::mutate(Condition = paste0(type, "_", Condition)) %>%
    distinct(gene, Condition, .keep_all = T) %>%
    dplyr::select(gene, log2FC, Condition) %>%
    pivot_wider(names_from = Condition, values_from = log2FC)
}


goseq_cluster <- function(inputF, myCluster, arcog_table = pfu_arcog, annotation_table = pfu_annotation){
  
  complete_set <- inputF %>%
    dplyr::distinct(gene,.keep_all = T) %>%
    dplyr::filter(!is.na(gene))
  
  # > all genes 
  assayed.genes <- complete_set %>%
    pull(gene) 
  
  # > get differentially expressed genes
  de.genes      <- complete_set %>%
    dplyr::filter(cluster %in% myCluster) %>%
    pull(gene) 
  
  # > get de-genes and background set in a list
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  
  # > add length information for calculations
  lengthGenes <- complete_set %>%
    left_join(annotation_table %>%
                distinct(locus_tag, .keep_all = T), 
              by = c("gene" = "locus_tag")) %>%
    distinct(gene, .keep_all = T) %>%
    dplyr::mutate(width = Length) %>%
    pull(width)
  
  # > calc arcog enrichment
  pwf <- goseq::nullp(gene.vector, bias.data = lengthGenes,'ensGene',plot.fit=FALSE)
  
  # > add arcog identifier
  category_mapping <- arcog_table %>%
    left_join(pfu_annotation %>%
                dplyr::select(locus_tag, old_locus_tag) %>%
                distinct(locus_tag, .keep_all = T), 
              by = c("new" = "old_locus_tag")) %>%
    dplyr::rename(arCOG = COG.x) %>%
    dplyr::filter(locus_tag %in% names(gene.vector)) %>%
    dplyr::select(locus_tag, arCOG) %>%
    as.data.frame()
  
  category.vector <- category_mapping$arCOG
  names(category.vector) <- as.factor(category_mapping$locus_tag)
  
  goseq_results <- goseq::goseq(pwf, gene2cat = category_mapping, use_genes_without_cat=TRUE) %>%
    as_tibble() %>%
    mutate(cluster = myCluster) %>%
    left_join(arcog_info, by = c("category" = "arcog_cat")) %>%
    dplyr::mutate(expected = numInCat*sum(numDEInCat)/sum(numInCat),
                  deviation_from_expected = numDEInCat/expected)
  
  return(goseq_results)
}

calc_z_pca_perCol <- function(inP, type, cond){
  
  inP %>%
    dplyr::filter(Condition %in% cond) %>%
    distinct(gene, Condition, .keep_all = T) %>%
    group_by(gene) %>%
    ungroup() %>%
    dplyr::select(gene, log2FC, Condition) %>%
    dplyr::mutate(Condition = paste0(type, "_", Condition)) %>%
    group_by(gene) %>%
    mutate(log2FC_scaled = (log2FC - mean(log2FC, na.rm = T))/sd(log2FC, na.rm = T)) %>% 
    dplyr::select(-log2FC) %>%
    pivot_wider(names_from = Condition, values_from = log2FC_scaled) %>%
    drop_na() %>%
    ungroup() 
}


# data ----
## annotation ====
pfu_arcog  <- vroom(here("data/genome/arcog_pfu_table_new.txt"))
arcog_info <- vroom(here("data/genome/funclass.tab.txt"), col_names = F) %>%
  dplyr::select(-2) %>%
  dplyr::rename(category = 1, category_name = 2) %>%
  dplyr::mutate(big_category = ifelse(category %in% 1:4, category, NA),
                big_category_name = ifelse(category %in% 1:4, category_name, NA)) %>%
  fill(big_category,.direction = "down") %>%
  fill(big_category_name,.direction = "down") %>%
  dplyr::filter(!category %in% 1:4) %>%
  dplyr::rename(arcog_cat = 1, 
                arcog_cat_name = 2,
                arcog_cat_big = 3,
                arcog_cat_big_name = 4) %>%
  dplyr::mutate(plot_order = c("03","04","02","01","05",
                               "06","07","08","09","10","11","12","13","14","15","16",
                               "17","18","19","20","21","22","23","24",
                               "25","26"))


## cds ids ====
pfu_gff <- read.gff(here("data/genome/GCF_008245085.1_ASM824508v1_genomic.gff")) %>%
  dplyr::filter(type == "gene") %>%
  dplyr::mutate(locus_tag = str_split_fixed(str_split_fixed(attributes, ";old_",2)[,1],"locus_tag=",2)[,2],
                old_locus_tag = str_split_fixed(attributes, "old_locus_tag=",2)[,2],
                width = abs(start-end)) %>%
  dplyr::select(locus_tag,old_locus_tag,start, end, strand,width)

pfu_annotation <- pfu_gff %>%
  left_join(pfu_arcog, by = c("old_locus_tag" = "new")) %>%
  left_join(arcog_info, by = c("COG.x" = "arcog_cat")) 

## Read in your data as a data frame ====
comp_table_rna <- vroom(here("data/Rdata/comp_table_RNA.tsv"))
comp_table_ms  <- vroom(here("data/Rdata/comp_table_MS.tsv"))

## Condition selection ====
cs_cond <- paste0("CS", c("1","2","3","R"))
hs_cond <- paste0("HS", c("1","2","R"))

## high varianve genes ====
rna_sign_genes <- comp_table_rna %>%
  dplyr::filter(Condition %in% hs_cond) %>%
  group_by(gene) %>%
  dplyr::filter(sum(padj < 0.05) >= 1) %>%
  distinct(gene) %>%
  drop_na() %>%
  pull(gene)

ms_sign_genes <- comp_table_ms %>%
  dplyr::filter(Condition %in% hs_cond) %>%
  group_by(gene) %>%
  dplyr::filter(sum(padj < 0.05) >= 1) %>%
  distinct(gene) %>%
  drop_na() %>%
  pull(gene)

## select all significant genes from MS analysis
sign_genes <- ms_sign_genes

## Calculate the z-scores for each gene's foldchanges ====
zscore_rna_hs <- calc_z_pca_perCol(comp_table_rna, "rna", hs_cond)
zscore_ms_hs  <- calc_z_pca_perCol(comp_table_ms, "ms", hs_cond)

# analysis ----
## PCA ====
### pre-select significant genes ####
rna_hs <- zscore_rna_hs %>% dplyr::filter(gene %in% sign_genes)
ms_hs  <- zscore_ms_hs %>% dplyr::filter(gene %in% sign_genes)

in_t_hs <- left_join(rna_hs,ms_hs) %>%
  drop_na()

pca_hs <- prcomp(in_t_hs[,-1], scale = F)

hs_dist_matrix <- dist(pca_hs$x, method = "euclidean")
hs_hc <- hclust(hs_dist_matrix, method = "ward.D2")

### n of cluster? ####
fviz_nbclust(x = as.matrix(hs_dist_matrix), hcut, method = "silhouette", k.max = 8) 


## set cluster size ####
hs_clusters <- cutree(hs_hc, k = 5)

hs_n_cluster <- max(hs_clusters)

hs_clusters1 <- in_t_hs %>%
  mutate(cluster = as.factor(hs_clusters)) %>%
  pivot_longer(c(-gene,-cluster), names_to = "Condition", values_to = "log2FC") %>%
  dplyr::mutate(tech = str_split_fixed(Condition, "_", 2)[,1],
                set = str_split_fixed(Condition, "_", 2)[,2]) 

hs_clusters2 <- hs_clusters1 %>%
  group_by(cluster, set, tech) %>%
  summarise(q1 = quantile(log2FC, prob = c(.25)),
            q2 = quantile(log2FC, prob = c(.5)),
            q3 = quantile(log2FC, prob = c(.75)),
            mean = mean(log2FC))

hs_compare_set <- data.table(cluster = as.factor(1:hs_n_cluster),
                             set = rep("C",hs_n_cluster),
                             tech = rep(x = c("rna", "ms"),each = hs_n_cluster),
                             q1 = 0, q2 = 0, q3 = 0, mean = 0) %>%
  as_tibble()

# check z score plot ----
ggplot(data = bind_rows(hs_clusters2,hs_compare_set),
       aes(x = set, y = q2, group = tech, color = tech, fill = tech)) +
  facet_grid(rows = vars(cluster),
             scales = "free") +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.4, linewidth = NA) +
  geom_hline(yintercept = 0, linetype = "solid") +
  theme_light() +
  scale_y_continuous(limits = c(-1.3,1.3), expand = c(0,0)) +
  geom_line(size = 1) +
  geom_point(shape = 21, color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_manual(values = c("#EEBD51", "#A2A2A2")) +
  scale_fill_manual(values = c("#EEBD51", "#A2A2A2")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed", color = "#B7B7B7")) +
  guides(fill = "none", color = "none")


# Goseq analysis ----
## run goseq ====
### input gene cluster ####
in_go_hs <- in_t_hs %>%
  mutate(cluster = as.factor(hs_clusters)) %>%
  dplyr::select(gene, cluster)

heat_shock_cluster <- in_go_hs

### loop ####
set1 <- data.table()
hs_cluster_table <- data.table()
for(i in 1:hs_n_cluster){
  set1 <- goseq_cluster(in_go_hs, i)
  hs_cluster_table <- bind_rows(hs_cluster_table,set1)
}

## plot ====
plot_heatmap_arcog_simple(hs_cluster_table, "clusters HS", 0.05)

# cluster analysis using different metrics ----
## counts ====
load(here("data/Rdata/normalized_counts_rem_bad_table.Rdata"))
load(here("data/Rdata/normalized_ibaq_MS_remBad.Rdata"))

counts_table_wide <- bind_rows(normalized_ibaq_MS_remBad %>% 
                                 dplyr::rename(counts = 3) %>%
                                 dplyr::mutate(type = "ibaq"),
                               normalized_counts_rem_bad_table %>% 
                                 dplyr::rename(counts = 3) %>%
                                 dplyr::mutate(type = "rna")) %>%
  dplyr::mutate(counts = log10(counts)) %>%
  distinct(gene, type, Condition, .keep_all = T) %>%
  pivot_wider(names_from = type, values_from = counts) 

counts_wide_cluster <- counts_table_wide %>%
  dplyr::filter(Condition %in% c(hs_cond, "Ctrl")) %>%
  left_join(hs_clusters1 %>% distinct(gene, cluster)) %>%
  pivot_longer(ibaq:rna,names_to = "tech", values_to = "log2FC") %>%
  dplyr::rename(set = Condition) %>%
  dplyr::filter(!is.na(cluster)) %>%
  group_by(cluster, set, tech) %>%
  dplyr::filter(!is.na(log2FC)) %>%
  summarise(q1 = quantile(log2FC, prob = c(.25)),
            q2 = quantile(log2FC, prob = c(.5)),
            q3 = quantile(log2FC, prob = c(.75)),
            mean = mean(log2FC, na.rm = T))


ggplot(data = counts_wide_cluster,
       aes(x = set, y = q2, group = tech, color = tech, fill = tech)) +
  facet_grid(rows = vars(cluster),cols = vars(tech),
             scales = "free") +
  geom_hline(yintercept = 2, color = "black", linetype = "solid") +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.4, linewidth = NA) +
  scale_y_continuous(limits = c(1, 3.6)) +
  theme_light() +
  geom_line(size = 1) +
  geom_point(shape = 21, color = "black", size = 2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_manual(values = c("#EEBD51", "#A2A2A2")) +
  scale_fill_manual(values = c("#EEBD51", "#A2A2A2")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

## log2 foldchanges ====
rna_wide_cs <- fc_to_wide(comp_table_rna, "rna", hs_cond)
ms_wide_cs  <- fc_to_wide(comp_table_ms, "ms", hs_cond)

fc_table_wide <- bind_rows(comp_table_rna %>%
                             dplyr::mutate(tech = "rna"),
                           comp_table_ms %>%
                             dplyr::mutate(tech = "ms")) 


fc_wide_cluster <- fc_table_wide %>%
  dplyr::filter(gene %in% sign_genes,
                Condition %in% c(hs_cond)) %>%
  left_join(hs_clusters1 %>% distinct(gene, cluster)) %>%
  dplyr::filter(!is.na(cluster)) %>%
  dplyr::rename(set = Condition) %>%
  group_by(cluster, set, tech) %>%
  dplyr::filter(!is.na(log2FC)) %>%
  summarise(q1 = quantile(log2FC, prob = c(.25)),
            q2 = quantile(log2FC, prob = c(.5)),
            q3 = quantile(log2FC, prob = c(.75)),
            mean = mean(log2FC, na.rm = T))

ggplot(data = bind_rows(fc_wide_cluster,hs_compare_set),
       aes(x = set, y = q2, group = tech, color = tech, fill = tech)) +
  facet_grid(rows = vars(cluster),cols = vars(tech),
             scales = "free_y") +
  geom_hline(yintercept = 0, color = "black", linetype = "solid") +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.4, linewidth = NA) +
  theme_light() +
  geom_line(size = 1) +
  geom_point(shape = 21, color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_manual(values = c("#EEBD51", "#A2A2A2")) +
  scale_fill_manual(values = c("#EEBD51", "#A2A2A2")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linetype = "dashed", color = "#B7B7B7")) +
  guides(fill = "none", color = "none")
