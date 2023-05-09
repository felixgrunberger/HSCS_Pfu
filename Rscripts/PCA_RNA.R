### normalization ####
vsd <- vst(dds, blind=FALSE)
ad <- plotPCA(vsd,intgroup = c("Condition", "bio_rep"))
ggad <- as_tibble(ad$data) 

ggplot(ggad, 
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
  xlab(ad$labels$x) +
  ylab(ad$labels$y)

### remove replicates after PCA ####
remove_rna_samples <- data.table(Condition = c("HS1", "HS2", "HS3_rec", "CS4_rec"),
                                 bio_rep = c(1, 4, 1, 2),
                                 remove = T) %>%
  left_join(sample_ids) %>%
  dplyr::select(sample_ids) %>%
  deframe()

## Start again from PCA ====
colDataRaw2 <- data.table(sample_ids = sample_design) %>%
  left_join(sample_ids) %>%
  dplyr::filter(!is.na(Condition),
                !sample_ids %in% remove_rna_samples) 

## DESeq object ====
### calc DESeq2 ####
dds2 <- DESeq2::DESeqDataSetFromMatrix(countData = counts[,colnames(counts) %in% colDataRaw2$sample_ids], 
                                       colData = colDataRaw2,
                                       design = ~Condition)

dds2 <- dds2[rowSums(counts(dds2)) > 1,]
dds_adjusted2 <- estimateSizeFactors(dds2)

### get normalized counts table ####
normalized_counts2 <- counts(dds_adjusted2, normalized = T)
vroom_write(normalized_counts2 %>%
              as.data.frame() %>%
              rownames_to_column("gene"), 
            here("data/Rdata/deseq_normalized_counts.tsv"))

### normalization ####
vsd2 <- vst(dds2, blind = FALSE)
ad2 <- plotPCA(vsd2,intgroup = c("Condition", "bio_rep"))
ggad2 <- as_tibble(ad2$data) 

pdf(here("figures/Plots/230307_PCA_RNA.pdf"),width = 6, height = 4)
ggplot(ggad2, 
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
  scale_x_continuous(expand = c(.2,0)) +
  scale_y_continuous(expand = c(.2,0)) +
  theme(panel.grid.minor = element_blank()) +
  xlab(ad2$labels$x) +
  ylab(ad2$labels$y) 
dev.off()

### save unnormalized counts ####
counts_remBad <- counts[,colnames(counts) %in% colDataRaw2$sample_ids]
save(counts_remBad,
     file = here("data/Rdata/featurecounts_rnaseq_remBad.Rdata"))

### save unnormalized counts ####
featurecounts_rem_bad_table <- counts_remBad %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "set", values_to = "counts") %>%
  left_join(sample_ids, by = c("set" = "sample_ids")) %>%
  left_join(pfu_gff_cds, multiple = "all") %>%
  group_by(Condition, gene) %>%
  summarise(Mean_counts = mean(counts, na.rm = T))
save(featurecounts_rem_bad_table,
     file = here("data/Rdata/featurecounts_rem_bad_table.Rdata"))

### save normalized counts ####
normalized_counts_rem_bad_table <- normalized_counts2 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "set", values_to = "counts") %>%
  left_join(sample_ids, by = c("set" = "sample_ids")) %>%
  left_join(pfu_gff %>%
              dplyr::rename(gene = locus_tag), multiple = "all") %>%
  mutate(rpk = counts/(width/1000)) %>%
  group_by(set) %>%
  mutate(scaling_factor = sum(rpk, na.rm = T)/1000000) %>%
  ungroup() %>%
  mutate(tpm = rpk/scaling_factor) %>%
  group_by(Condition, gene) %>%
  summarise(Mean_tpm = mean(tpm, na.rm = T))
save(normalized_counts_rem_bad_table,
     file = here("data/Rdata/normalized_counts_rem_bad_table.Rdata"))