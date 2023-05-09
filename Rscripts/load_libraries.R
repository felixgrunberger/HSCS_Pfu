# > libraries and functions needed for the analysis < #

# libraries ----
packages <- c("tidyverse", "ape", "data.table", "Rsubread", "readxl", "vroom", "limma","tictoc","ggsci","openxlsx",
              "DESeq2", "ggforce", "Biostrings", "ggseqlogo", "cowplot", "scico", "patchwork","ggridges",
              "ggalluvial", "janitor", "RColorBrewer", "colorspace", "ggpubr", "cowplot", "reshape2")
invisible(lapply(packages, require, character.only = TRUE))

# colors ----
batlow_custom <- c(scico::scico(8, palette = 'batlow')[c(1,3:5)],
                   "grey70",
                   scico::scico(8, palette = 'batlow')[6:8])

batlow_custom_re <- rev(c("grey70",
                                  scico::scico(8, palette = 'batlow')[c(1,3:5)],
                                  scico::scico(8, palette = 'batlow')[6:8]))
BrBG5 <- colorspace::divergingx_hcl(n = 101, palette = "BrBG", rev = F)[c(10,30,50,70,90)]
BrBG7 <- colorspace::divergingx_hcl(n = 7, palette = "BrBG", rev = F)
BrBG11 <- colorspace::divergingx_hcl(n = 11, palette = "BrBG", rev = F)


# functions ----
read_in_gff <- function(input_file){
  read.gff(input_file) %>%
    dplyr::filter(!type %in% c("exon", "gene", "region", "origin of replication")) %>%
    as_tibble() %>%
    dplyr::mutate(start_feature = start, end_feature = end,strand_feature = strand) %>%
    dplyr::mutate(Parent = str_split_fixed(str_split_fixed(attributes, ";Parent=",2)[,2],";Dbxref",2)[,1],
                  ecogene = str_split_fixed(str_split_fixed(attributes, ",GeneID", 2)[,1], "EcoGene:",2)[,2],
                  short_gene = str_split_fixed(str_split_fixed(attributes, ";locus_tag", 2)[,1], "gene=",2)[,2],
                  id_name = ifelse(type %in% "repeat_region", str_split_fixed(str_split_fixed(attributes, ";Note=", 2)[,1], "ID=", 2)[,2],
                                   ifelse(type %in% "pseudogene", str_split_fixed(str_split_fixed(attributes, ";Dbxref=", 2)[,1], "ID=", 2)[,2],
                                          ifelse(type %in% "sequence_feature", str_split_fixed(str_split_fixed(attributes, ";Dbxref=", 2)[,1], "ID=", 2)[,2],
                                                 ifelse(type %in% "mobile_genetic_element", str_split_fixed(str_split_fixed(attributes, ";gbkey=", 2)[,1], "ID=", 2)[,2],
                                                        str_split_fixed(str_split_fixed(attributes, ";Parent=", 2)[,1], "ID=", 2)[,2])))),
                  locus_name = ifelse(type %in% c("CDS","mobile_genetic_element", "ncRNA", "recombination_feature"), str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], ";", 2)[,1],
                                      ifelse(type ==  "pseudogene",  str_split_fixed(str_split_fixed(attributes, ";gene_biotype", 2)[,1], "gene=", 2)[,2],
                                             ifelse(type == "repeat_region", str_split_fixed(str_split_fixed(attributes, ";gbkey", 2)[,1], "Note=", 2)[,2],
                                                    ifelse(type %in% "sequence_feature", str_split_fixed(str_split_fixed(attributes, ";locus_tag=", 2)[,1], "gene=", 2)[,2],
                                                           ifelse(type %in% "mobile_genetic_element", str_split_fixed(attributes, "insertion sequence:", 2)[,2],
                                                                  ifelse(type == "rRNA", str_split_fixed(str_split_fixed(attributes, ";product=", 2)[,2], " ", 2)[,1], 
                                                                         ifelse(type == "tRNA", str_split_fixed(attributes, ";product=", 2)[,2], NA ))))))),
                  width = abs(start_feature - end_feature)) %>%
    dplyr::select(seqid, id_name, locus_name, start_feature, end_feature, strand_feature, Parent, type, width, ecogene, short_gene) %>%
    mutate(gene = str_split_fixed(Parent,"-",2)[,2])
}