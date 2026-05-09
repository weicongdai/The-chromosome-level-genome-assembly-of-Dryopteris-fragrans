library(tidyverse)
library(ggsci)

pro_TE <- read.delim("pro_TE.inter",header = F) %>%
  dplyr::select(1,8,4,9) %>%
  set_names(c("Chr","GeneID","TE_type","TE_len")) %>%
  dplyr::mutate(TE_type = case_when(
    TE_type %in% "Gypsy_LTR_retrotransposon" ~ "LTR/Gypsy",
    TE_type %in% "Copia_LTR_retrotransposon" ~ "LTR/Copia",
    TE_type %in% c("DNA", "DNA?") ~ "DNA transposon",
    TE_type %in% c("L1_LINE_retrotransposon", "L2_LINE_retrotransposon", 
                   "RTE_LINE_retrotransposon", "LINE_element", "CR1_LINE_retrotransposon") ~ "LINE",
    TE_type %in% "helitron" ~ "helitron",
    TE_type %in% "SINE_element" ~ "SINE",
    TE_type %in% c("pararetrovirus", "Bel_Pao_LTR_retrotransposon", "LTR_retrotransposon") ~ "LTR/Other",
    TE_type %in% "Penelope_retrotransposon" ~ "Penelope_retrotransposon",  # 或其他类别
    TRUE ~ "Other"
  )) %>%
  dplyr::mutate(GeneID = str_replace_all(.$GeneID, "TU", "model"))

pro_TE_stat <- pro_TE %>%
  group_by(GeneID) %>%
  summarise(len = sum(TE_len))

gene_info <- read.delim("D:/Desktop/data/gene_info.txt") %>%
  dplyr::select(c(-11)) %>%
  left_join(., pro_TE_stat, by = c("id" = "GeneID")) %>%
  relocate(len, .after = Ratio)

colnames(gene_info)[11] <- "pro_TElen"

gene_info <- gene_info %>%
  mutate(pro_TElen = if_else(is.na(pro_TElen),0,pro_TElen))

gene_info <- gene_info %>%
  mutate(subgroup = case_when(
    subgroup == "short-intron genes with low TE proportion" ~ "SLG",
    subgroup == "short-intron genes with high TE proportion" ~ "SHG",
    subgroup == "long-intron genes with TE insertion" ~ "LG",
    subgroup == "short-intron TE-free genes" ~ "SG",
    .default = subgroup
  ))
write.table(gene_info, file = "D:/Desktop/data/gene_info.txt", sep = "\t", 
            col.names = T, row.names = F, quote = F)

pro_TE_stat <- pro_TE %>%
  group_by(GeneID, TE_type) %>%
  summarise(total_len = sum(TE_len)) %>%
  left_join(., dplyr::select(gene_info, id, subgroup), 
            by = c("GeneID" = "id"))

write.table(pro_TE_stat, file = "D:/Desktop/data/08.Pro/v2/pro_TE_stat.txt", sep = "\t", 
            col.names = T, row.names = F, quote = F)

pdf("D:/Desktop/data/08.Pro/v2/pro_TE_stat.pdf",width = 12,height = 4)
ggplot(filter(pro_TE_stat, !is.na(subgroup)), 
       aes(x = subgroup, y = total_len)) +
  geom_boxplot(aes(fill = subgroup),show.legend = F) +
  facet_wrap(.~TE_type,ncol = 4) +
  scale_fill_npg() +
  labs(x = NULL,y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12,colour = "black")
  )
dev.off()