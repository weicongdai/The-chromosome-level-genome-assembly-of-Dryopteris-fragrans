library(tidyverse)
library(ggsci)
library(rstatix)
# 各 subgroup 的可变剪切数量
sub_as <- read.delim("D:/Desktop/data/gene_info.txt",sep = "\t") %>%
  select(subgroup, Num.of.mRNA) %>%
  mutate(
    subgroup = case_when(subgroup == "short-intron genes with low TE proportion" ~ "SLG",
                         subgroup == "short-intron genes with high TE proportion" ~ "SHG",
                         subgroup == "long-intron genes with TE insertion" ~ "LG",
                         subgroup == "short-intron TE-free genes" ~ "SG",
                         .default = subgroup)
  ) %>%
  filter(!is.na(Num.of.mRNA)) %>%
  mutate(
    subgroup = factor(subgroup, levels = c("LG","SHG","SLG","SG"))
  )

# 绘图
pdf("D:/Desktop/data/16.AS/sub_as.pdf",height = 3,width = 4)
ggplot(sub_as,aes(x = subgroup, y = Num.of.mRNA)) +
  geom_boxplot(aes(fill = subgroup),width = 0.7, outliers = F, show.legend = F) +
  scale_fill_npg() +
  labs(x = NULL,y = "AS Number") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, colour = "black", face = "bold"),
    axis.title = element_text(size = 15, colour = "black", face = "bold")
  )
dev.off()
#--------------------------------------------------------------------------------
# 各 subgroup 的发生AS的数量
sub_as_stat <- sub_as %>%
  mutate(as_type = if_else(Num.of.mRNA == 1, true = "NASG", false = "ASG")) %>%
  mutate(as_type = factor(as_type, levels = c("NASG","ASG"))) %>%
  group_by(subgroup, as_type) %>%
  summarise(
    num = n()
  )
# 绘图
pdf("D:/Desktop/data/16.AS/sub_as_stat.pdf",height = 3,width = 4)
ggplot(sub_as_stat,aes(x = subgroup,y = num)) +
  geom_col(aes(fill = as_type),position = "fill") +
  scale_fill_npg() +
  labs(x = NULL,y = "Ratio") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, colour = "black", face = "bold"),
    axis.title = element_text(size = 15, colour = "black", face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 9, colour = "black", face = "bold")
  )
dev.off()
#--------------------------------------------------------------------------------
# 各 subgroup 的内含子数量
sub_intron <- read.delim("D:/Desktop/data/gene_info.txt",sep = "\t") %>%
  select(subgroup, Num.of.Intron) %>%
  mutate(
    subgroup = case_when(subgroup == "short-intron genes with low TE proportion" ~ "SLG",
                         subgroup == "short-intron genes with high TE proportion" ~ "SHG",
                         subgroup == "long-intron genes with TE insertion" ~ "LG",
                         subgroup == "short-intron TE-free genes" ~ "SG",
                         .default = subgroup)
  ) %>%
  filter(!is.na(Num.of.Intron)) %>%
  mutate(
    subgroup = factor(subgroup, levels = c("LG","SHG","SLG","SG"))
  )

# 绘图
pdf("D:/Desktop/data/16.AS/sub_intron.pdf",height = 3,width = 4)
ggplot(sub_intron,aes(x = subgroup, y = Num.of.Intron)) +
  geom_boxplot(aes(fill = subgroup),width = 0.7, outliers = F, show.legend = F) +
  scale_fill_npg() +
  labs(x = NULL,y = "Intron Number") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, colour = "black", face = "bold"),
    axis.title = element_text(size = 15, colour = "black", face = "bold")
  )
dev.off()
#--------------------------------------------------------------------------------