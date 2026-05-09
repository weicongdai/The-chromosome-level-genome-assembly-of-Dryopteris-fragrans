library(clusterProfiler)
library(org.Df.eg.db)

gene_id <- read.delim("expansion.geneid") %>%
  pull(ID)


gene <- read.delim("all_gene.id",header = F) %>%
  pull(var = 1)

ego <- enrichGO(gene = gene_id,
                universe = gene, #背景基因集,即所有基因
                OrgDb = org.Df.eg.db,
                keyType = "GID",
                ont = "ALL",  #MF,CC,BP
                pvalueCutoff = 0.05, #若富集结果太少,则改为1
                qvalueCutoff = 0.05,
                readable = FALSE)

ego <- simplify(ego, cutoff = 0.7, #去除相似度高于0.7的GO
                by = "p.value",
                select_fun = min)
ego <- as.data.frame(ego)

write.csv(ego,file = "D:/Desktop/data/12.系统发育与基因家族进化/expansion.GO.csv",
          quote = F, row.names = F)

expan_go <- read.csv("expan_GO_filter.csv")

pdf("expan_GO.pdf",height = 8,width = 8)
ggplot(data = expan_go, aes(x = RichFactor,y = Description)) +
  geom_point(shape = 21, aes(size = Count, fill = pvalue)) +
  labs(y = NULL) +
  scale_fill_gradientn(
    colors = c("#1E88E5", "#64B5F6", "#FFFFFF", "#FFC107", "#FF5722")) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10,color = "black")
  )
dev.off()