library(ggplot2)
library(dplyr)
library(hrbrthemes)

F1 <- read.csv("condition_2_D_vs_1_E_gen.csv", header=TRUE, row.names=1)
F2 <- read.csv("condition_3_P_vs_2_D_gen.csv", header=TRUE, row.names=1)


F1_filtered <- F1[F1$hgnc_id%in% F2$hgnc_id,]
F1_filtered$phase <- "1"


F2_filtered <- F2[F2$hgnc_id %in% F1$hgnc_id,]
F2_filtered$phase <- "2"



Merge1 <- merge(F1_filtered, F2_filtered, by = "hgnc_id") %>% 
  mutate(diff_log2 = abs(log2FoldChange.x - log2FoldChange.y)) %>% 
  arrange(diff_log2) 

# write.csv(Merge1, "gen_compare_diff1_diff2.csv")
# 
# 
# ggplot(Merge1) +
#   geom_segment( aes(x = GENEID.x, xend = GENEID.x, y = log2FoldChange.x, yend = log2FoldChange.y)) +
#   geom_point( aes(x = GENEID.x, y = log2FoldChange.x), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
#   geom_point( aes(x = GENEID.x, y = log2FoldChange.y), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
#   coord_flip()+
#   theme_ipsum() +
#   theme(
#     legend.position = "none",
#   ) +
#   xlab("Gens") +
#   ylab("Difference in log2 counts")



##Comparison
common_elements <- intersect(k1, k2)

 
### Enrishment Analysis

gor <- enrichGO(gene = common_elements,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pvalueCutoff = 0.05,
                readable = TRUE)
g1 <- data.frame(gor)

write.csv(g1, "541GO.csv")

keggr <- enrichKEGG(gene = common_elements,
                    organism = "hsa",
                    keyType = "ncbi-geneid",
                    pvalueCutoff = 0.05)



keggr <- setReadable(keggr, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
k1 <- data.frame(keggr)
write.csv(g1, "541KEGG.csv")



pathway_list <- keggr$ID[keggr$qvalue < 0.05]
pathview(gene.data = gene_list, pathway.id = pathway, species = "hsa")





