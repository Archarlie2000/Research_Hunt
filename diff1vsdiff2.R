df1_KEGG <- read.csv("condition_2_D_vs_1_E_kegg.csv", header = TRUE)[, -1]
df2_KEGG <- read.csv("condition_3_P_vs_2_D_kegg.csv", header = TRUE)[, -1]

df1_GEN <- read.csv("condition_2_D_vs_1_E_gen.csv", header = TRUE)[, -1]
df2_GEN <- read.csv("condition_3_P_vs_2_D_gen.csv", header = TRUE)[, -1]


common_elements <- intersect(df1_KEGG$ID, df2_KEGG$ID)
common_df <- df1_KEGG[common_df$ID %in% common_elements, ]



target_gene <- sapply(strsplit(common_df$geneID, "/"), as.list)
df1_GEN_filtered <- df1_GEN[df1_GEN$GENEID %in% target_gene[[1]], ]
df2_GEN_filtered <- df2_GEN[df1_GEN$GENEID %in% target_gene[[1]], ]

Holistic_view <-  merge(df1_GEN_filtered, df2_GEN_filtered, by = "GENEID")




list_1 <- read.csv("genlist_1.csv")
list_2 <- read.csv("genlist_2.csv")
