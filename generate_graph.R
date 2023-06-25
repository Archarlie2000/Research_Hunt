

library(gage)
library(gageData)
library(pathview)

ENRTEZID_ID <- mapIds(org.Hs.eg.db, keys = Gene_ID, keytype = "GENENAME" , column = "ENTREZID")
names(FoldChanges) = ENRTEZID_ID
head(FoldChanges)
keggrespathways = data.frame(id = rownames(kegres$greater), kegres$greater) %>%
  tibble::as.tibble() %>%
  filter(row_number() <=35) %>%
  .$id %>%
  as.character()
keggrespathways




keggres = gage(exprs = LogChange, gsets = kegg.sets.hm, same.dir = TRUE)


tmp = sapply("hsa04914", function(pid) pathview(gene.data = logchange1, pathway.id = pid, species="hsa"))

colnames(genlist_1) <- c("X", "gene_list")
colnames(genlist_2) <- c("X", "gene_list")


logchange1 <- genlist_1$gene_list
logchange2 <- genlist_2$gene_list

genID1 <- genlist_1$X
genID2 <- genlist_2$X

names(logchange1) <- genID1
names(logchange2) <- genID2

#############################

tmp = sapply("04914", function(pid) pathview(gene.data = logchange2, pathway.id = pid, species="hsa"))

tmp = sapply("04914", function(pid) pathview(gene.data = logchange1, pathway.id = pid, species="hsa"))
