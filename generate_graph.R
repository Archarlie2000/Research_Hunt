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

tmp = sapply(keggresids, function(pid) pathview(gene.data = LogChange, pathway.id = pid, species="hsa"))