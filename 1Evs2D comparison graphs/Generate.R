

for (pathway in pathway_list[1:2]) {
  pathview(gene.data = gene_list, pathway.id = pathway, species = "hsa")
}


pathview(gene.data = gene_list, pathway.id = "hsa01250", species = "hsa")

for (pathway in pathway_list[4:24]) {
  pathview(gene.data = gene_list, pathway.id = pathway, species = "hsa")
}
