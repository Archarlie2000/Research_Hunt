for (pathway in pathway_list) {
  pathview(gene.data = gene_list, pathway.id = pathway, species = "hsa")
}

pathview(gene.data = gene_list, pathway.id = "hsa04814", species = "hsa")

for (pathway in pathway_list[6:8]) {
  pathview(gene.data = gene_list, pathway.id = pathway, species = "hsa")
}
