

library(gage)
library(gageData)
library(pathview)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(STRINGdb)
options(repos = BiocManager::repositories())


get_genelist <- function(rna_results){

  # Add HGNC and uniprot names to the database
  ensembl_m <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                       dataset = "nvison_gene_ensembl")
  hgnc_m <- getBM(filters = "ensembl_gene_id",
                  attributes = c("ensembl_gene_id","hgnc_id"),
                  values = rna_results$Row.names, 
                  mart = ensembl_m)
  
  
  # Merge gen names to gen ID
  rna_results <- merge(rna_results, hgnc_m, by.x = "Row.names", by.y = "ensembl_gene_id")
  
  # Drop all empty rows and anything not significant
  # Drop all the wierd stuffs like N/A or blanks
  rna_results <- rna_results[rna_results$pvalue < 0.05,]
  rna_results <- rna_results[!apply(rna_results == "", 1, any), ,]
  rna_results <- rna_results[!is.na(rna_results$pvalue),]
  rna_results <- rna_results[!is.na(rna_results$hgnc_id),]
  
  #write.csv(rna_results, "condition_2_D_vs_1_E_gen.csv")
  
  
  # Convert HGNC IDs into EntrezGene IDs
  ensembl_h <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                       dataset = "hsapiens_gene_ensembl")
  
  hgnc_h <- getBM(filters = "hgnc_id",
                  attributes = c("hgnc_id", "entrezgene_id"),
                  values = rna_results$hgnc_id, 
                  mart = ensembl_h)
  
  hgnc_duplicated <- hgnc_h[duplicated(hgnc_h$hgnc_id),]
  
  hgnc_h <- hgnc_h[!duplicated(hgnc_h$hgnc_id),]
  
  rna_results <- merge(rna_results, hgnc_h, by.x = "hgnc_id", by.y = "hgnc_id")
  
  # Create an ordered list of the rna_results
  rna_list <- rna_results[order(rna_results$log2FoldChange, decreasing = TRUE),]
  rna_list <- rna_list[!is.na(rna_list$entrezgene_id),c("entrezgene_id", "log2FoldChange")]
  
  # Distill rna_results list into simplified gene_list
  gene_list <- rna_list$log2FoldChange
  names(gene_list) <- rna_list$entrezgene_id
  return(gene_list)
}


#############################
selectedgene <- "04913"


rna_results <- read.csv("condition_2_D_vs_1_E.csv", header=TRUE, row.names=1)
logchange1 <- get_genelist(rna_results)


OG <- getwd()
dir.create("E2D")
dir.create("D2P")

setwd(paste(getwd(), "E2D", sep = "/"))
tmp = sapply(selectedgene, 
             function(pid) pathview(gene.data = logchange1, 
                                    pathway.id = pid, 
                                    species="hsa"))

setwd(paste(getwd(), "D2P", sep = "/"))
tmp = sapply(selectedgene, 
             function(pid) pathview(gene.data = logchange2, 
                                    pathway.id = pid, 
                                    species="hsa"))

setwd(getwd())




setwd(paste(getwd(), "my_folder", sep = "/"))



