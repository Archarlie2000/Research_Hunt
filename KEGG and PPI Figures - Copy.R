library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(STRINGdb)

##### DATA PREPARATION #####

#Set the directory with your files

# Load dataset


rna_results <- read.csv("condition_3_P_vs_2_D.csv", header=TRUE, row.names=1)

# Add HGNC and uniprot names to the database
ensembl_m <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                   dataset = "nvison_gene_ensembl")

hgnc_m <- getBM(filters = "ensembl_gene_id",
               attributes = c("ensembl_gene_id","hgnc_id"),
               values = rna_results$Row.names, 
               mart = ensembl_m)

rna_results <- merge(rna_results, hgnc_m, by.x = "Row.names", by.y = "ensembl_gene_id")
# rna_results <- merge(rna_results, uniprot, by.x = "Row.names", by.y = "ensembl_gene_id")

# Drop all empty rows and anything not significant
rna_results <- rna_results[rna_results$pvalue < 0.05,]
rna_results <- rna_results[!apply(rna_results == "", 1, any), ,]
rna_results <- rna_results[!is.na(rna_results$pvalue),]
rna_results <- rna_results[!is.na(rna_results$hgnc_id),]


write.csv(rna_results, "condition_3_P_vs_2_D_gen.csv")

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

##### GO ENRICHMENT #####

# Create GO enrichment analysis (over representation)
gor <- enrichGO(gene = names(gene_list),
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pvalueCutoff = 0.05,
                readable = TRUE)

# Save the GO enrichment analysis

write.csv(as.data.frame(gor), "condition_3_P_vs_2_D_geneontology.csv")


##### KEGG ENRICHMENT #####

# Create KEGG pathway gene enrichment analysis and make the file readable
# options(clusterProfiler.download.method = "wininet") #If KEGG isn't receiving data from web run this command
keggr <- enrichKEGG(gene = names(gene_list),
                 organism = "hsa",
                 keyType = "ncbi-geneid",
                 pvalueCutoff = 0.05)
keggr <- setReadable(keggr, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# Save the KEGG pathway gene enrichment analysis

write.csv(as.data.frame(keggr), "condition_3_P_vs_2_D_kegg.csv")

##### KEGG VISUAL PATHWAY #####

# Get the list of significant KEGG pathways
pathway_list <- keggr$ID[keggr$qvalue < 0.05]

# Prepare directory for pathway output
dir.create(file.path(dir,paste0("pathways_",rna_file)), recursive = TRUE)
setwd(file.path(dir,paste0("pathways_",rna_file)))

# Make KEGG pathways using the gene_list expression values from significant pathways
for (pathway in pathway_list) {
  pathview(gene.data = gene_list, pathway.id = pathway, species = "hsa")
}

##### PROTEIN PROTEIN INTERATION GRAPH #####

#
string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=200, network_type="functional", input_directory="")
example1_mapped <- string_db$map( rna_list, "entrezgene_id", removeUnmapped = TRUE )
