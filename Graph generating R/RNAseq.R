library(tximport)
library(readr)
library(DESeq2)
library(apeglm)
library(Glimma)

#Set the directory with your files
dir <- "~/Hunt Lab/Mink Transcriptome/Results/genes"

#Select the samples you will run
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)

#Deletes the samples that will not be used in the *pairwise* comparison
samples <- samples[!grepl("3_P", samples$condition),] #Estrus Diapause Comparison
#samples <- samples[!grepl("1_E", samples$condition),] #Diapause Pregnancy Comparison

#Imports the sample.genes.results from RSEM
files <- file.path(dir,samples$experiment)

#say the names of your samples
names(files) <- samples$sample

#import using tximport
tx2gene <- read.table(file.path(dir, "tx2geneMink.txt"), sep = '\t', header = TRUE)
txi.rsem <- tximport(files, type = "rsem", tx2gene = tx2gene)

#remove genes with 0 length
zero_length = (apply(txi.rsem$length, 1, min) == 0)
txi.rsem$length = txi.rsem$length[!zero_length,]
txi.rsem$abundance = txi.rsem$abundance[!zero_length,]
txi.rsem$counts = txi.rsem$counts[!zero_length,]

#transform tximport object to DESeq2 object
sampleTable <- data.frame(condition = factor(samples$condition))
rownames(sampleTable) <- colnames(txi.rsem$counts)
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)

#start DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")

#Give gene name to all transcripts
tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME),]
csvfile <- merge(as.data.frame(res), tx2gene, by.x="row.names", by.y= "TXNAME", all.x=TRUE, all.y = FALSE)
writedir <- file.path(dir,paste(resultsNames(dds)[2],".csv", sep = ""))
write.csv(csvfile, writedir, col.names = TRUE)

