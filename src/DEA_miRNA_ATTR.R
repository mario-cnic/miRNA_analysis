############## Installation ################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("DESeq2")
# install.packages("readxl")
# install.packages("tidyverse")



############## Preprocessing ###############
library(DESeq2, quietly = TRUE)
library(readxl)
library(tidyverse,quietly = TRUE)

counts_xlsx<-"C:/Users/lgonzalezp/Desktop/miRNA_analysis/data/12.1.miRNAExp/Readcount_TPM.xls"
counts_csv_path <-"C:/Users/lgonzalezp/Desktop/miRNA_analysis/data/readcounts.csv"
# pre_counts <- read_xls(counts_xlsx)
# No puedo leer archivo, por que
# file.exists(counts_xlsx)
# format_from_ext(counts_xlsx)
# format_from_signature(counts_xlsx)

## Create readcount file from Novogene
pre_counts2 <- read.delim(counts_xlsx)
# Variable transformation
pre_counts2 <- pre_counts2 %>% 
  select(ends_with(".readcount")) %>% 
  rename_with(~gsub(".readcount","",.x,fixed = TRUE))

write.csv(pre_counts2,counts_csv_path, row.names = FALSE)
help(write.csv)
## Create metadata file
meta_csv_path <- "C:/Users/lgonzalezp/Desktop/miRNA_analysis/data/meta.csv"
samples <- colnames(pre_counts2)[2:13]
condition <- sapply(samples,function(x) ifelse(grepl('ATTR',x,fixed = TRUE),'disease','control'))

metadata <- data.frame(condition)
metadata <- rownames_to_column(metadata,"ID")

write.csv(metadata, meta_csv_path,row.names = FALSE)

################ Diff Exp #####################
readcounts <- read_csv(counts_csv_path,) %>% column_to_rownames("sRNA")
head(readcounts)
metadata <- read_csv(meta_csv_path) %>% mutate(condition = factor(condition))
head(metadata)
deseq <- DESeqDataSetFromMatrix(countData = readcounts,
                              colData = metadata,
                              design = ~ condition)
deseq <- DESeq(deseq)
results <- results(deseq)
resultsNames(deseq) #????
results_table <- data.frame(results(deseq))

#normalized_counts <- data.frame(counts(deseq,normalized=TRUE))



####Volcano plot
library(ggplot2)

results_table$Minuslog10padj <- -log10(results_table$padj) #necesaria esta columna
# para hacer el volcano plot, que usa valores de p value ajustado en base 10

logical_na <- !is.na(results_table$padj)
filtered_results <-results_table[logical_na,]

significant <- filtered_results$padj < 0.05 #& abs(filtered_results$log2FoldChange) > 2 #coger aquellas valores significativos
# y con valor de foldchange mayor que 2, que es con lo que vamos a poner colorines como valor de significancia (is.infinite es quitar los valores infinitos)
significant_results <- filtered_results[significant,] #vector que solo contiene los valores significativos 

positive <- significant_results$log2FoldChange >0
significant_positive <- significant_results[positive,]

negative <- significant_results$log2FoldChange <0
significant_negative <- significant_results[negative,]

trans <- -log10(0.05)

significant_positive$microRNA <- row.names (significant_positive)

ggplot() + 
  geom_point(data= filtered_results, aes(x=log2FoldChange, y = Minuslog10padj)) + 
  geom_point(data = significant_positive,aes( x=log2FoldChange, y = Minuslog10padj, color=microRNA))+
  geom_point(data = significant_negative,aes( x=log2FoldChange, y = Minuslog10padj), color="red")+
  geom_line(data=filtered_results,aes(x=log2FoldChange,y=trans),color="blue")
  

lista <- arrange (significant_results, desc(abs(significant_results$log2FoldChange))) %>% 
  arrange (significant_results, significant_results$padj)

#write.table (lista, file= "microRNAs.csv", sep=",", dec=".", col.names=TRUE, row.names=TRUE )
write.csv (lista, file= "microRNAs.csv")


