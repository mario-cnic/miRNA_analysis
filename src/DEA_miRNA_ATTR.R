


############## Installation ################
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("DESeq2")
install.packages("readxl")
install.packages("tidyverse")

############## Preprocessing ###############
library(DESeq2, quietly = TRUE)
library(readxl)
library(tidyverse,quietly = TRUE)

counts_xlsx <- "~/Proyectos/miRNA_analysis/data/Novogene/X204SC24050930-Z01-F002_release/X204SC24050930-Z01-F002/02.Result_X204SC24050930-Z01-F002_Homo_sapiens/Result_X204SC24050930-Z01-F002_Homo_sapiens/Result_X204SC24050930-Z01-F002_Homo_sapiens/12.DiffExprAnalysis/12.1.miRNAExp/Readcount_TPM.xls"
counts_csv_path <- "data/readcounts.csv"
pre_counts <- read_xls(counts_xlsx)
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

## Create metadata file
meta_csv_path <- "data/meta.csv"
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
dds <- DESeqDataSetFromMatrix(countData = readcounts,
                              colData = metadata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
tail(res)
resultsNames(dds)
