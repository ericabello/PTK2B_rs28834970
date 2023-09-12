library("RColorBrewer")
library("pheatmap")
library("DESeq2")
library("ggplot2")
library("genefilter")
library(tidyverse)
library(BiocParallel)
register(MulticoreParam(4))
options(stringsAsFactors = F)
#library("org.Hs.eg.db")
library("dplyr")
#library(DiffBind)
#library("clusterProfiler")
#library("ChIPseeker")
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(EnhancedVolcano)


setwd("~/Documents/data_PTK2B_microglia/CR_June 2022/analysis_CR_CEBPB_101022/raw_counts_101022/")
#ANALYSIS BEDCOVERAGE COUNTS + DESEQ2 BOTH BATCHES. excluded F4 as per RNAseq
#merge df counts and rename columns to input into dds object (as countdata)
initial_filename <- "micro_PTK2B_B11_CR_C_CR_CEBPb_micro_PTK2B.counts";
full_df <- read.table(paste(initial_filename, sep="/") , header = F, stringsAsFactors = F, sep = "\t")
colnames(full_df) <- c("chr", "start", "end", initial_filename)
full_df$position <- paste(paste(full_df$chr, full_df$start, sep=":"), full_df$end, sep="-");
full_df <- full_df[c("position", initial_filename)]
for (filename in c("micro_PTK2B_K1_CR_C_CR_CEBPb_micro_PTK2B.counts","micro_PTK2B_K2_CR_C_CR_CEBPb_micro_PTK2B.counts","micro_PTK2B_C11_CR_C_CR_CEBPb_micro_PTK2B.counts","micro_PTK2B_D3_CR_C_CR_CEBPb_micro_PTK2B.counts","micro_PTK2B_C8_CR_C_CR_CEBPb_micro_PTK2B.counts")) {
  counts_df <- read.table(paste(filename, sep="/") , header = F, stringsAsFactors = F, sep = "\t") 
  full_df[filename] = counts_df[,4]
}

head(full_df)

rownames(full_df) <- full_df[,1]
colnames(full_df)[2:7] <- c("WT_B11","WT_K1","WT_K2","HOM_C11","HOM_D3","HOM_C8")
full_df <- full_df[2:7]
#filtered out peaks that had readcount less than 10
counts_filter <- function(x, n) { expr <- full_df %>% filter(rowMeans(.) > n)
return (expr)}
full_df <- counts_filter(full_df, 10)
head(full_df)

#make coldata and dds objects
coldata <- data.frame(row.names = colnames(full_df), genotype_PTK2B=as.factor(c(rep("WT",3) ,rep("HOM", 3))))
coldata

dds <- DESeqDataSetFromMatrix(countData = full_df,
                              colData = coldata,
                              design = ~genotype_PTK2B)
#estimate size factors
dds <- estimateSizeFactors(dds)
#tried without filtering for low read counts because peak calling was already filtering step and should have selected peaks with high read count (aka significant peaks). filter dds by low read counts. 
#idx <- rowSums( counts(dds, normalized=TRUE) >= 20 ) >= 3
#dds <- dds[idx,]

dds$genotype_PTK2B <- relevel(dds$genotype_PTK2B, ref = "WT")

dds <- DESeq(dds)

#calculate normalised counts
counts_normalised <- as.data.frame(counts(dds, normalized=TRUE))
counts_normalised$peak <- rownames(counts_normalised)
write.table(counts_normalised, file="../normalised_counts_CR_CEBP_micro.csv", quote=F, row.names=F)

#extract results
res <- results(dds)
summary(res)


#save all res
write.csv(as.data.frame(res), 
          file="../CR_CEBP_deseq2_results_micro.csv")

res_up <- subset(res,res$log2FoldChange>0.5 & padj<0.1)
res_up
write.csv(res_up, file="../res_up_0.5_0.05_deseq2_CR_CEBP_micro_PTK2BPeak.csv")
#res_down <- subset(res,res$log2FoldChange<(-0.5) & res$pvalue<0.05 )
#res_down
#write.csv(res_down, file="../res_down_0.5_0.05_deseq2_CR_CEBP_micro.csv")


#volcano plot to show the only statistically significant peak is PTK2B variant peak
res_plot <- res1 %>% mutate(diffExpr = case_when(padj < 0.05 & log2FoldChange > 0.5 ~ "UP", padj < 0.05 & log2FoldChange < -0.5 ~ "DOWN")) %>% replace_na(list(diffExpr="Not significant")) 
res_plot <- res_plot %>% mutate(delabel = case_when(res_plot$diffExpr != "Not significant" ~ rownames(res_plot)))
ggplot(data=res_plot, aes(x=log2FoldChange, y=-log10(padj), label=delabel)) + geom_point(size=2,aes(col=diffExpr)) + xlim(-4,4) + ylim(0,2) + scale_color_manual(values=c("blue", "red")) +theme_minimal()+geom_text_repel()
ggsave('../volcanoPlot_CEBP_CR_micro_biggerpoints.jpeg',device= "jpeg")
