library("RColorBrewer")
library("pheatmap")
library("DESeq2")
library("ggplot2")
library("genefilter")
library(tidyverse)
library(BiocParallel)
register(MulticoreParam(4))
options(stringsAsFactors = F)
library("org.Hs.eg.db")
library("dplyr")
#library(DiffBind)
#library("clusterProfiler")
#library("ChIPseeker")
#library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("limma")


setwd("/Volumes/Samsung_T5/Sanger_laptop_300823/Documents/ptk2b_kolf_B11_C8_D3/ATACseq/final_counts_Nov2022/")
#ANALYSIS BEDCOVERAGE COUNTS + DESEQ2 BOTH BATCHES. excluded F4 as per RNAseq
#merge df counts and rename columns to input into dds object (as countdata)
initial_filename <- "B11_1_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts";
full_df <- read.table(paste(initial_filename, sep="/") , header = F, stringsAsFactors = F, sep = "\t")
colnames(full_df) <- c("chr", "start", "end", initial_filename)
full_df$position <- paste(paste(full_df$chr, full_df$start, sep=":"), full_df$end, sep="-");
full_df <- full_df[c("position", initial_filename)]
for (filename in c("B11_2_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","B11_3_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","KOLF2_1_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","KOLF2_2_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","KOLF2_3_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts", "C11_1_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","C11_2_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","C11_3_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts", "D3_1_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","D3_2_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","D3_3_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","C8_1_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","C8_2_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts","C8_3_macro_PTK2B_ATAC_ATAC_macro_PTK2B.counts")) {
  counts_df <- read.table(paste(filename, sep="/") , header = F, stringsAsFactors = F, sep = "\t") 
  full_df[filename] = counts_df[,4]
}

head(full_df)

rownames(full_df) <- full_df[,1]
colnames(full_df)[2:16] <- c("WT_B11_1","WT_B11_2","WT_B11_3"," KOLF2_1"," KOLF2_2"," KOLF2_3", "HOM_C11_1","HOM_C11_2","HOM_C11_3", "HOM_D3_1","HOM_D3_2","HOM_D3_3","HOM_C8_1","HOM_C8_2","HOM_C8_3")
full_df <- full_df[2:16]
#make coldata and dds objects
coldata <- data.frame(row.names = colnames(full_df), genotype_PTK2B=as.factor(c(rep("WT",6) ,rep("HOM", 9))), batch=as.factor(c(rep("2", 6), rep("1", 3), rep("2", 6))))
coldata
dds <- DESeqDataSetFromMatrix(countData = full_df,
                              colData = coldata,
                              design = ~ batch+ genotype_PTK2B)
#estimate size factors
dds <- estimateSizeFactors(dds)
#tried without filtering for low read counts because peak calling was already filtering step and should have selected peaks with high read count (aka significant peaks). filter dds by low read counts. 
#idx <- rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 3
#dds <- dds[idx,]

dds$genotype_PTK2B <- relevel(dds$genotype_PTK2B, ref = "WT")

dds <- DESeq(dds)

#calculate normalised counts
counts_normalised <- as.data.frame(counts(dds, normalized=TRUE))
counts_normalised$peak <- rownames(counts_normalised)
write.table(counts_normalised, file="normalised_counts_ATAC_macro.csv", quote=F, row.names=F)

#extract results
res <- results(dds)
summary(res)
#save all res
write.csv(as.data.frame(res), 
          file="ATAC_deseq2_results_macro.csv")

res_up <- subset(res,res$log2FoldChange>0.5 & padj<0.05 )
res_up
write.csv(res_up, file="res_up_0.5_0.05_deseq2_ATAC_macro.csv")
res_down <- subset(res,res$log2FoldChange<(-0.5) & res$padj<0.05 )
res_down
write.csv(res_down, file="res_down_0.5_0.05_deseq2_ATAC_macro.csv")



legend <- coldata
legend

#batch correction on vsd tranformed counts
vsd <- vst(dds, blind = FALSE)
#for custom pca plot return df of PC1 and 2
df <- plotPCA(vsd, "genotype_PTK2B", returnData = TRUE)

a <- c("1","2","3","4", "5")
df$hiPSC_clone <- c(rep(a, each=3))
#custom PCA plot
#ggplot(data= df, mapping = aes(x=PC1 , y=PC2))+geom_point(size=4, mapping = aes(color=genotype_PTK2B, shape=hiPSC_clone))+theme_classic()+theme(aspect.ratio = 1)

#to change shapes of each sample manually and % PC1 and 2 labels!!!
ggplot(data= df, mapping = aes(x=PC1 , y=PC2))+geom_point(size=4, mapping = aes(color=genotype_PTK2B, shape=hiPSC_clone))+ scale_color_manual(values=c("orange", "blue"))+scale_shape_manual(values=c(17, 18, 16,19, 15,20))+theme_classic()+theme(aspect.ratio = 1,legend.title = element_text(color = "black", size = 8),legend.text = element_text(color = "black", size = 6), legend.key.size = unit(0.3, "cm"),legend.margin=margin(0,0,0,0),
                                                                                                                                                                                                                                                   legend.box.margin=margin(0,0,0,0), legend.box.spacing = margin(10,10,10,10))+ xlab("PC1 80% variance") + ylab("PC2 10% variance")
#ggsave("/Users/eb19/Documents/ptk2b_kolf_B11_C8_D3/RNAseq/all_samples_noF4/PCA_vsdBatchCorr_macro_noNames.jpeg", device= "jpeg")

ggsave('PCA_ATAC_vsdBatchCorr_macro_noNames.pdf', device = "pdf")

#annotate peaks down 0.5
library("clusterProfiler")
library("ChIPseeker")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
as.data.frame(rownames(res_down)) -> df
separate(df, 1, c("chr", "start"), sep = ":") -> df
separate(df, 2, c("start", "end"), sep = "-") -> df

makeGRangesFromDataFrame(df ) -> df_down
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
res_anno <- annotatePeak(df_down,TxDb = txdb, annoDb="org.Hs.eg.db")
as.GRanges(res_anno)
write.csv(as.data.frame(res_anno), file="ATAC_macro_down05005_peaks_annotated.csv")

#annotate peaks up
as.data.frame(rownames(res_up)) -> df
separate(df, 1, c("chr", "start"), sep = ":") -> df
separate(df, 2, c("start", "end"), sep = "-") -> df

makeGRangesFromDataFrame(df ) -> df_up
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
res_anno <- annotatePeak(df_up,TxDb = txdb, annoDb="org.Hs.eg.db")
as.GRanges(res_anno)
write.csv(as.data.frame(res_anno), file="ATAC_macro_up05005_peaks_annotated.csv")

#annotate peaks with adjp< 0.05
res_005 <- subset(res,res$padj<0.05 )
as.data.frame(rownames(res_005)) -> df
separate(df, 1, c("chr", "start"), sep = ":") -> df
separate(df, 2, c("start", "end"), sep = "-") -> df

makeGRangesFromDataFrame(df ) -> df_down
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
res_anno <- annotatePeak(df_down,TxDb = txdb, annoDb="org.Hs.eg.db")
as.GRanges(res_anno)
write.csv(as.data.frame(res_anno), file="ATAC_macro_adjp005_peaks_annotated.csv")

