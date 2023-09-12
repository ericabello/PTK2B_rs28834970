library("RColorBrewer")
library("pheatmap")
library("DESeq2")
library("ggplot2")
library("genefilter")
library(tidyverse)
library(BiocParallel)
#register(MulticoreParam(4))
#options(stringsAsFactors = F)
library("org.Hs.eg.db")
library(dplyr)
library("ggrepel")
options(ggrepel.max.overlaps = Inf) 
library(gprofiler2)
library(stringr)
library(data.table)
library(gridExtra)
library(gplots)
library(gridExtra)

setwd("/Users/eb19/Documents/data_PTK2B_microglia/RNAseq/final files feb 2021/")
#import counts table from FeatureCounts output
read.table("../counts/counts_ptk2b_micro_271119.txt" , header = T, stringsAsFactors = F, sep = "\t") -> fc_table
head(fc_table)



#fix col names samples and divide into 2 df, one with gene meta data and one with counts per sample
colnames(fc_table)[c(7:24)] <- c("WT_B11_1","WT_B11_2","WT_B11_3","KOLF_1","KOLF_2","KOLF_3", "HOM_C11_1","HOM_C11_2","HOM_C11_3", "HOM_D3_1","HOM_D3_2","HOM_D3_3","HOM_C8_1","HOM_C8_2","HOM_C8_3")
remove_dot_Geneid <- function(x) {x$Geneid <- gsub("\\.[\\d]+", "", x$Geneid, perl=T)
return(x) }
fc_table <- remove_dot_Geneid(fc_table)

gene.meta = fc_table[,c(1,6)] 
counts_macro = fc_table[,-c(1:6)]
rownames(counts_macro) <- fc_table[,1]
head(counts_macro)


#make legend df with factors
legend <- data.frame(row.names = colnames(counts_macro), genotype_PTK2B=as.factor(c(rep("WT", 6), rep("HOM",9))))
legend


#calculate tpms from raw counts FC. a is the counts df
tpms <- function (a) { rpk = apply(a, MARGIN=2, FUN=function(x) x / gene.meta$Length) * 1e3
tpm = apply(rpk, MARGIN=2, FUN=function(x) x / sum(x)) * 1e6
tpm = tpm %>% as.data.frame()
return(tpm)}
tpm <- tpms(counts_macro)

#filter TPMs. x is tpm df and n is mean filter (minimum average TPM in all samples to consider gene as expressed)
tpm_filter <- function(x, n) { expr <- tpm %>% filter(rowMeans(.) > n)
return (expr)}
tpm_expressed <- tpm_filter(tpm, 1)
head(tpm_expressed)
#add symnbols to tpm df and save them
annots <- select(org.Hs.eg.db, keys=rownames(tpm_expressed),
                 columns=c("SYMBOL"), keytype="ENSEMBL")

annots = annots[!duplicated(annots["ENSEMBL"]),]
tpm_expressed$ENSEMBL <- rownames(tpm_expressed)
tpm_expressed <- inner_join(annots, tpm_expressed, by= "ENSEMBL")
tpm_expressed$ENSEMBL -> rownames(tpm_expressed)
#write.csv(tpm_expressed, 
#file="/Users/eb19/Documents/data_PTK2B_microglia/RNAseq/final files feb 2021/PTK2B_microglia_TPMs_noF4.csv")


#do DE for genes with mean at least 1 or TPM (filter raw counts by tpm_expressed)
counts_macro <- counts_macro %>% filter(rownames(counts_macro) %in% rownames(tpm_expressed))

#make deseq2 df and DE analysis 
coldata <- legend 
coldata
dds <- DESeqDataSetFromMatrix(countData = counts_macro,
                              colData = coldata,
                              design = ~genotype_PTK2B)
dds <- estimateSizeFactors(dds)
dds$genotype_PTK2B <- relevel(dds$genotype_PTK2B, ref = "WT")

dds <- DESeq(dds)

counts_normalised <- counts(dds, normalized=TRUE) %>% as.data.frame()
#write.csv(counts_normalised, 
#file="/Users/eb19/Documents/data_PTK2B_microglia/RNAseq/final files feb 2021/PTK2B_microglia_normalisedCounts_DESeq2_noF4.csv")

#results df (p vales)
res <- results(dds)
summary(res)


#customised PCA to save
#for custom pca plot return df of PC1 and 2
rlog <- rlog(dds, blind=TRUE)
p <- rlog %>% plotPCA("genotype_PTK2B")
p <- p + geom_text(aes_string(label = "name"), color = "black") + theme(text = element_text(size = 14))+theme_classic()+theme(aspect.ratio = 1)
plot(p)
#for custom pca plot return df of PC1 and 2
df <- plotPCA(rlog, "genotype_PTK2B", returnData = TRUE)
a <- c("1","2","3","4", "5")
df$hiPSC_clone <- c(rep(a, each=3))
#custom PCA plot
#ggplot(data= df, mapping = aes(x=PC1 , y=PC2))+geom_point(size=4, mapping = aes(color=genotype_PTK2B, shape=hiPSC_clone))+theme_classic()+theme(aspect.ratio = 1)
#to change shapes of each sample manually
ggplot(data= df, mapping = aes(x=PC1 , y=PC2))+geom_point(size=4, mapping = aes(color=genotype_PTK2B, shape=hiPSC_clone))+ scale_color_manual(values=c("orange", "blue"))+scale_shape_manual(values=c(17, 18, 16,19, 15,20))+theme_classic()+theme(aspect.ratio = 1,legend.title = element_text(color = "black", size = 8),legend.text = element_text(color = "black", size = 6), legend.key.size = unit(0.3, "cm"),legend.margin=margin(0,0,0,0),
                                                                                                                                                                                                                                                   legend.box.margin=margin(0,0,0,0), legend.box.spacing = margin(10,10,10,10))+ xlab("PC1 66% variance") + ylab("PC2 27% variance")
ggsave("/Users/eb19/Documents/data_PTK2B_microglia/RNAseq/final files feb 2021/PCA_rog_micro_noF4_noNames.jpeg", device= "jpeg")


#annotate genes in res and make resultTable df with symbols and save it
resultTable <- res %>% as.data.frame()
annots <- select(org.Hs.eg.db, keys=rownames(resultTable),
                 columns=c("SYMBOL","GENENAME"), keytype="ENSEMBL")

annots = annots[!duplicated(annots["ENSEMBL"]),]
resultTable$ENSEMBL <- rownames(resultTable)
resultTable <- inner_join(annots, resultTable, by= "ENSEMBL")
write.csv(resultTable, 
          file="/Users/eb19/Documents/data_PTK2B_microglia/RNAseq/final files feb 2021/PTK2B_microglia_resultsDESEq2Annotated_noF4.csv")


#vulcano plot with colours logfc 0.5/-0.5 
#resultTable$diffExpr <- "NO"
#resultTable$diffExpr[resultTable$log2FoldChange > 1 & resultTable$padj < 0.05] <- "UP"
#resultTable$diffExpr[resultTable$log2FoldChange < -(1) & resultTable$padj < 0.05] <- "DOWN"
#resultTable$delabel <- NA
#resultTable$delabel[resultTable$diffExpr != "NO"] <- resultTable$SYMBOL[resultTable$diffExpr != "NO"]
res_plot <- resultTable %>% mutate(diffExpr = case_when(padj < 0.05 & log2FoldChange > 0.5 ~ "UP", padj < 0.05 & log2FoldChange < -0.5 ~ "DOWN")) %>% replace_na(list(diffExpr="NO")) 
res_plot <- res_plot %>% mutate(delabel = case_when(res_plot$diffExpr != "NO"~ res_plot$SYMBOL))
ggplot(data=res_plot, aes(x=log2FoldChange, y=-log10(padj), color= diffExpr)) + geom_point(size=1,aes(col=diffExpr)) + xlim(-10,10) + scale_color_manual(values=c("blue", "black", "red")) +theme_minimal()
ggsave('/Users/eb19/Documents/data_PTK2B_microglia/RNAseq/final files feb 2021/vulcanoPlot_noF4_micro.jpeg',device= "jpeg", width = 8)

#filters results by pvalue and logFC after dropping rows of genes with p Valuse NA. saves csv of genes up and down with logFC and padj, saves csv with summary how many up and down, makes 2 df (one up and one down)
filter_DE_genes <- function (x, n, p) {
  df <- x %>% drop_na(padj) %>% filter((log2FoldChange > n | log2FoldChange < -n) & padj < p)  %>% mutate(diffExpr = if_else(log2FoldChange >= n, 'UP', 'DOWN'))
  fileName=paste0("/Users/eb19/Documents/data_PTK2B_microglia/RNAseq/final files feb 2021/PTK2B_micro_genes_up_down_p005_FC",n,".csv")
  write.csv(df, file=fileName)
  count <- df %>% count(diffExpr)
  fileName=paste0("/Users/eb19/Documents/data_PTK2B_microglia/RNAseq/final files feb 2021/PTK2B_res_micro_summary_p005_FC",n,".csv")
  write.csv(count, file=fileName)
  up <- df %>% filter(log2FoldChange > n ) %>% dplyr::select(.,ENSEMBL,log2FoldChange,padj)
  down <- df %>% filter(log2FoldChange < n ) %>% dplyr::select(.,ENSEMBL,log2FoldChange,padj)
  output <- list(up=up, down=down)
  return(output)
}
#apply filter function for LogFc 1
logFC1 <- filter_DE_genes(resultTable, 1, 0.05)
#extract data and sort them so ordered for GSEA
up_1 <- data.frame(logFC1$up) %>% arrange(desc(log2FoldChange))
down_1 <- as.data.frame(logFC1$down) %>% arrange(log2FoldChange)

#apply filter function for LogFc 0.5
logFC05 <- filter_DE_genes(resultTable, 0.5, 0.05)
#extract data and sort them so ordered for GSEA
up_0.5 <- as.data.frame(logFC05$up) %>% arrange(desc(log2FoldChange))
down_0.5 <- as.data.frame(logFC05$down) %>% arrange(log2FoldChange)

#GSEA of DE genes
library(gprofiler2)
#library(clusterProfiler)
#extract list of expressed genes as background
#from filtered results df extract list of ordered gene symbols -TRY BOTH ENSEMBLE AND SYMBOLS AND SEE IF I GET THE SAME WITH NAs(pseudogenes)
genes_expressed <- resultTable$ENSEMBL
genes_up_1 <- up_1$ENSEMBL
genes_down_1 <- down_1$ENSEMBL
genes_up_05 <- up_0.5$ENSEMBL
genes_down_05 <- down_0.5$ENSEMBL
#run GProfiler: re run with "as_short_link = TRUE" to get web link only to share
#function to run GProfiler with set paramenters on different gene sets
###evcodes= TRUE (before (try false to see if it finds more genes))not modified to FALSE to include more genes, genes in GO are only more stringent ones compared to GProfiler run on web server)
run_GProfiler <- function (x,y) {df <- gost(query = x, 
                                            organism = "hsapiens", ordered_query = TRUE, 
                                            significant = TRUE, exclude_iea = TRUE, 
                                            evcodes = TRUE, 
                                            user_threshold = 0.05, correction_method = "bonferroni", 
                                            domain_scope = "custom_annotated", custom_bg = y, 
                                            as_short_link = FALSE) 
return (df)}
#run GProfiler to all genes DE logfc 0.5 or 1 (and extract result, $meta will give you metadata of analysis)
GProfiler_up_1 <- run_GProfiler(genes_up_1, genes_expressed)$result
GProfiler_down_1 <- run_GProfiler(genes_down_1, genes_expressed)$result

GProfiler_up_05 <- run_GProfiler(genes_up_05, genes_expressed)$result
GProfiler_down_05 <- run_GProfiler(genes_down_05, genes_expressed)$result


#save GProfiler results for all genes DE
save_csv_GProf <- function(x) {df_name = deparse(substitute(x))
fileName=paste0("/Users/eb19/Documents/data_PTK2B_microglia/RNAseq/final files feb 2021/PTK2B_microglia_noF4_",df_name,".csv")
fwrite(x, file=fileName)}
save_csv_GProf(GProfiler_down_05)
save_csv_GProf(GProfiler_down_1)
save_csv_GProf(GProfiler_up_05)
save_csv_GProf(GProfiler_up_1)

