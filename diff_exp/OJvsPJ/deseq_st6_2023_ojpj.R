#Pooja Singh 
#AUg2020

library("DESeq2")

sessionInfo()
packageVersion("DESeq2")
library(EnhancedVolcano)
library(rtracklayer)
library(dplyr)
library("RColorBrewer")
library("gplots")


##########

setwd("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/")

#read in data
countData <- as.matrix(read.csv("../gene_count_matrix_prepDE.csv", row.names="gene_id"))

## after analysing PCAs, drop outlier samples"
drop <- c("Pp.L.OJ.3")
countData <- countData[ , !(colnames(countData) %in% drop)]


#######subselect 
colData <- read.csv("../phenotypes.txt", sep="\t", row.names=1)
colData <- colData[!rownames(colData) %in% drop, ]
colData$condition <- paste0(colData$lake, "_", colData$jaw)

#Check all sample IDs in colData are also in CountData and match their orders

all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))



dds1 <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ lib + condition)
dds1 <- dds1[ rowSums(counts(dds1)) > 10, ]
dds1 <- DESeq(dds1) 


png("pjoj_dispersion_plot.png", 1000, 1000)
plotDispEsts(dds1)
dev.off()



#matrix of normalized counts

vsd1 <- varianceStabilizingTransformation(dds1)
matrix1 <- (assay(vsd1))
write.table(matrix1, file = "DESeq_NormalizedCounts_vsd_pjoj.txt", sep = "\t", quote=F)

#########################################################################


#########################################################################LV

#make contrasts
resultsNames(dds1)

res_LV <- results(dds1, contrast=c("condition", "LV_PJ", "LV_OJ"), alpha=0.05)
res_LV_sig <- res_LV[which(res_LV$padj < 0.05),]
res_LV_sig <- res_LV_sig[order(res_LV_sig$padj),] 
write.table(res_LV, file = "LV_PJOJ_DESeq_fullresults.txt", sep = "\t", quote=F)
write.table(res_LV_sig, file = "LV_PJOJ_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)

##########################################################################
## plot VOLCANO plots with gene annotation

#################read anno so violin plots can be annotated
gff <- readGFF("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/O_niloticus_UMD_NMBU.99.gff3.genes")
gff$name_full <-  sapply(strsplit(gff$description, split=" [", fixed = TRUE), `[`, 1)
#write.table(gff, "O_niloticus_UMD_NMBU.99.gff3.genes.clean.txt", quote=F, row.names=F, sep="\t")
 
############ LV OJ VOLCANO
res_LV$gene_id <- sapply(strsplit(rownames(res_LV), split="|", fixed = TRUE), `[`, 2)
res_LV$our_id <- sapply(strsplit(rownames(res_LV), split="|", fixed = TRUE), `[`, 1)
res_LV_anno <- merge(data.frame(res_LV), gff, by="gene_id", all.x = TRUE)
res_LV_anno <- res_LV_anno %>%  mutate(name_full = coalesce(name_full,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LV_anno <- res_LV_anno %>%  mutate(Name = coalesce(Name,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LV_anno <- res_LV_anno[order(res_LV_anno$padj),] 
write.table(res_LV_anno, "LV_PJOJ_DESeq_fullresults.anno.txt", quote=F, row.names=F, sep="\t")
res_LV_anno_sig <- res_LV_anno[which(res_LV_anno$padj < 0.05),]
write.table(res_LV_anno_sig, "LV_PJOJ_DESeq_fullresults_0.05_anno.txt", quote=F, row.names=F, sep="\t")

lvtop <- head(res_LV_anno, n=5)

svg("LV_PJOJ_DESeq_fullresults.volcano_manuscript.svg")
EnhancedVolcano(res_LV_anno, lab = res_LV_anno$Name, x = 'log2FoldChange', y = 'pvalue', 
                pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, 
                gridlines.minor = FALSE, drawConnectors = TRUE, widthConnectors = 0.20, selectLab = lvtop$Name)
dev.off()




#########################################################################LM

#make contrasts
resultsNames(dds1)

res_LM <- results(dds1, contrast=c("condition", "LM_PJ", "LM_OJ"), alpha=0.05)
res_LM_sig <- res_LM[which(res_LM$padj < 0.05),]
res_LM_sig <- res_LM_sig[order(res_LM_sig$padj),] 
write.table(res_LM, file = "LM_PJOJ_DESeq_fullresults.txt", sep = "\t", quote=F)
write.table(res_LM_sig, file = "LM_PJOJ_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)

############ LM OJ VOLCANO
res_LM$gene_id <- sapply(strsplit(rownames(res_LM), split="|", fixed = TRUE), `[`, 2)
res_LM$our_id <- sapply(strsplit(rownames(res_LM), split="|", fixed = TRUE), `[`, 1)
res_LM_anno <- merge(data.frame(res_LM), gff, by="gene_id", all.x = TRUE)
res_LM_anno <- res_LM_anno %>%  mutate(name_full = coalesce(name_full,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LM_anno <- res_LM_anno %>%  mutate(Name = coalesce(Name,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LM_anno <- res_LM_anno[order(res_LM_anno$padj),] 
write.table(res_LM_anno, "LM_PJOJ_DESeq_fullresults.anno.txt", quote=F, row.names=F, sep="\t")
res_LM_anno_sig <- res_LM_anno[which(res_LM_anno$padj < 0.05),]
write.table(res_LM_anno_sig, "LM_PJOJ_DESeq_fullresults_0.05_anno.txt", quote=F, row.names=F, sep="\t")

lmtop <- head(res_LM_anno, n=10)

svg("LM_PJOJ_DESeq_fullresults.volcano_manuscript.svg")
EnhancedVolcano(res_LM_anno, lab = res_LM_anno$Name, x = 'log2FoldChange', y = 'pvalue', 
                pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, 
                gridlines.minor = FALSE, drawConnectors = TRUE, widthConnectors = 0.20, selectLab = lmtop$Name)
dev.off()



#########################################################################LT

#make contrasts
resultsNames(dds1)

res_LT <- results(dds1, contrast=c("condition", "LT_PJ", "LT_OJ"), alpha=0.05)
res_LT_sig <- res_LT[which(res_LT$padj < 0.05),]
res_LT_sig <- res_LT_sig[order(res_LT_sig$padj),] 
write.table(res_LT, file = "LT_PJOJ_DESeq_fullresults.txt", sep = "\t", quote=F)
write.table(res_LT_sig, file = "LT_PJOJ_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)

##########################################################################

############ LT OJ VOLCANO
res_LT$gene_id <- sapply(strsplit(rownames(res_LT), split="|", fixed = TRUE), `[`, 2)
res_LT$our_id <- sapply(strsplit(rownames(res_LT), split="|", fixed = TRUE), `[`, 1)
res_LT_anno <- merge(data.frame(res_LT), gff, by="gene_id", all.x = TRUE)
res_LT_anno <- res_LT_anno %>%  mutate(name_full = coalesce(name_full,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LT_anno <- res_LT_anno %>%  mutate(Name = coalesce(Name,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LT_anno <- res_LT_anno[order(res_LT_anno$padj),] 
write.table(res_LT_anno, "LT_PJOJ_DESeq_fullresults.anno.txt", quote=F, row.names=F, sep="\t")
res_LT_anno_sig <- res_LT_anno[which(res_LT_anno$padj < 0.05),]
write.table(res_LT_anno_sig, "LT_PJOJ_DESeq_fullresults_0.05_anno.txt", quote=F, row.names=F, sep="\t")

lttop <- head(res_LT_anno, n=5)

svg("LT_PJOJ_DESeq_fullresults.volcano_manuscript.svg")
EnhancedVolcano(res_LT_anno, lab = res_LT_anno$Name, x = 'log2FoldChange', y = 'pvalue', 
                pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, 
                gridlines.minor = FALSE, drawConnectors = TRUE, widthConnectors = 0.20, selectLab = lttop$Name)
dev.off()

