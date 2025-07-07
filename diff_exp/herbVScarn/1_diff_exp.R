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

setwd("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023")

load("deseq2.RData")

#read in data
countData <- as.matrix(read.csv("gene_count_matrix_prepDE.csv", row.names="gene_id"))

## after analysing PCAs, drop outlier samples"
drop <- c("Pp.L.OJ.3")
countData <- countData[ , !(colnames(countData) %in% drop)]


##subset

countData1 <- countData[,grepl("OJ", colnames(countData))]
countData2 <- countData[,grepl("PJ", colnames(countData))]


#######subselect 
colData <- read.csv("phenotypes.txt", sep="\t", row.names=1)

colData <- colData[!rownames(colData) %in% drop, ]
colData$condition <- paste0(colData$lake, "_", colData$diet1)

colData1 <- subset(colData, colData$jaw == "OJ")
colData2 <- subset(colData, colData$jaw == "PJ")


#Check all sample IDs in colData are also in CountData and match their orders

all(rownames(colData1) %in% colnames(countData1))
countData1 <- countData1[, rownames(colData1)]
all(rownames(colData1) == colnames(countData1))


all(rownames(colData2) %in% colnames(countData2))
countData2 <- countData2[, rownames(colData2)]
all(rownames(colData2) == colnames(countData2))

#oj
dds1 <- DESeqDataSetFromMatrix(countData = countData1, colData = colData1, design = ~ lib + condition)
dds1 <- dds1[ rowSums(counts(dds1)) > 10, ]
dds1 <- DESeq(dds1) 

#pj

dds2 <- DESeqDataSetFromMatrix(countData = countData2, colData = colData2, design = ~ lib + condition)
dds2 <- dds2[ rowSums(counts(dds2)) > 10, ]
dds2 <- DESeq(dds2) 


#png("dispersion_plot_OJ.png", 1000, 1000)
plotDispEsts(dds1)
#dev.off()


png("dispersion_plot_PJ.png", 1000, 1000)
plotDispEsts(dds2)
dev.off()


#matrix of normalized counts

vsd1 <- varianceStabilizingTransformation(dds1)
matrix1 <- (assay(vsd1))
#write.table(matrix1, file = "DESeq_NormalizedCounts_vsd_OJ.txt", sep = "\t", quote=F)

vsd2 <- varianceStabilizingTransformation(dds2)
matrix2 <- (assay(vsd2))
#write.table(matrix2, file = "DESeq_NormalizedCounts_vsd_PJ.txt", sep = "\t", quote=F)

#########################################################################


#########################################################################

#make contrasts
resultsNames(dds1)
resultsNames(dds2)

## OJ
res_LV <- results(dds1, contrast=c("condition", "LV_herbivore", "LV_carnivore"), alpha=0.05)
res_LV_sig <- res_LV[which(res_LV$padj < 0.05),]
res_LV_sig <- res_LV_sig[order(res_LV_sig$padj),] 
#write.table(res_LV, file = "LV_OJ_herbcarn_DESeq_fullresults.txt", sep = "\t", quote=F)
#write.table(res_LV_sig, file = "LV_OJ_herbcarn_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)

res_LM <- results(dds1, contrast=c("condition", "LM_herbivore", "LM_carnivore"), alpha=0.05)
res_LM_sig <- res_LM[which(res_LM$padj < 0.05),]
res_LM_sig <- res_LM_sig[order(res_LM_sig$padj),] 
#write.table(res_LM, file = "LM_OJ_herbcarn_DESeq_fullresults.txt", sep = "\t", quote=F)
#write.table(res_LM_sig, file = "LM_OJ_herbcarn_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)


res_LT <- results(dds1, contrast=c("condition", "LT_herbivore", "LT_carnivore"), alpha=0.05)
res_LT_sig <- res_LT[which(res_LT$padj < 0.05),]
res_LT_sig <- res_LT_sig[order(res_LT_sig$padj),] 
#write.table(res_LT, file = "LT_OJ_herbcarn_DESeq_fullresults.txt", sep = "\t", quote=F)
#write.table(res_LT_sig, file = "LT_OJ_herbcarn_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)


##########################################################################
## plot VOLCANO plots with gene annotation

#################read anno so violin plots can be annotated
gff <- readGFF("O_niloticus_UMD_NMBU.99.gff3.genes")
gff$name_full <-  sapply(strsplit(gff$description, split=" [", fixed = TRUE), `[`, 1)
#write.table(gff, "O_niloticus_UMD_NMBU.99.gff3.genes.clean.txt", quote=F, row.names=F, sep="\t")
 
############ LV OJ VOLCANO
res_LV$gene_id <- sapply(strsplit(rownames(res_LV), split="|", fixed = TRUE), `[`, 2)
res_LV$our_id <- sapply(strsplit(rownames(res_LV), split="|", fixed = TRUE), `[`, 1)
res_LV_anno <- merge(data.frame(res_LV), gff, by="gene_id", all.x = TRUE)
res_LV_anno <- res_LV_anno %>%  mutate(name_full = coalesce(name_full,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LV_anno <- res_LV_anno %>%  mutate(Name = coalesce(Name,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LV_anno <- res_LV_anno[order(res_LV_anno$padj),] 
write.table(res_LV_anno, "LV_OJ_herbcarn_DESeq_fullresults.anno.txt", quote=F, row.names=F, sep="\t")
res_LV_anno_sig <- res_LV_anno[which(res_LV_anno$padj < 0.05),]
write.table(res_LV_anno_sig, "LV_OJ_herbcarn_DESeq_fullresults_0.05_anno.txt", quote=F, row.names=F, sep="\t")


svg("LV_OJ_herbcarn_DESeq_fullresults.volcano.svg", width = 9, height = 9)
EnhancedVolcano(res_LV_anno, lab = res_LV_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()

pdf("LV_OJ_herbcarn_DESeq_fullresults.volcano.pdf")
EnhancedVolcano(res_LV_anno, lab = res_LV_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()


############ LM OJ VOLCANO
res_LM$gene_id <- sapply(strsplit(rownames(res_LM), split="|", fixed = TRUE), `[`, 2)
res_LM$our_id <- sapply(strsplit(rownames(res_LM), split="|", fixed = TRUE), `[`, 1)
res_LM_anno <- merge(data.frame(res_LM), gff, by="gene_id", all.x = TRUE)
res_LM_anno <- res_LM_anno %>%  mutate(name_full = coalesce(name_full,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LM_anno <- res_LM_anno %>%  mutate(Name = coalesce(Name,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LM_anno <- res_LM_anno[order(res_LM_anno$padj),] 
write.table(res_LM_anno, "LM_OJ_herbcarn_DESeq_fullresults.anno.txt", quote=F, row.names=F, sep="\t")
res_LM_anno_sig <- res_LM_anno[which(res_LM_anno$padj < 0.05),]
write.table(res_LM_anno_sig, "LM_OJ_herbcarn_DESeq_fullresults_0.05_anno.txt", quote=F, row.names=F, sep="\t")

svg("LM_OJ_herbcarn_DESeq_fullresults.volcano.svg", width = 9, height = 9)
EnhancedVolcano(res_LM_anno, lab = res_LM_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()

pdf("LM_OJ_herbcarn_DESeq_fullresults.volcano.pdf")
EnhancedVolcano(res_LM_anno, lab = res_LM_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()


############ LT OJ VOLCANO
res_LT$gene_id <- sapply(strsplit(rownames(res_LT), split="|", fixed = TRUE), `[`, 2)
res_LT$our_id <- sapply(strsplit(rownames(res_LT), split="|", fixed = TRUE), `[`, 1)
res_LT_anno <- merge(data.frame(res_LT), gff, by="gene_id", all.x = TRUE)
res_LT_anno <- res_LT_anno %>%  mutate(name_full = coalesce(name_full,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LT_anno <- res_LT_anno %>%  mutate(Name = coalesce(Name,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LT_anno <- res_LT_anno[order(res_LT_anno$padj),] 
write.table(res_LT_anno, "LT_OJ_herbcarn_DESeq_fullresults.anno.txt", quote=F, row.names=F, sep="\t")
res_LT_anno_sig <- res_LT_anno[which(res_LT_anno$padj < 0.05),]
write.table(res_LT_anno_sig, "LT_OJ_herbcarn_DESeq_fullresults_0.05_anno.txt", quote=F, row.names=F, sep="\t")

#res_LT_anno <- read.table("LT_OJ_herbcarn_DESeq_fullresults.anno.txt", header=T, sep="\t", fill=TRUE)
res_LT_anno <- read.table("LT_OJ_herbcarn_DESeq_fullresults.anno.txt", header=T, sep="\t")
svg("LT_OJ_herbcarn_DESeq_fullresults.anno_withoutENSONIG.svg", width = 9, height = 9)
EnhancedVolcano(res_LT_anno, lab = ifelse(grepl("^ENSONIG", res_LT_anno$Name), NA, res_LT_anno$Name), x = 'log2FoldChange', 
                y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, 
                gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()

pdf("LT_OJ_herbcarn_DESeq_fullresults.anno.pdf")
EnhancedVolcano(res_LT_anno, lab = res_LT_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()

##########################################################################
##PJ
res_LV <- results(dds2, contrast=c("condition", "LV_herbivore", "LV_carnivore"), alpha=0.05)
res_LV_sig <- res_LV[which(res_LV$padj < 0.05),]
res_LV_sig <- res_LV_sig[order(res_LV_sig$padj),] 
write.table(res_LV, file = "LV_PJ_herbcarn_DESeq_fullresults.txt", sep = "\t", quote=F)
write.table(res_LV_sig, file = "LV_PJ_herbcarn_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)


res_LM <- results(dds2, contrast=c("condition", "LM_herbivore", "LM_carnivore"), alpha=0.05)
res_LM_sig <- res_LM[which(res_LM$padj < 0.05),]
res_LM_sig <- res_LM_sig[order(res_LM_sig$padj),] 
write.table(res_LM, file = "LM_PJ_herbcarn_DESeq_fullresults.txt", sep = "\t", quote=F)
write.table(res_LM_sig, file = "LM_PJ_herbcarn_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)


res_LT <- results(dds2, contrast=c("condition", "LT_herbivore", "LT_carnivore"), alpha=0.05)
res_LT_sig <- res_LT[which(res_LT$padj < 0.05),]
res_LT_sig <- res_LT_sig[order(res_LT_sig$padj),] 
write.table(res_LT, file = "LT_PJ_herbcarn_DESeq_fullresults.txt", sep = "\t", quote=F)
write.table(res_LT_sig, file = "LT_PJ_herbcarn_DESeq_fullresults_0.05.txt", sep = "\t", quote=F)

##########################################################################
## plot VOLCANO plots with gene annotation


############ LV PJ VOLCANO
res_LV$gene_id <- sapply(strsplit(rownames(res_LV), split="|", fixed = TRUE), `[`, 2)
res_LV$our_id <- sapply(strsplit(rownames(res_LV), split="|", fixed = TRUE), `[`, 1)
res_LV_anno <- merge(data.frame(res_LV), gff, by="gene_id", all.x = TRUE)
res_LV_anno <- res_LV_anno %>%  mutate(name_full = coalesce(name_full,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LV_anno <- res_LV_anno %>%  mutate(Name = coalesce(Name,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LV_anno <- res_LV_anno[order(res_LV_anno$padj),] 
write.table(res_LV_anno, "LV_PJ_herbcarn_DESeq_fullresults.anno.txt", quote=F, row.names=F, sep="\t")
res_LV_anno_sig <- res_LV_anno[which(res_LV_anno$padj < 0.05),]
write.table(res_LV_anno_sig, "LV_PJ_herbcarn_DESeq_fullresults_0.05_anno.txt", quote=F, row.names=F, sep="\t")

svg("LV_PJ_herbcarn_DESeq_fullresults.anno.svg", width = 9, height = 9)
EnhancedVolcano(res_LV_anno, lab = res_LV_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()

pdf("LV_PJ_herbcarn_DESeq_fullresults.anno.pdf")
EnhancedVolcano(res_LV_anno, lab = res_LV_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()



############ LM PJ VOLCANO
res_LM$gene_id <- sapply(strsplit(rownames(res_LM), split="|", fixed = TRUE), `[`, 2)
res_LM$our_id <- sapply(strsplit(rownames(res_LM), split="|", fixed = TRUE), `[`, 1)
res_LM_anno <- merge(data.frame(res_LM), gff, by="gene_id", all.x = TRUE)
res_LM_anno <- res_LM_anno %>%  mutate(name_full = coalesce(name_full,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LM_anno <- res_LM_anno %>%  mutate(Name = coalesce(Name,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LM_anno <- res_LM_anno[order(res_LM_anno$padj),] 
write.table(res_LM_anno, "LM_PJ_herbcarn_DESeq_fullresults.anno.txt", quote=F, row.names=F, sep="\t")
res_LM_anno_sig <- res_LM_anno[which(res_LM_anno$padj < 0.05),]
write.table(res_LM_anno_sig, "LM_PJ_herbcarn_DESeq_fullresults_0.05_anno.txt", quote=F, row.names=F, sep="\t")

svg("LM_PJ_herbcarn_DESeq_fullresults.anno.svg", width = 9, height = 9)
EnhancedVolcano(res_LM_anno, lab = res_LM_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()

pdf("LM_PJ_herbcarn_DESeq_fullresults.anno.pdf")
EnhancedVolcano(res_LM_anno, lab = res_LM_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()


############ LT PJ VOLCANO
res_LT$gene_id <- sapply(strsplit(rownames(res_LT), split="|", fixed = TRUE), `[`, 2)
res_LT$our_id <- sapply(strsplit(rownames(res_LT), split="|", fixed = TRUE), `[`, 1)
res_LT_anno <- merge(data.frame(res_LT), gff, by="gene_id", all.x = TRUE)
res_LT_anno <- res_LT_anno %>%  mutate(name_full = coalesce(name_full,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LT_anno <- res_LT_anno %>%  mutate(Name = coalesce(Name,gene_id)) ## any missing gene names will be replaced with ENSEMBLE IDS
res_LT_anno <- res_LT_anno[order(res_LT_anno$padj),] 
write.table(res_LT_anno, "LT_PJ_herbcarn_DESeq_fullresults.anno.txt", quote=F, row.names=F, sep="\t")
res_LT_anno_sig <- res_LT_anno[which(res_LT_anno$padj < 0.05),]
write.table(res_LT_anno_sig, "LT_PJ_herbcarn_DESeq_fullresults_0.05_anno.txt", quote=F, row.names=F, sep="\t")

svg("LT_PJ_herbcarn_DESeq_fullresults.anno.svg", width = 9, height = 9)
EnhancedVolcano(res_LT_anno, lab = res_LT_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()

pdf("LT_PJ_herbcarn_DESeq_fullresults.anno.pdf")
EnhancedVolcano(res_LT_anno, lab = res_LT_anno$Name, x = 'log2FoldChange', y = 'pvalue', pCutoff=0.01, FCcutoff=2, labSize=4, pointSize=1, gridlines.major = FALSE, gridlines.minor = FALSE, col = c("grey50", "grey50", "grey50", "red2"))
dev.off()


save.image(file="deseq2.RData")

# 
# 
# ##########################################################################
# 
# #plot MA plot
# 
# png("LV_DEseq_MAplot_0.05.png")
# plotMA(res_LV, main="DESeq2_0.05", ylim=c(-2,2))
# dev.off()
# 
# png("LM_DEseq_MAplot_0.05.png")
# plotMA(res_LM, main="DESeq2_0.05", ylim=c(-2,2))
# dev.off()
# 
# png("LT_DEseq_MAplot_0.05.png")
# plotMA(res_LT, main="DESeq2_0.05", ylim=c(-2,2))
# dev.off()
# 
# ##########################################################################
# 
# #heatmaps of sample to samples distances: clustering
# 
# png("DEseq_clustering_heatmap_ofsamples_rld.png", 1000, 1000)
# sampleDists <- dist(t(assay(rld)))
# library("RColorBrewer") 
# sampleDistMatrix <- as.matrix(sampleDists) 
# rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-") 
# colnames(sampleDistMatrix) <- NULL 
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
# pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, col=colors)
# dev.off()
# 
# 
# png("DEseq_clustering_heatmap_ofsamples_vsd.png", 1000, 1000)
# sampleDists <- dist(t(assay(vsd)))
# library("RColorBrewer") 
# sampleDistMatrix <- as.matrix(sampleDists) 
# rownames(sampleDistMatrix) <- paste(vsd$lake, vsd$type, sep="-") 
# colnames(sampleDistMatrix) <- NULL 
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
# pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, col=colors, cex=0.8)
# dev.off()
# 
# #######################################################################
# 
# 
# #principal component plot of samples
# 
# png("PCA_samples_vsd_lib.png", width = 10, height = 8, units = 'in', res = 600, pointsize=8)
# plotPCA(vsd, intgroup=c("paired"))
# dev.off()
# 
# #PCA sample name
# 
# png("PCA_samples_vsd_outliers_identify.png", width = 10, height = 8, units = 'in', res = 600, pointsize=8)
# z <- plotPCA(vsd, intgroup=c("lib"))
# z + geom_label(aes(label = name))
# dev.off()
# 
# 
# data <- plotPCA(vsd, intgroup=c("lib"), returnData=TRUE) 
# percentVar <- round(100 * attr(data, "percentVar")) 
# ggplot(data, aes(PC1, PC2, color=lib)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
# 
# 
# data <- plotPCA(vsd, intgroup=c("jaw", "lake"), returnData=TRUE) 
# percentVar <- round(100 * attr(data, "percentVar")) 
# ggplot(data, aes(PC1, PC2, color=lake, shape=jaw)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
# 
# 
# 
# #OJ PJA
# 
# data <- plotPCA(vsd1, intgroup=c("diet1", "species"), returnData=TRUE) 
# percentVar <- round(100 * attr(data, "percentVar")) 
# ggplot(data, aes(PC1, PC2, color=species, shape=diet1)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
# 
# data <- plotPCA(vsd1, intgroup=c("lake", "diet1"), returnData=TRUE) 
# percentVar <- round(100 * attr(data, "percentVar")) 
# ggplot(data, aes(PC1, PC2, color=diet1, shape=lake)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
# 
# # PJ PCA
# 
# 
# data <- plotPCA(vsd2, intgroup=c("lake", "species"), returnData=TRUE) 
# percentVar <- round(100 * attr(data, "percentVar")) 
# ggplot(data, aes(PC1, PC2, color=species, shape=lake)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance"))
# 
# data <- plotPCA(vsd2, intgroup=c("species", "diet1"), returnData=TRUE) 
# percentVar <- round(100 * attr(data, "percentVar")) 
# ggplot(data, aes(PC1, PC2, color=species, shape=diet1)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")
