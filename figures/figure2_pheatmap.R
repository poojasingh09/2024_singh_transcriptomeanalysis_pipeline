#singhpoo.singh09@gmail.com
#Singh et al
#plot heatmaps!

#udated september 2021
library(RColorBrewer)
library(tidyverse)
library(WGCNA)
library(pheatmap)
library(tidyr)
library(dplyr)
packageVersion("WGCNA")
packageVersion("pheatmap")
sessionInfo()
library(cowplot)
library(data.table)
library(viridis)



annotation_colors = list(lake=c(LT="goldenrod2", LM="#009E73", LV="#E49EC2", RV="#56B4E9"), jaw=c(OJ="black", PJ="grey"))


## colours: http://sape.inf.usi.ch/quick-reference/ggplot2/colour


########## read in input data ##########

anno <- read.table("/Users/singhpoo/Desktop/Papers/2019_Postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/phenotypes.txt", header=T)
anno1 <- data.frame(cbind(anno$lake, anno$jaw))
rownames(anno1) <- anno$ids
colnames(anno1) <- c("lake", "jaw")

anno1$jaw <- as.factor(anno1$jaw) ## super important for pheatmap or it wont read the colours
anno1$lake <- as.factor(anno1$lake)


g <- read.table("/Users/singhpoo/Desktop/Papers/2019_Postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/DESeq_NormalizedCounts_vsd.txt", header=T, row.names=1)
t <- read.table("/Users/singhpoo/Desktop/Papers/2019_Postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/DESeq_NormalizedCounts_vsd.txt", header=T, row.names=1)

g1 <- g[rowSums(g<5)<100, ]
t1 <- t[rowSums(t<5)<100, ]

#g1 <- sample_n(g1, 2000)
#t1 <- sample_n(t1, 2000)


##### FIGURE 2 heatmaps! ###################



#### gene

g1t <- t(g)


d <- dist(g1t, method = "euclidean")
#d <- dist(cor(g1), method = "euclidean")

fitg <- hclust(d, method="ward.D2")

col1 = c("black", "grey")
col2 = c("lightblue1", "brown", "cyan4", "pink")
col3 = c("red", "green", "blue")
col4 = c("darkgreen", "lightgreen", "orange", "white", "yellow")

pheno_data <- read.table("/Users/singhpoo/Desktop/Papers/2019_Postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/phenotypes.txt", header=T)

#diet2 <- pheno_data$diet2
diet1 <- pheno_data$diet1
lake <- pheno_data$lake
jaw <- pheno_data$jaw


#diet21 <- labels2colors(diet2,colorSeq=col4)
diet11 <- labels2colors(diet1, naColor="white",colorSeq=col3)
lake1 <- labels2colors(lake, naColor="white",colorSeq=col2)
jaw1 <- labels2colors(jaw, naColor="white",colorSeq=col1)

lab <- cbind(jaw1, lake1)
lab <- cbind(lab, diet11)
#lab <- cbind(lab, diet21)

colnames(lab) <- c("jaw","lake","diet1")
rownames(lab) <- pheno_data$ids

lab1 <- lab[match(rownames(g1t), rownames(lab)), ]
#lab2 <- lab[match(rownames(t1t), rownames(lab)), ]

svg("/Users/singhpoo/Desktop/Papers/2019_Postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/hclust_gene_expression.svg", h=6, w=10)
plotDendroAndColors(fitg, lab1, cex.dendroLabels = 0.4,colorHeight=0.1,autoColorHeight=F)
dev.off()


svg("Figure2A_pheatmap.svg")
p1 <- pheatmap(cor(g1), cex=0.9, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="ward.D2", color = viridis(12), annotation_colors=annotation_colors, annotation_row=anno1, fontsize_row = 8, fontsize_col=8)
p1
dev.off()




####trancript

t1t <- t(t1)

d <- dist(t1t, method = "euclidean")
#d <- dist(cor(t1), method = "euclidean")

fitt <- hclust(d, method="ward.D2")

col1 = c("black", "grey")
col2 = c("lightblue1", "brown", "cyan4", "pink")
col3 = c("red", "green", "blue")
col4 = c("darkgreen", "lightgreen", "orange", "white", "yellow")

pheno_data <- read.table("/Users/singhpoo/Desktop/Papers/2019_Postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/phenotypes.txt", header=T)

#diet2 <- pheno_data$diet2
diet1 <- pheno_data$diet1
lake <- pheno_data$lake
jaw <- pheno_data$jaw


#diet21 <- labels2colors(diet2,colorSeq=col4)
diet11 <- labels2colors(diet1, naColor="white",colorSeq=col3)
lake1 <- labels2colors(lake, naColor="white",colorSeq=col2)
jaw1 <- labels2colors(jaw, naColor="white",colorSeq=col1)

lab <- cbind(jaw1, lake1)
lab <- cbind(lab, diet11)
#lab <- cbind(lab, diet21)

colnames(lab) <- c("jaw","lake","diet1")
rownames(lab) <- pheno_data$ids

#lab1 <- lab[match(rownames(g1t), rownames(lab)), ]
lab2 <- lab[match(rownames(t1t), rownames(lab)), ]


pdf("/Users/singhpoo/Desktop/Papers/2019_Postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/hclust_transcript_expression.pdf", h=6, w=10)
plotDendroAndColors(fitt, lab2, cex.dendroLabels = 0.4,colorHeight=0.1,autoColorHeight=F)
dev.off()

svg("Figure2B_pheatmap.svg")
p2 <- pheatmap(cor(t1), cex=0.9, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="ward.D2", color = viridis(12), annotation_colors=annotation_colors, annotation_row=anno1, fontsize_row = 8, fontsize_col=8)
p2
dev.off()


## BOTH FIGURE FINAL

plot_list=list()
plot_list[['p1']]=p1[[4]]
plot_list[['p2']]=p2[[4]]

svg("Figure2_rowsum100_pheatmap.svg")
grid.arrange(grobs=plot_list, ncol=2)
dev.off()


