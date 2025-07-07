#pooja.singh09@gmail.com
#Singh et al
#plot heatmaps!

#udated september 2024
library(RColorBrewer)
library(tidyverse)
#library(WGCNA)
library(pheatmap)
library(tidyr)
library(dplyr)
packageVersion("WGCNA")
packageVersion("pheatmap")
sessionInfo()
library(cowplot)
library(data.table)
library(viridis)
library(gridExtra)
#install.packages("wesanderson")
library(wesanderson)


## figure 2
setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/manuscript')

#annotation_colors = list(lake=c(LT="goldenrod2", LM="#009E73", LV="#E49EC2", RV="#56B4E9"), jaw=c(OJ="black", PJ="grey"))
annotation_colors = list(lake=c(LT="#FAD510", LM="#00A08A", LV="#FD6467", RV="#56B4E9"), jaw=c(OJ="black", PJ="grey"))



## colours: http://sape.inf.usi.ch/quick-reference/ggplot2/colour


########## read in input data ##########

anno <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/phenotypes.txt", header=T)
anno1 <- data.frame(cbind(anno$lake, anno$jaw))
rownames(anno1) <- anno$ids
colnames(anno1) <- c("lake", "jaw")

anno1$jaw <- as.factor(anno1$jaw) ## super important for pheatmap or it wont read the colours
anno1$lake <- as.factor(anno1$lake)


anno2 <- data.frame(cbind(anno$lake, anno$jaw, anno$species))
rownames(anno2) <- anno$ids
colnames(anno2) <- c("lake", "jaw", "species")

anno1$jaw <- as.factor(anno1$jaw) ## super important for pheatmap or it wont read the colours
anno1$lake <- as.factor(anno1$lake)

anno2$jaw <- as.factor(anno2$jaw) ## super important for pheatmap or it wont read the colours
anno2$lake <- as.factor(anno2$lake)
anno2$species <- as.factor(anno2$species)


g <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/DESeq_NormalizedCounts_vsd.txt", header=T, row.names=1)

##### FIGURE 2 heatmaps! ###################



#### gene


gcorr <-cor(g,method="spearman")


#svg("Figure2A_pheatmap.svg")
p1 <- pheatmap(gcorr, cex=0.9, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="ward.D2", color = viridis(256), annotation_colors=annotation_colors, annotation_row=anno1, fontsize_row = 8, fontsize_col=8)
p1
#dev.off()

p1 <- pheatmap(gcorr, cex=0.9, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="ward.D2", color = viridis(256), annotation_colors=annotation_colors, annotation_row=anno1, fontsize_row = 8, fontsize_col=8)
p1



#svg("Figure2A_pheatmap_supp.svg")
p1supp <- pheatmap(gcorr, cex=0.9, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="ward.D2", color = viridis(256), annotation_row=anno2, fontsize_row = 8, fontsize_col=8, annotation_colors=annotation_colors)
p1supp
#dev.off()

####PSI FROM SUPPA

a <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/suppa/filtered.gffcmp.all.annotated_isoform.psi", header=T)

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
a[is.nan.data.frame(a)] <- 0

acorr <-cor(a,method="spearman")

#svg("Figure2B_pheatmap_psi.svg")
p2 <- pheatmap(acorr, cex=0.9, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="ward.D2", color = viridis(256), annotation_colors=annotation_colors, annotation_row=anno1, fontsize_row = 8, fontsize_col=8)
p2
#dev.off()


#svg("Figure2B_pheatmap_psi_supp.svg")
p2supp <- pheatmap(acorr, cex=0.9, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="ward.D2", color = viridis(256), annotation_row=anno2, fontsize_row = 8, fontsize_col=8, annotation_colors=annotation_colors)
p2supp
#dev.off()


## PSI  rho of genes with more than one isoform only

gene_names <- sapply(strsplit(rownames(a), ";"), `[`, 1)
gene_counts <- table(gene_names)
genes_1iso <- names(gene_counts[gene_counts == 1])
print(genes_1iso)
a_filtered <- a[!(gene_names %in% genes_1iso), ]
dim(a_filtered)
dim(a)


is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
a_filtered[is.nan.data.frame(a_filtered)] <- 0

a_filteredcorr <-cor(a_filtered,method="spearman")

svg("Figure2B_pheatmap_psi_onlysplicedgenes.svg")
p2f <- pheatmap(a_filteredcorr, cex=0.9, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="ward.D2", color = viridis(256), annotation_colors=annotation_colors, annotation_row=anno1, fontsize_row = 8, fontsize_col=8)
p2f
dev.off()


## BOTH FIGURE FINAL

plot_list=list()
plot_list[['p1']]=p1[[4]]
plot_list[['p2']]=p2[[4]]

#svg("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/manuscript/Figure2_pheatmap_psi.svg")
grid.arrange(grobs=plot_list, ncol=1)
#dev.off()

## BOTH FIGURE FINAL 2025

plot_list=list()
plot_list[['p1']]=p1[[4]]
plot_list[['p2f']]=p2f[[4]]

svg("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/manuscript/Figure2_pheatmap_psi_onlygeneswithmorethan1iso_supp.svg", h=10, w=8)
grid.arrange(grobs=plot_list, ncol=2)
dev.off()


## supp mat wiht species labels

plot_list=list()
plot_list[['p1']]=p1[[4]]
plot_list[['p2f']]=p2f[[4]]

svg("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/manuscript/Figure2_pheatmap_psi_supp_withlabels.svg", h=24, w=16)
grid.arrange(grobs=plot_list, ncol=2)
dev.off()


### boxplot of correlations

psi <- unlist(as.list(data.frame(acorr)))
gene <- unlist(as.list(data.frame(gcorr)))

svg("Figure2C_boxplot.svg", h=3, w=6)
boxplot(gene, psi, outline=F, ylim=c(0,1), horizontal=TRUE)
dev.off()
t.test(psi, gene, alternative = c("two.sided")) #p-value < 2.2e-16
data <- cbind(data.frame(psi), data.frame(gene))

## boxplot of correlations for genes with more than 1 isoforms

psi <- unlist(as.list(data.frame(a_filteredcorr)))
gene <- unlist(as.list(data.frame(gcorr)))

svg("Figure2C_boxplot_morethan1isoform.svg", h=3, w=6)
boxplot(gene, psi, outline=F, ylim=c(0,1), horizontal=TRUE)
dev.off()
t.test(psi, gene, alternative = c("two.sided")) #p-value < 2.2e-16
data <- cbind(data.frame(psi), data.frame(gene))




### alternative isoforms corr vs main isoform
main_iso <- subset(a_filtered, grepl("(.1)$", rownames(a_filtered)))
alternative_iso <- subset(a_filtered, !grepl("(.1)$", rownames(a_filtered)))

main_isocorr <-cor(main_iso,method="spearman")
alternative_isocorr <-cor(alternative_iso,method="spearman")

main <- unlist(as.list(data.frame(main_isocorr)))
alt <- unlist(as.list(data.frame(alternative_isocorr)))

svg("2025_s1_boxplot_alternativeisoforms.svg", h=6, w=5)
boxplot(psi, main, alt, outline=F, ylim=c(0,1), horizontal=F, names=c("all isoforms","main isoforms","alternative isoforms"), 
        ylab="Spearman's rho correlation PSI")
dev.off()


t.test(psi, alt, alternative = c("greater"))
t.test(psi, main, alternative = c("greater"))
t.test(main, alt, alternative = c("greater"))



### boxplot per lake
is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
a_filtered[is.nan.data.frame(a_filtered)] <- 0
a_filteredcorr <-cor(a_filtered,method="spearman")

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
a[is.nan.data.frame(a)] <- 0
accorr <-cor(a,method="spearman")

lv <- anno[anno$lake == "LV",]
lvg <- g[,colnames(g) %in% lv$ids]
lva <- a[,colnames(a) %in% lv$ids]
lvaf <- a_filtered[,colnames(a_filtered) %in% lv$ids]

lm <- anno[anno$lake == "LM",]
lmg <- g[,colnames(g) %in% lm$ids]
lma <- a[,colnames(a) %in% lm$ids]
lmaf <- a_filtered[,colnames(a_filtered) %in% lm$ids]

lt <- anno[anno$lake == "LT",]
ltg <- g[,colnames(g) %in% lt$ids]
lta <- a[,colnames(a) %in% lt$ids]
ltaf <- a_filtered[,colnames(a_filtered) %in% lt$ids]


rv <- anno[anno$lake == "RV",]
rvg <- g[,colnames(g) %in% rv$ids]
rva <- a[,colnames(a) %in% rv$ids]
rvaf <- a_filtered[,colnames(a_filtered) %in% rv$ids]


lvcorr <-cor(lvg,method="spearman")
lmcorr <-cor(lmg,method="spearman")
ltcorr <-cor(ltg,method="spearman")
rvcorr <-cor(rvg,method="spearman")

genelv <- unlist(as.list(data.frame(lvcorr)))
genelm <- unlist(as.list(data.frame(lmcorr)))
genelt <- unlist(as.list(data.frame(ltcorr)))
generv <- unlist(as.list(data.frame(rvcorr)))


lvcorra <-cor(lva,method="spearman")
lmcorra <-cor(lma,method="spearman")
ltcorra <-cor(lta,method="spearman")
rvcorra <-cor(rva,method="spearman")

psilv <- unlist(as.list(data.frame(lvcorra)))
psilm <- unlist(as.list(data.frame(lmcorra)))
psilt <- unlist(as.list(data.frame(ltcorra)))
psirv <- unlist(as.list(data.frame(rvcorra)))

lvcorraf <-cor(lvaf,method="spearman")
lmcorraf <-cor(lmaf,method="spearman")
ltcorraf <-cor(ltaf,method="spearman")
rvcorraf <-cor(rvaf,method="spearman")

psilvf <- unlist(as.list(data.frame(lvcorraf)))
psilmf <- unlist(as.list(data.frame(lmcorraf)))
psiltf <- unlist(as.list(data.frame(ltcorraf)))
psirvf <- unlist(as.list(data.frame(rvcorraf)))

svg("Figure2B_supp_boxplot_lakes.svg")
boxplot(generv,genelv, genelm, genelt, psirv, psilv, psilm, psilt, psirvf, psilvf, psilmf, psiltf, outline=F, 
        ylim=c(0,1), horizontal=F,  names=c("NR", "LV","LM","LT","NR", "LV","LM","LT","NR", "LV","LM","LT"), 
        ylab="Spearman's rho correlation", col=c("grey", "grey", "grey", "grey", "white", "white", "white","white", "beige", "beige","beige","beige"))
text(x=3,y=0.9, "gene expression")
text(x=6,y=0.72, "PSI all genes")
text(x=10, y=0.7, "PSI only spliced genes")
dev.off()

t.test(psilt, psilm, alternative = c("less"))


