

##pooja.singh09@gmail.com
## modularity gene expression tanglegram
## caculate mean of gene expression and make NJ distance tree on 1-spearmans rho and then compare OJ vs PJ phylogeny using tanglegram

library(edgeR)
library(ouch)
library(reshape2)
library(ggplot2)
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(ape))
library(phylogram)
library(dplyr)


setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution')

######################################### mean gene exp per species read in

oj <- read.table("gene_tpm_all_samples_onil_OJ_mean_norm.txt", header=T, comment.char="", row.names=1, stringsAsFactors = F)
pj <- read.table("gene_tpm_all_samples_onil_PJ_mean_norm.txt", header=T, comment.char="", row.names=1, stringsAsFactors = F)

######### spearmans correlation

ojcor <- cor(oj, method="spearman", use = "pairwise.complete.obs")
pjcor <- cor(pj, method="spearman", use = "pairwise.complete.obs")

ojdist <- 1-ojcor
pjdist <- 1-pjcor

pjtree <- nj(pjdist)
pjtree <- root(pjtree, 1)

ojtree <- nj(ojdist)
ojtree <- root(ojtree, 1)


dend1 <- as.dendrogram (ojtree)
dend2 <- as.dendrogram (pjtree)


## plot the tanglegram
# Align and plot two dendrograms side by side

svg("tanglegram_gen_exp_norm.svg", h=10, w=7) 
# Assuming dend1 and dend2 are your dendrogram objects
# Modify the text size of the labels
dend1 <- set(dend1, "labels_cex", 2.5)  # Adjust this value as needed
dend2 <- set(dend2, "labels_cex", 2.5)  # Adjust this value as needed

# Set graphical parameters for axis label size
par(cex.axis = 2)  # Adjust this value to increase x-axis label size

# Create the tanglegram
dendlist(dend1, dend2) %>%
  untangle(method = "step1side") %>%
  tanglegram(
    lwd = 1,
    edge.lwd = 1
  )

# Reset graphical parameters to default
par(cex.axis = 1)
dev.off()     


######################################### mean isoform exp per species read in

oj <- read.table("transcript_tpm_all_samples_onil_OJ_mean_norm.txt", header=T, comment.char="", row.names=1, stringsAsFactors = F)
pj <- read.table("transcript_tpm_all_samples_onil_PJ_mean_norm.txt", header=T, comment.char="", row.names=1, stringsAsFactors = F)

######### spearmans correlation

ojcor <- cor(oj, method="spearman", use = "pairwise.complete.obs")
pjcor <- cor(pj, method="spearman", use = "pairwise.complete.obs")

ojdist <- 1-ojcor
pjdist <- 1-pjcor

pjtree <- nj(pjdist)
pjtree <- root(pjtree, 1)

ojtree <- nj(ojdist)
ojtree <- root(ojtree, 1)


dend1 <- as.dendrogram (ojtree)
dend2 <- as.dendrogram (pjtree)


## plot the tanglegram
# Align and plot two dendrograms side by side

svg("tanglegram_transcript_exp_norm.svg") 
dendlist(dend1, dend2) %>% untangle(method = "step1side") %>% tanglegram( lwd = 1, edge.lwd =1)
dev.off()  



###################################################################################################################################################### non-norm

######################################### mean gene exp per species read in

oj <- read.table("gene_tpm_all_samples_onil_OJ_mean.txt", header=T, comment.char="", row.names=1, stringsAsFactors = F)
pj <- read.table("gene_tpm_all_samples_onil_PJ_mean.txt", header=T, comment.char="", row.names=1, stringsAsFactors = F)


######### spearmans correlation

ojcor <- cor(oj, method="spearman", use = "pairwise.complete.obs")
pjcor <- cor(pj, method="spearman", use = "pairwise.complete.obs")

ojdist <- 1-ojcor
pjdist <- 1-pjcor

pjtree <- nj(pjdist)
pjtree <- root(pjtree, 1)

ojtree <- nj(ojdist)
ojtree <- root(ojtree, 1)


dend1 <- as.dendrogram (ojtree)
dend2 <- as.dendrogram (pjtree)


## plot the tanglegram
# Align and plot two dendrograms side by side

svg("tanglegram_gen_eexp.svg") 
dendlist(dend1, dend2) %>% untangle(method = "step1side") %>% tanglegram( lwd = 1, edge.lwd =1)
dev.off()     


######################################### mean isoform exp per species read in

oj <- read.table("transcript_tpm_all_samples_onil_OJ_mean.txt", header=T, comment.char="", row.names=1, stringsAsFactors = F)
pj <- read.table("transcript_tpm_all_samples_onil_PJ_mean.txt", header=T, comment.char="", row.names=1, stringsAsFactors = F)

######### spearmans correlation

ojcor <- cor(oj, method="spearman", use = "pairwise.complete.obs")
pjcor <- cor(pj, method="spearman", use = "pairwise.complete.obs")

ojdist <- 1-ojcor
pjdist <- 1-pjcor

pjtree <- nj(pjdist)
pjtree <- root(pjtree, 1)

ojtree <- nj(ojdist)
ojtree <- root(ojtree, 1)


dend1 <- as.dendrogram (ojtree)
dend2 <- as.dendrogram (pjtree)


## plot the tanglegram
# Align and plot two dendrograms side by side

svg("tanglegram_transcript_exp.svg") 
dendlist(dend1, dend2) %>% untangle(method = "step1side") %>% tanglegram( lwd = 1, edge.lwd =1)
dev.off()



############# tanglegram of genes convergent for OJ and PJ

a <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/PJOJ_convergent_genes_up_down.ourid.genecount")

oj <- a[,grepl("OJ", colnames(a))]
dim(oj)

pj <- a[,grepl("PJ", colnames(a))]
dim(pj)


ojcor <- cor(oj, method="spearman", use = "pairwise.complete.obs")
pjcor <- cor(pj, method="spearman", use = "pairwise.complete.obs")

ojdist <- 1-ojcor
pjdist <- 1-pjcor

pjtree <- nj(pjdist)
pjtree <- root(pjtree, 1)

ojtree <- nj(ojdist)
ojtree <- root(ojtree, 1)


dend1 <- as.dendrogram (ojtree)
dend2 <- as.dendrogram (pjtree)


## plot the tanglegram
# Align and plot two dendrograms side by side

svg("tanglegram_gen_exp_norm.svg") 
dendlist(dend1, dend2) %>% untangle(method = "step1side") %>% tanglegram( lwd = 1, edge.lwd =1)
dev.off()     



