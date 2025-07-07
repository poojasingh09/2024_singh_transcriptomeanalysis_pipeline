#pooja.singh09@gmail.com
#Singh et al
#plot pca, hclust, and heatmap
#udated september 2021

library(tidyverse)
library(WGCNA)
library(pheatmap)
library(tidyr)
library(dplyr)
packageVersion("WGCNA")
packageVersion("pheatmap")
sessionInfo()
library(cowplot)



annapalette<- c("LM"="goldenrod2", "LT"="#009E73", "LV"="#E49EC2", "RV"="#56B4E9")

## colours: http://sape.inf.usi.ch/quick-reference/ggplot2/colour


########## read in input data ##########

anno <- read.table("/Users/pooja/Desktop/Papers/2018_Postdoc_papers/RNAseq_newgrant/allst26/analysis/new_ref/02_DEG/phenotypes.txt", header=T)
anno1 <- data.frame(cbind(anno$lake, anno$jaw))
rownames(anno1) <- anno$ids
colnames(anno1) <- c("lake", "jaw")

g <- read.table("/Users/pooja/Desktop/Papers/2018_Postdoc_papers/RNAseq_newgrant/allst26/analysis/new_ref/02_DEG/deseq_new/DESeq_NormalizedCounts_vsd.txt", header=T, row.names=1)
t <- read.table("/Users/pooja/Desktop/Papers/2018_Postdoc_papers/RNAseq_newgrant/allst26/analysis/new_ref/03_DTU/DESeq_NormalizedCounts_vsd.txt", header=T, row.names=1)

g1 <- g[rowSums(g<5)<100, ]
t1 <- t[rowSums(t<5)<100, ]

#g1 <- sample_n(g1, 2000)
#t1 <- sample_n(t1, 2000)


##### FIGURE 1A PCA ###################

g1t <- g %>% as.matrix()  %>% t() # transpose the matrix so that rows = samples and columns = variables



gene_pca <- prcomp(g1t)
gene_pc_eigenvalues <- gene_pca$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
gene_pc_eigenvalues <- tibble(PC = factor(1:length(gene_pc_eigenvalues)), 
                         variance = gene_pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
gene_pc_eigenvalues

gene_pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

gene_pc_scores <- data.frame(gene_pca$x)


anno1$sample <- rownames(anno1)

# get the PC scores from prcomp object
gene_pca_plot <- gene_pca$x %>% 
  # convert it to a tibble
  as_tibble(rownames = "sample") %>% 
  # join with "anno1" table
  merge(anno1, by = "sample") %>% 
  # make the plot
  ggplot(aes(x = PC1, y = PC2, 
             color = lake, shape = jaw)) +
  xlab("PC1 (6.8 %)")+
  ylab("PC2 (4.2 %)")+
  geom_point(size=2, alpha = 0.80) 

#plot1 <- gene_pca_plot + scale_color_manual(values=c("goldenrod2", "seagreen3", "tomato2", "steelblue2")) + theme_classic()
plot1 <- gene_pca_plot + scale_color_manual(values=annapalette) + theme_classic()



##### FIGURE 1B PCA ###################

t1t <- t %>% as.matrix()  %>% t() # transpose the matrix so that rows = samples and columns = variables



transcript_pca <- prcomp(t1t)
transcript_pc_eigenvalues <- transcript_pca$sdev^2

# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
transcript_pc_eigenvalues <- tibble(PC = factor(1:length(transcript_pc_eigenvalues)), 
                         variance = transcript_pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
transcript_pc_eigenvalues

transcript_pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

transcript_pc_scores <- data.frame(transcript_pca$x)


anno1$sample <- rownames(anno1)

# get the PC scores from prcomp object
transcript_pca_plot <- transcript_pca$x %>% 
  # convert it to a tibble
  as_tibble(rownames = "sample") %>% 
  # join with "anno1" table
  merge(anno1, by = "sample") %>% 
  # make the plot
  ggplot(aes(x = PC1, y = PC2, 
             color = lake, shape = jaw)) +
  xlab("PC1 (7.5 %)")+
  ylab("PC2 (4.4 %)")+
  geom_point(size=2, alpha = 0.80) 

#plo2 <- transcript_pca_plot + scale_color_manual(values=c("goldenrod2", "seagreen3", "tomato2", "steelblue2")) + theme_classic()
plot2 <- transcript_pca_plot + scale_color_manual(values=annapalette) + theme_classic()



### FINAL FIG 1 ###


svg("figure1_PCA.svg", h=5, w=12)
plot_grid(plot1, plot2, labels = c('A', 'B'), label_size = 12)
dev.off()

pdf("figure1_PCA.pdf", h=5, w=12)
plot_grid(plot1, plot2, labels = c('A', 'B'), label_size = 12)
dev.off()





