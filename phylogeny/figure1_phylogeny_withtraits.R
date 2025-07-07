## st26
## pooja.singh09@gmail.com
## figure 1 phylogeny


setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/12_phylogeny')

## libs

library(diversitree)
library(ape)
library(ggtree)
library(ggnewscale)
library(geiger)
library(treeio)
library(ggplot2)
library(phylobase)
library(adephylo)
library(wesanderson)
library(dplyr)
library(gridExtra)
library(grid)

## read in tree

tree <- ape::read.nexus("RAxML_bipartitions.all.filt.L.vcf.SNP.gz.recode.mergedOnil.fmiss.raxmlout.T10.clean.tre")
plot(tree,show.node=TRUE)
listTips(tree)
is.rooted(tree) 


## annotat
anno <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/phenotypes.txt", header=T)

dat1 <- anno[anno$jaw == "PJ",]
dat0 <- unique(data.frame(dat1$diet1, dat1$species))
dim(dat0)
rownames(dat0) <- dat0$dat1.species
colnames(dat0) <- c("diet","species")
dat <- dat0

##Plot
circ <- ggtree(tree, layout = "rectangular") +
  geom_point(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 50, 
                 x = x, y = y, size = as.numeric(label)), 
             shape = 21, fill = "black", color = "black", alpha = 0.4) +
  scale_size_continuous(range = c(1, 3)) +  # Reduce circle sizes
  theme(legend.position = "right")

print(circ)

circ_p1 <- gheatmap(circ, dat0[1], offset = 0.01, width = .06, 
                    colnames_angle = 95, colnames_offset_y = 0.2, colnames = F) + 
                    scale_fill_manual(values = c( herbivore = "green",  carnivore= "red", omnivore="grey90"), na.value = "grey90")


circ_p1

circ_p1 + 
  theme(legend.position = "none") + 
  geom_tiplab(aes(x = 0.33, y = y), fontface = "italic", size = 3, inherit.aes = FALSE)

ggsave("st26_phylogeny_fig1.svg", device = svglite::svglite,  width=10, height=10)


##Plot circ
circ <- ggtree(tree, layout = "circular", branch.length = "branch.length") 

print(circ)

circ_p1 <- gheatmap(circ, dat0[1], offset = 0.01, width = .06, 
                    colnames_angle = 95, colnames_offset_y = 0.2, colnames = F) + 
  scale_fill_manual(values = c( herbivore = "green",  carnivore= "red", omnivore="grey90"), na.value = "grey90")


circ_p1

circ_p1 + 
  theme(legend.position = "none") + 
  geom_tiplab(aes(x = 0.33, y = y), fontface = "italic", size = 3, inherit.aes = FALSE)

ggsave("st26_phylogeny_fig1_circ.svg", device = svglite::svglite,  width=10, height=10)


