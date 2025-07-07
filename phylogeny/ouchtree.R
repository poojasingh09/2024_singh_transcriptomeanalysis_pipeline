library(ouch)
library(ape)

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/12_phylogeny')
filename <- "RAxML_bipartitions.all.filt.L.vcf.SNP.gz.recode.mergedOnil.fmiss.raxmlout.T10.clean.tre"
tree <- ape::read.tree(filename)
ouchtree <- ape2ouch(tree, scale = TRUE, branch.lengths = tree$edge.length)

sink('st26_phylogeny_ouch.txt')
ouchtree
sink()


