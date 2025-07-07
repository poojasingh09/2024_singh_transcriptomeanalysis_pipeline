


#pooja.singh09@gmail.com
# running evee tools

## input files
# index files should have 2 cols, col1: species and col2: sample id
# gene expression files should have OG, geneID, expIND1, expIND2 etc
# sepcies tree convert newick to ouch format: 

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/selection/EVEEtools')

############### GENE EXPRESSION

#1 get residuals with Oniloticus as reference

Rscript --vanilla ../EVEE-Tools/residuals.R gene_tpm_all_samples_onil_OJevee.tsv gene_tpm_all_samples_onil_OJevee.index.txt On OJ_residuals_gene.txt &
Rscript --vanilla ../EVEE-Tools/residuals.R gene_tpm_all_samples_onil_PJevee.tsv gene_tpm_all_samples_onil_PJevee.index.txt On PJ_residuals_gene.txt &


########## in R prepare phylogeny to ouch format
########## the species column in the tree must be NA for all inner nodes. 

library(ouch)
library(ape)
filename <- "/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/12_phylogeny/RAxML_bipartitions.all.filt.L.vcf.SNP.gz.recode.mergedOnil.fmiss.raxmlout.T10.clean.tre"
tree <- ape::read.tree(filename)
ouchtree <- ape2ouch(tree, scale = FALSE, branch.lengths = tree$edge.length)

ouch.df <- as(ouchtree, "data.frame")
colnames(ouch.df)[3:4] <- c("time", "species")
leaf.idx <- grep("[A-z]+", ouch.df$species)
ouch.df$species[-leaf.idx] <- NA
write.table(ouch.df, "st26_phylogeny_ouch1.txt", sep="\t", quote=FALSE, row.names=FALSE)

######## exit R

#2 fit evolutionary mean and variance for each gene
# Significance of OU fit ("qvalues"), evolutionary mean ("thetas") and variance ("var") under OU model, and rate of divergence under Brownian motion model ("brownSigmas") for each gene.
# significant qvalues indicate gene is evolving non-neutrally


 &
Rscript --vanilla ../../EVEE-Tools/fitOUModel.R gene_tpm_all_samples_onil_PJevee.tsv gene_tpm_all_samples_onil_PJevee.index.txt st26_phylogeny_ouch.txt PJ_ou_stats_gene.txt &


Rscript --vanilla ../../EVEE-Tools/fitOUModel.R gene_tpm_all_samples_onil_OJevee.tsv gene_tpm_all_samples_onil_OJevee.index.txt RAxML_bipartitions.all.filt.L.vcf.SNP.gz.recode.mergedOnil.fmiss.raxmlout.T10.clean.CAGE.ouch.tre.txt OJ_ou_stats_a.txt

#3 identify diff exp across tree
# This script is a wrapper script for ouch to fit multivariate OU distributions to expression data to identifying differential gene expression across a tree 
# we implement 5 regimes: null, lt vs all, lm vs all, lv vs all, herb vs others, carn vs others, lacustrine vs riverine
# regime file should have 1 row for every row in the phylogeny file

Rscript --vanilla ../../EVEE-Tools/ouRegimes.R gene_tpm_all_samples_onil_OJevee.tsv gene_tpm_all_samples_onil_OJevee.index.txt st26_phylogeny_ouch.txt regimes_2.txt OJ_regimes_diffexp_regime2_gene &
Rscript --vanilla ../../EVEE-Tools/ouRegimes.R gene_tpm_all_samples_onil_PJevee.tsv gene_tpm_all_samples_onil_PJevee.index.txt st26_phylogeny_ouch.txt regimes_2.txt PJ_regimes_diffexp_regime2_gene &


#3.1 to identify regimes plot node labels in R for the phlyogeny

library(ape)
library(ouch)

ouch <- read.delim("st26_phylogeny_ouch1.txt")
tree = ouchtree(ouch$node, ouch$ancestor, ouch$time, labels = as.character(ouch$species))
plot(tree, node.names=TRUE, text_opts=list(cex=0.4))
svg('st26_phylogeny_ouch1.svg')
plot(tree)
dev.off()


############## PSI

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/selection/EVEEtools')

#1 get residuals with Oniloticus as reference

Rscript --vanilla ../../EVEE-Tools/residuals.R isoform_psi_all_samples_onil_OJevee.tsv isoform_psi_all_samples_onil_OJevee.index.txt On OJ_residuals_psi.txt &
Rscript --vanilla ../../EVEE-Tools/residuals.R isoform_psi_all_samples_onil_PJevee.tsv isoform_psi_all_samples_onil_PJevee.index.txt On PJ_residuals_psi.txt &
  
  
######## exit R

#2 fit evolutionary mean and variance for each transcript
# Significance of OU fit ("qvalues"), evolutionary mean ("thetas") and variance ("var") under OU model, and rate of divergence under Brownian motion model ("brownSigmas") for each transcript.
# significant qvalues indicate transcript is evolving non-neutrally

Rscript --vanilla ../../EVEE-Tools/fitOUModel.R isoform_psi_all_samples_onil_OJevee.tsv isoform_psi_all_samples_onil_OJevee.index.txt st26_phylogeny_ouch.txt OJ_ou_stats_psi.txt &
Rscript --vanilla ../../EVEE-Tools/fitOUModel.R isoform_psi_all_samples_onil_PJevee.tsv isoform_psi_all_samples_onil_PJevee.index.txt st26_phylogeny_ouch.txt PJ_ou_stats_psi.txt &
  
  
#3 identify diff exp across tree
# This script is a wrapper script for ouch to fit multivariate OU distributions to expression data to identifying differential transcript expression across a tree 
# we implement 5 regimes: null, lt vs all, lm vs all, lv vs all, herb vs others, carn vs others, lacustrine vs riverine
# regime file should have 1 row for every row in the phylogeny file

Rscript --vanilla ../../EVEE-Tools/ouRegimes.R isoform_psi_all_samples_onil_OJevee.tsv isoform_psi_all_samples_onil_OJevee.index.txt st26_phylogeny_ouch.txt regimes_2.txt OJ_regimes_diffexp_regime2_psi &
Rscript --vanilla ../../EVEE-Tools/ouRegimes.R isoform_psi_all_samples_onil_PJevee.tsv isoform_psi_all_samples_onil_PJevee.index.txt st26_phylogeny_ouch.txt regimes_2.txt PJ_regimes_diffexp_regime2_psi &
  
  
  #3.1 to identify regimes plot node labels in R for the phlyogeny
  
library(ape)
library(ouch)

ouch <- read.delim("st26_phylogeny_ouch.txt")
tree = ouchtree(ouch$node, ouch$ancestor, ouch$time, labels = as.character(ouch$species))
plot(tree, node.names=TRUE, text_opts=list(cex=0.2))




############## ISOFORM EXPRESSION


#1 get residuals with Oniloticus as reference

Rscript --vanilla ../../EVEE-Tools/residuals.R transcript_tpm_all_samples_onil_OJevee.tsv transcript_tpm_all_samples_onil_OJevee.index.txt On OJ_residuals_transcript.txt &
Rscript --vanilla ../../EVEE-Tools/residuals.R transcript_tpm_all_samples_onil_PJevee.tsv transcript_tpm_all_samples_onil_PJevee.index.txt On PJ_residuals_transcript.txt &
  
  
  ######## exit R
  
  #2 fit evolutionary mean and variance for each transcript
  # Significance of OU fit ("qvalues"), evolutionary mean ("thetas") and variance ("var") under OU model, and rate of divergence under Brownian motion model ("brownSigmas") for each transcript.
  # significant qvalues indicate transcript is evolving non-neutrally
  
  Rscript --vanilla ../../EVEE-Tools/fitOUModel.R transcript_tpm_all_samples_onil_OJevee.tsv transcript_tpm_all_samples_onil_OJevee.index.txt st26_phylogeny_ouch.txt OJ_ou_stats_transcript.txt &
  Rscript --vanilla ../../EVEE-Tools/fitOUModel.R transcript_tpm_all_samples_onil_PJevee.tsv transcript_tpm_all_samples_onil_PJevee.index.txt st26_phylogeny_ouch.txt PJ_ou_stats_transcript.txt &
  
  
  #3 identify diff exp across tree
  # This script is a wrapper script for ouch to fit multivariate OU distributions to expression data to identifying differential transcript expression across a tree 
  # we implement 5 regimes: null, lt vs all, lm vs all, lv vs all, herb vs others, carn vs others, lacustrine vs riverine
  # regime file should have 1 row for every row in the phylogeny file
  
  Rscript --vanilla ../../EVEE-Tools/ouRegimes.R transcript_tpm_all_samples_onil_OJevee.tsv transcript_tpm_all_samples_onil_OJevee.index.txt st26_phylogeny_ouch.txt regimes_2.txt OJ_regimes_diffexp_regime2_transcript &
  Rscript --vanilla ../../EVEE-Tools/ouRegimes.R transcript_tpm_all_samples_onil_PJevee.tsv transcript_tpm_all_samples_onil_PJevee.index.txt st26_phylogeny_ouch.txt regimes_2.txt PJ_regimes_diffexp_regime2_transcript &
  
  
  #3.1 to identify regimes plot node labels in R for the phlyogeny
  
  library(ape)
library(ouch)

ouch <- read.delim("st26_phylogeny_ouch.txt")
tree = ouchtree(ouch$node, ouch$ancestor, ouch$time, labels = as.character(ouch$species))
plot(tree, node.names=TRUE, text_opts=list(cex=0.2))



