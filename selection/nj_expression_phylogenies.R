

## pooja.singh09@gmail.com
## gene expression vs isoform expression tree
####1 gene expression OJ vs PJ



library(edgeR)
library(ouch)
library(reshape2)
library(ggplot2)
library(ape)

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution')

############################################################################################################# gene exp
####### OJ gene expression treee


exp <- read.table("gene_tpm_all_samples_onil_OJevee.tsv", header=T, comment.char="", row.names=1, stringsAsFactors = F)
index <- read.table("gene_tpm_all_samples_onil_OJevee.index.txt", header=F, comment.char="")
phylogeny <- read.delim("st26_phylogeny_ouch.txt")


##get average of species
data=exp[2:ncol(exp)]
data_mean = c()
for (i in unique(index[,1])) {
        curSpecies = data[,index[,1] == i]
        if (sum(index[,1] == i) == 1) {
                data_mean = cbind(data_mean, curSpecies)
        } else {
                data_mean = cbind(data_mean, apply(curSpecies, 1, function(x) mean(x, na.rm=T)))
        }
}
colnames(data_mean) = unique(index[,1])


##normalize data##
data_mean.noNA = data_mean
data_mean.noNA[is.na(data_mean)] = 0
scale = calcNormFactors(data_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]
data_mean.norm = sapply(seq(1, length(scale)), function(x) data_mean[,x] / scale[x])
data_mean.norm = log10(data_mean.norm+0.01)
colnames(data_mean.norm) = colnames(data_mean)
print(head(data_mean.norm))

# calculate gene exp distance between samples

dist <- cor(data_mean.norm, method="spearman", use="complete.obs")
dist1 <- 1-dist
ojtree <- ape::nj(dist1)
ojtree <- root(ojtree, 1)

# calculate tree

svg("gene_expression_nj_tree_oj.svg")
plot(ojtree)
add.scale.bar(0.1,0.4,cex=0.6)
dev.off()

ape::write.tree(ojtree, file='gene_expression_nj_tree_oj.tre') 


#######  PJ  gene expression treee

exp <- read.table("gene_tpm_all_samples_onil_PJevee.tsv", header=T, comment.char="", row.names=1, stringsAsFactors = F)
index <- read.table("gene_tpm_all_samples_onil_PJevee.index.txt", header=F, comment.char="")
phylogeny <- read.delim("st26_phylogeny_ouch.txt")


##get average of species
data=exp[2:ncol(exp)]
data_mean = c()
for (i in unique(index[,1])) {
        curSpecies = data[,index[,1] == i]
        if (sum(index[,1] == i) == 1) {
                data_mean = cbind(data_mean, curSpecies)
        } else {
                data_mean = cbind(data_mean, apply(curSpecies, 1, function(x) mean(x, na.rm=T)))
        }
}
colnames(data_mean) = unique(index[,1])




##normalize data##
data_mean.noNA = data_mean
data_mean.noNA[is.na(data_mean)] = 0
scale = calcNormFactors(data_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]
data_mean.norm = sapply(seq(1, length(scale)), function(x) data_mean[,x] / scale[x])
data_mean.norm = log10(data_mean.norm+0.01)
colnames(data_mean.norm) = colnames(data_mean)
print(head(data_mean.norm))

# calculate gene exp distance between samples

dist <- cor(data_mean.norm, method="spearman", use="complete.obs")
dist1 <- 1-dist
pjtree <- ape::nj(dist1)
pjtree <- root(pjtree, 1)

# calculate tree

svg("gene_expression_nj_tree_pj.svg")
plot(pjtree)
add.scale.bar(0.1,0.4,cex=0.6)
dev.off()

ape::write.tree(pjtree, file='gene_expression_nj_tree_pj.tre')


######################################################################################################################################################## isoform exp


all <- read.table("transcript_tpm_all_samples_onil.tsv", header=T,  stringsAsFactors = F)
oj <- all[,grepl("OJ", colnames(all))]
pj <- all[,grepl("PJ", colnames(all))]

oj1 <- cbind(all$On20d3,oj)
oj1 <- cbind(all$Transcript_ID,oj1)
colnames(oj1)[1] <- "GENE_NAME"
colnames(oj1)[2] <- "On"
rownames(oj1) <- paste0("OG",1:nrow(oj1))

pj1 <- cbind(all$On20d3,pj)
pj1 <- cbind(all$Transcript_ID,pj1)
colnames(pj1)[1] <- "GENE_NAME"
colnames(pj1)[2] <- "On"
rownames(pj1) <- paste0("OG",1:nrow(pj1))

#######OJ
exp <- oj1
index <- read.table("gene_tpm_all_samples_onil_OJevee.index.txt", header=F, comment.char="")
phylogeny <- read.delim("st26_phylogeny_ouch.txt")

##get average of species
data=exp[2:ncol(exp)]
data_mean = c()
for (i in unique(index[,1])) {
        curSpecies = data[,index[,1] == i]
        if (sum(index[,1] == i) == 1) {
                data_mean = cbind(data_mean, curSpecies)
        } else {
                data_mean = cbind(data_mean, apply(curSpecies, 1, function(x) mean(x, na.rm=T)))
        }
}
colnames(data_mean) = unique(index[,1])


##normalize data##
data_mean.noNA = data_mean
data_mean.noNA[is.na(data_mean)] = 0
scale = calcNormFactors(data_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]
data_mean.norm = sapply(seq(1, length(scale)), function(x) data_mean[,x] / scale[x])
data_mean.norm = log10(data_mean.norm+0.01)
colnames(data_mean.norm) = colnames(data_mean)
print(head(data_mean.norm))


# calculate gene exp distance between samples

dist <- cor(data_mean.norm, method="spearman", use="complete.obs")
dist1 <- 1-dist
ojtree <- ape::nj(dist1)
ojtree <- root(ojtree, 1)

# calculate tree

svg("isoform_expression_nj_tree_oj.svg")
plot(ojtree)
add.scale.bar(0.1,0.4,cex=0.6)
dev.off()

ape::write.tree(ojtree, file='isoform_expression_nj_tree_oj.tre')


#######PJ

exp <- pj1
index <- read.table("gene_tpm_all_samples_onil_PJevee.index.txt", header=F, comment.char="")
phylogeny <- read.delim("st26_phylogeny_ouch.txt")


##get average of species
data=exp[2:ncol(exp)]
data_mean = c()
for (i in unique(index[,1])) {
        curSpecies = data[,index[,1] == i]
        if (sum(index[,1] == i) == 1) {
                data_mean = cbind(data_mean, curSpecies)
        } else {
                data_mean = cbind(data_mean, apply(curSpecies, 1, function(x) mean(x, na.rm=T)))
        }
}
colnames(data_mean) = unique(index[,1])


##normalize data##
data_mean.noNA = data_mean
data_mean.noNA[is.na(data_mean)] = 0
scale = calcNormFactors(data_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]
data_mean.norm = sapply(seq(1, length(scale)), function(x) data_mean[,x] / scale[x])
data_mean.norm = log10(data_mean.norm+0.01)
colnames(data_mean.norm) = colnames(data_mean)
print(head(data_mean.norm))


# calculate gene exp distance between samples

dist <- cor(data_mean.norm, method="spearman", use="complete.obs")
dist1 <- 1-dist
pjtree <- ape::nj(dist1)
pjtree <- root(pjtree, 1)

# calculate tree

svg("isoform_expression_nj_tree_pj.svg")
plot(pjtree)
add.scale.bar(0.1,0.4,cex=0.6)
dev.off()

ape::write.tree(pjtree, file='isoform_expression_nj_tree_pj.tre')
