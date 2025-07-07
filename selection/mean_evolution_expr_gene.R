################################################################################################# mean gene expression change over evolutionary time

## pooja.singh09@gmail.com
## april 2023
## investigateing evolution change in gene expression
## caculate mean of gene expression and correlate with evolution time (branch length)

library(edgeR)
library(ouch)
library(reshape2)
library(ggplot2)

setwd("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution")
getwd()
#######OJ
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
#write.table(data_mean,"gene_tpm_all_samples_onil_OJ_mean.txt", quote=F, sep="\t")

##normalize data##
data_mean.noNA = data_mean
data_mean.noNA[is.na(data_mean)] = 0
scale = calcNormFactors(data_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]
data_mean.norm = sapply(seq(1, length(scale)), function(x) data_mean[,x] / scale[x])
data_mean.norm = log10(data_mean.norm+0.01)
colnames(data_mean.norm) = colnames(data_mean)
print(head(data_mean.norm))
#write.table(data_mean.norm,"gene_tpm_all_samples_onil_OJ_mean_norm.txt", quote=F, sep="\t")

# calculate gene exp distance between samples

dist <- cor(data_mean.norm, method="spearman", use="complete.obs")
dist1 <- 1-dist
exp_dist_on <- dist1[,1]
exp_dist_on <- data.frame(exp_dist_on)
#exp_dist_on  <- exp_dist_on[order(row.names(exp_dist_on)), ]

p <- phylogeny[46:66,]
p1 <- p[,c(2,4)]
row.names(p1) <- p1$species
p1 <- p1[order(row.names(p1)), ]

out <- merge(p1, exp_dist_on, by="row.names")
out1 <- out[order(out$time), ]
out_oj <- out1



#######PJ
exp <- read.table("gene_tpm_all_samples_onil_PJevee.tsv", header=T, comment.char="", row.names=1, stringsAsFactors = F)
index <- read.table("gene_tpm_all_samples_onil_PJevee.index.txt", header=F, comment.char="")
#phylogeny <- read.delim("st26_phylogeny_ouch1.txt")


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
#write.table(data_mean,"gene_tpm_all_samples_onil_PJ_mean.txt", quote=F, sep="\t")

##normalize data##
data_mean.noNA = data_mean
data_mean.noNA[is.na(data_mean)] = 0
scale = calcNormFactors(data_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]
data_mean.norm = sapply(seq(1, length(scale)), function(x) data_mean[,x] / scale[x])
data_mean.norm = log10(data_mean.norm+0.01)
colnames(data_mean.norm) = colnames(data_mean)
print(head(data_mean.norm))
#write.table(data_mean.norm,"gene_tpm_all_samples_onil_PJ_mean_norm.txt", quote=F, sep="\t")

# calculate gene exp distance between samples

dist <- cor(data_mean.norm, method="spearman", use="complete.obs")
dist1 <- 1-dist
exp_dist_on <- dist1[,1]
exp_dist_on <- data.frame(exp_dist_on)
#exp_dist_on  <- exp_dist_on[order(row.names(exp_dist_on)), ]

p <- phylogeny[46:66,]
p1 <- p[,c(2,4)]
row.names(p1) <- p1$species
p1 <- p1[order(row.names(p1)), ]


out <- merge(p1, exp_dist_on, by="row.names")
out1 <- out[order(out$time), ]
out_pj <- out1

## plot 1-p gene expression distance vs time (branch length)
library(ggplot2)
out3 <- merge(out_oj, out_pj, by="species")

Data = data.frame(out3$time.x, out3$exp_dist_on.x, out3$exp_dist_on.y)
colnames(Data) <- c("time","oj","pj")

svg(filename="gene_expression_divergence_over_time.svg")
ggplot(Data) +
	geom_smooth(aes(Data$time,oj),colour="black",fill="black") +
	geom_smooth(aes(Data$time,pj),colour="grey",fill="grey") +
	xlab('Evolutionary time (branch length)') + 
	ylab('Gene expression divergence (1-p)') + 
  	theme_bw()
dev.off()

data <- melt(Data[,2:3], measure.vars=c('oj', 'pj'))

svg(filename="gene_expression_divergence_tissues.svg")
ggplot(data) +
	geom_boxplot(aes(x=variable, y=value, color=variable)) +
	xlab('Tissue') + 
	ylab('Gene expression divergence across species (1-p)') + 
  	theme_bw()
dev.off()



#################################################################################################  mean transcript expression change over evolutionary time


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
#phylogeny <- read.delim("st26_phylogeny_ouch.txt")

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
#write.table(data_mean,"transcript_tpm_all_samples_onil_OJ_mean.txt", quote=F, sep="\t")

##normalize data##
data_mean.noNA = data_mean
data_mean.noNA[is.na(data_mean)] = 0
scale = calcNormFactors(data_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]
data_mean.norm = sapply(seq(1, length(scale)), function(x) data_mean[,x] / scale[x])
data_mean.norm = log10(data_mean.norm+0.01)
colnames(data_mean.norm) = colnames(data_mean)
print(head(data_mean.norm))
#write.table(data_mean.norm,"transcript_tpm_all_samples_onil_OJ_mean_norm.txt", quote=F, sep="\t")

# calculate isoform exp distance between samples

dist <- cor(data_mean.norm, method="spearman", use="complete.obs")
dist1 <- 1-dist
exp_dist_on <- dist1[,1]
exp_dist_on <- data.frame(exp_dist_on)
#exp_dist_on  <- exp_dist_on[order(row.names(exp_dist_on)), ]

p <- phylogeny[46:66,]
p1 <- p[,c(2,4)]
row.names(p1) <- p1$species
p1 <- p1[order(row.names(p1)), ]

out <- merge(p1, exp_dist_on, by="row.names")
out1 <- out[order(out$time), ]
out_oj <- out1



#######PJ
exp <- pj1
index <- read.table("gene_tpm_all_samples_onil_PJevee.index.txt", header=F, comment.char="")
#phylogeny <- read.delim("st26_phylogeny_ouch.txt")


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
#write.table(data_mean,"transcript_tpm_all_samples_onil_PJ_mean.txt", quote=F, sep="\t")

##normalize data##
data_mean.noNA = data_mean
data_mean.noNA[is.na(data_mean)] = 0
scale = calcNormFactors(data_mean.noNA, refColumn = 1, method="TMM")
scale = scale / scale[1]
data_mean.norm = sapply(seq(1, length(scale)), function(x) data_mean[,x] / scale[x])
data_mean.norm = log10(data_mean.norm+0.01)
colnames(data_mean.norm) = colnames(data_mean)
print(head(data_mean.norm))
#write.table(data_mean.norm,"transcript_tpm_all_samples_onil_PJ_mean_norm.txt", quote=F, sep="\t")

# calculate isoform exp distance between samples

dist <- cor(data_mean.norm, method="spearman", use="complete.obs")
dist1 <- 1-dist
exp_dist_on <- dist1[,1]
exp_dist_on <- data.frame(exp_dist_on)
#exp_dist_on  <- exp_dist_on[order(row.names(exp_dist_on)), ]

p <- phylogeny[46:66,]
p1 <- p[,c(2,4)]
row.names(p1) <- p1$species
p1 <- p1[order(row.names(p1)), ]


out <- merge(p1, exp_dist_on, by="row.names")
out1 <- out[order(out$time), ]
out_pj <- out1

## plot 1-p gene expression distance vs time (branch length)
library(ggplot2)
out3 <- merge(out_oj, out_pj, by="species")

Data1 = data.frame(out3$time.x, out3$exp_dist_on.x, out3$exp_dist_on.y)
colnames(Data1) <- c("time","oj","pj")

svg(filename="transcript_expression_divergence_over_time.svg")
ggplot(Data1) +
	geom_smooth(aes(Data1$time,oj),colour="black",fill="black") +
	geom_smooth(aes(Data1$time,pj),colour="grey",fill="grey") +
	xlab('Evolutionary time (branch length)') + 
	ylab('Gene expression divergence (1-p)') + 
  	theme_bw()
dev.off()

data1 <- melt(Data1[,2:3], measure.vars=c('oj', 'pj'))

svg(filename="transcript_expression_divergence_tissues.svg")
ggplot(data1) +
	geom_boxplot(aes(x=variable, y=value, color=variable)) +
	xlab('Tissue') + 
	ylab('Gene expression divergence across species (1-p)') + 
  	theme_bw()
dev.off()



##########################plot together

colnames(Data1) <- c("evoltime","t_oj","t_pj")
colnames(Data) <- c("evoltime","g_oj","g_pj")
both <- merge(Data, Data1, by="evoltime")


svg(filename="both_expression_divergence_over_time_v2.svg")
ggplot(both) +
	geom_smooth(aes(evoltime,g_oj),colour="hotpink1",fill="hotpink1") +
	geom_smooth(aes(evoltime,g_pj),colour="darkgrey",fill="darkgrey") +
	geom_smooth(aes(evoltime,t_oj),colour="deeppink4",fill="deeppink4") +
	geom_smooth(aes(evoltime,t_pj), colour="lightgrey",fill="lightgrey") +
	xlab('Evolutionary time (scaled branch length)') + 
	ylab('Expression divergence (1-p)') + 
	ylim(0.2,0.5)+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  	theme_bw()
dev.off()

svg(filename="phylogeny.svg")
plot(tree)
add.scale.bar()
dev.off()

both1 <- melt(both[,2:5])

svg(filename="both_expression_divergence_tissues.svg")
ggplot(both1) +
	geom_boxplot(aes(x=variable, y=value, color=variable)) +
	xlab('Tissue') + 
	ylab('Expression divergence across species (1-p)') + 
	ylim(0,0.6)+
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  	theme_bw()
dev.off()

print(t.test(both$g_pj, both$t_pj))
print(t.test(both$g_oj, both$t_oj))
print(t.test(both$g_pj, both$g_oj))
print(t.test(both$t_pj, both$t_oj))


### plot variance under OU model 
