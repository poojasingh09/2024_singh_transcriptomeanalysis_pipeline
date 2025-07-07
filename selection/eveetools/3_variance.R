# pooja.singh
# nov 2023
# here i am wrangling the output of EVEE tools to identify the genes and transcripts that fit the BM or OU model


# libraries
library(ggplot2)


setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/selection/EVEEtools')

## read in data
oj_gene <- read.table("OJ_ou_stats_gene.txt")
oj_transcript <- read.table("OJ_ou_stats_transcript.txt")

pj_gene <- read.table("PJ_ou_stats_gene.txt")
pj_transcript <- read.table("PJ_ou_stats_transcript.txt")

head(oj_gene)


################ neutral vs selective

# lets first set up a data frame to store the data
out <- data.frame(matrix(nrow = 2, ncol = 4))
colnames(out) <- c("OJ_gene","OJ_transcript","PJ_gene","PJ_transcript")
rownames(out) <- c("neutral","selective") 

## brownian motion rejected (p < 0.05) and OU model significant

oj_gene_sig <- oj_gene[oj_gene$qvalues < 0.05,]
oj_gene_nonsig <- oj_gene[oj_gene$qvalues >= 0.05,]
out[1,1] <- nrow(oj_gene_nonsig) 
out[2,1] <- nrow(oj_gene_sig) 

oj_transcript_sig <- oj_transcript[oj_transcript$qvalues < 0.05,]
oj_transcript_nonsig <- oj_transcript[oj_transcript$qvalues >= 0.05,]
out[1,2] <- nrow(oj_transcript_nonsig) 
out[2,2] <- nrow(oj_transcript_sig) 


pj_gene_sig <- pj_gene[pj_gene$qvalues < 0.05,]
pj_gene_nonsig <- pj_gene[pj_gene$qvalues >= 0.05,]
out[1,3] <- nrow(pj_gene_nonsig) 
out[2,3] <- nrow(pj_gene_sig) 

pj_transcript_sig <- pj_transcript[pj_transcript$qvalues < 0.05,]
pj_transcript_nonsig <- pj_transcript[pj_transcript$qvalues >= 0.05,]
out[1,4] <- nrow(pj_transcript_nonsig) 
out[2,4] <- nrow(pj_transcript_sig) 



# lets first set up a data frame to store the data for stacked plots
out1 <- data.frame(matrix(nrow = 8, ncol = 3))
colnames(out1) <- c("category","neutral","number")
out1$category <- c("OJ_g","OJ_g","OJ_t", "OJ_t", "PJ_g","PJ_g","PJ_t", "PJ_t")
out1$neutral <- c("neutral","selective","neutral","selective", "neutral","selective", "neutral","selective")

out1[1,3] <- nrow(oj_gene_nonsig) 
out1[2,3] <- nrow(oj_gene_sig) 
out1[3,3] <- nrow(oj_transcript_nonsig) 
out1[4,3] <- nrow(oj_transcript_sig) 
out1[5,3] <- nrow(pj_gene_nonsig) 
out1[6,3] <- nrow(pj_gene_sig) 
out1[7,3] <- nrow(pj_transcript_nonsig) 
out1[8,3] <- nrow(pj_transcript_sig) 

# plot stacked
svg("neutral_selective.svg")
ggplot(out1, aes(fill=neutral, y=number, x=category)) + 
  geom_bar(position="stack", stat="identity") + coord_flip()
dev.off()

# unstacked
out1 <- data.frame(matrix(nrow = 8, ncol = 2))
#BM
out1[1,2] <- nrow(oj_gene_nonsig)/nrow(oj_gene) 
out1[2,2] <- nrow(pj_gene_nonsig)/nrow(pj_gene)
out1[3,2] <- nrow(oj_transcript_nonsig)/nrow(oj_transcript)  
out1[4,2] <- nrow(pj_transcript_nonsig)/nrow(pj_transcript) 

#OU
out1[5,2] <- nrow(oj_gene_sig)/nrow(oj_gene) 
out1[6,2] <- nrow(pj_gene_sig)/nrow(pj_gene) 
out1[7,2] <- nrow(oj_transcript_sig)/nrow(oj_transcript) 
out1[8,2] <- nrow(pj_transcript_sig)/nrow(pj_transcript) 

##categories
out1[1,1] <- "oj_gene_bm"
out1[2,1] <- "pj_gene_bm"
out1[3,1] <- "oj_transcript_bm"
out1[4,1] <- "pj_transcript_bm"
out1[5,1] <- "oj_gene_ou"
out1[6,1] <- "pj_gene_ou"
out1[7,1] <- "oj_transcript_ou" 
out1[8,1] <- "pj_transcript_ou"

colnames(out1) <- c("category","number")

out1$category <- factor(out1$category, levels = out1$category)

svg("figure_s1_neutral_selective.svg")
ggplot(out1, aes(y=number, x=category)) + 
  geom_bar(stat = "identity", width=0.4, ) + coord_flip() +
  theme_bw() + theme(aspect.ratio = 1/4) + ylab("proportion") + xlab("")
dev.off()

out2 <- out1[c(1,2,5,6),]

svg("figure_3_neutral_selective.svg")
ggplot(out2, aes(y=number, x=category)) + 
  geom_bar(stat = "identity", width=0.4, ) + coord_flip() +
  theme_bw() + theme(aspect.ratio = 1/4) + ylab("proportion") + xlab("")
dev.off()

############### variance all

oj_gene <- oj_gene[oj_gene$qvalues < 0.05,]
oj_transcript <- oj_transcript[oj_transcript$qvalues < 0.05,]
pj_gene <- pj_gene[pj_gene$qvalues < 0.05,]
pj_transcript <- pj_transcript[pj_transcript$qvalues < 0.05,]


svg("figure3_evolutionary_variance_oj_pj.svg")
boxplot(oj_gene$var, pj_gene$var, oj_transcript$var, pj_transcript$var, outline=F, ylab="evolutionary variance", ylim=c(0,1), col=c("hotpink1",  "lightgrey", "deeppink4","darkgrey"))
dev.off()   


t.test(oj_gene$var, oj_transcript$var) #p-value < 2.2e-16
t.test(pj_gene$var, pj_transcript$var) #p-value < 2.2e-16
t.test(oj_gene$var, pj_gene$var)  #p-value = 0.6477
t.test(oj_transcript$var, pj_transcript$var) #p-value = 0.8289


############### variance OJ vs PJ genes

lt <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/LT_PJOJ_DESeq_fullresults_0.05.ourid.only")
lm <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/LM_PJOJ_DESeq_fullresults_0.05.ourid.only")
lv <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/LV_PJOJ_DESeq_fullresults_0.05.ourid.only")
conv <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/PJOJ_convergent_genes_up_down.ourid")
alx <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/alx.ourid")

lt_oj_gene <- merge(oj_gene, lt, by.x="gene_name", by.y="V1")
lm_oj_gene <- merge(oj_gene, lm, by.x="gene_name", by.y="V1")
lv_oj_gene <- merge(oj_gene, lv, by.x="gene_name", by.y="V1")
conv_oj_gene <- merge(oj_gene, conv, by.x="gene_name", by.y="V1")
alx_oj_gene <- merge(oj_gene, alx, by.x="gene_name", by.y="V1")

lt_pj_gene <- merge(pj_gene, lt, by.x="gene_name", by.y="V1")
lm_pj_gene <- merge(pj_gene, lm, by.x="gene_name", by.y="V1")
lv_pj_gene <- merge(pj_gene, lv, by.x="gene_name", by.y="V1")
conv_pj_gene <- merge(pj_gene, conv, by.x="gene_name", by.y="V1")
alx_pj_gene <- merge(pj_gene, alx, by.x="gene_name", by.y="V1")



#svg("figure3_evolutionary_variance_oj_pj_DEGs.svg")
boxplot(conv_oj_gene$var, conv_pj_gene$var, alx_oj_gene$var, alx_pj_gene$var, outline=F, ylab="evolutionary variance", ylim=c(0,1), col=c("hotpink1",  "lightgrey", "deeppink4","darkgrey", "green"))
abline(h=median(na.omit(oj_gene$var)), lty=2)
#dev.off()   

#svg("figure3_evolutionary_variance_oj_pj_DEGs.svg")
boxplot(lt_pj_gene$var, lm_pj_gene$var, lv_pj_gene$var, conv_pj_gene$var, alx_pj_gene$var, outline=F, ylab="evolutionary variance", ylim=c(0,1), col=c("hotpink1",  "lightgrey", "deeppink4","darkgrey", "green"))
abline(h=median(na.omit(oj_gene$var)), lty=2)
#dev.off() 



countoj <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/gene_count_matrix_prepDE_speciesmean_oj.csv")
countoj$names <- rownames(countoj) 
countoj$names <-  sapply(strsplit(countoj$names, split="|", fixed = TRUE), `[`, 1)   
head(countoj)
alx_count <- countoj[countoj$names %in% alx$V1, ]

a <- t((alx_count[1,c(1:20)]))
barplot(a[,1])



t.test(oj_gene$var, oj_transcript$var) #p-value < 2.2e-16
t.test(pj_gene$var, pj_transcript$var) #p-value < 2.2e-16
t.test(oj_gene$var, pj_gene$var)  #p-value = 0.6477
t.test(oj_transcript$var, pj_transcript$var) #p-value = 0.8289



############### brownSigma: divergence under brownian motion model

boxplot(oj_gene$brownSigmas, oj_transcript$brownSigmas, pj_gene$brownSigmas, pj_transcript$brownSigmas, outline=F, ylab="evolutionary brownSigmasiance", col=c("red","blue","darkgrey","lightgrey"))

t.test(oj_gene$brownSigmas, oj_transcript$brownSigmas)
t.test(pj_gene$brownSigmas, pj_transcript$brownSigmas)
t.test(oj_gene$brownSigmas, pj_gene$brownSigmas)
t.test(oj_transcript$brownSigmas, pj_transcript$brownSigmas)

################# lists for GO enrichment
head(oj_gene)
dim(oj_gene)

g_oj <- quantile(na.omit(oj_gene_sig$var), prob=c(0.99, 0.01))
g_oj99 <- oj_gene_sig[na.omit(oj_gene_sig$var) > g_oj[1],]
g_oj1 <- oj_gene_sig[na.omit(oj_gene_sig$var) < g_oj[2],]
write.table(g_oj99, '/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/selection/oj_geneexp_OUsig_99variance.txt', quote=F, row.names=F)
write.table(g_oj1, '/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/selection/oj_geneexp_OUsig_1variance.txt', quote=F, row.names=F)
dim(g_oj1)
dim(g_oj99)

g_pj <- quantile(na.omit(pj_gene_sig$var), prob=c(0.99, 0.01))
g_pj99 <- pj_gene_sig[na.omit(pj_gene_sig$var) > g_oj[1],]
g_pj1 <- pj_gene_sig[na.omit(pj_gene_sig$var) < g_oj[2],]
write.table(g_pj99, '/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/selection/pj_geneexp_OUsig_99variance.txt', quote=F, row.names=F)
write.table(g_pj1, '/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/selection/pj_geneexp_OUsig_1variance.txt', quote=F, row.names=F)
dim(g_oj1)
dim(g_pj99)


g_oj <- quantile(na.omit(oj_transcript_sig$var), prob=c(0.99, 0.01))
g_oj99 <- oj_transcript_sig[na.omit(oj_transcript_sig$var) > g_oj[1],]
g_oj1 <- oj_transcript_sig[na.omit(oj_transcript_sig$var) < g_oj[2],]
write.table(g_oj99, '/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/selection/oj_transcriptexp_OUsig_99variance.txt', quote=F, row.names=F)
write.table(g_oj1, '/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/selection/oj_transcriptexp_OUsig_1variance.txt', quote=F, row.names=F)
dim(g_oj1)
dim(g_oj99)

g_pj <- quantile(na.omit(pj_transcript_sig$var), prob=c(0.99, 0.01))
g_pj99 <- pj_transcript_sig[na.omit(pj_transcript_sig$var) > g_oj[1],]
g_pj1 <- pj_transcript_sig[na.omit(pj_transcript_sig$var) < g_oj[2],]
write.table(g_pj99, '/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/selection/pj_transcriptexp_OUsig_99variance.txt', quote=F, row.names=F)
write.table(g_pj1, '/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/selection/pj_transcriptexp_OUsig_1variance.txt', quote=F, row.names=F)
dim(g_pj1)
dim(g_pj99)




getwd()

######## run GO enrichment here
#setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/selection/geneexp')



######## check variance of convergent genes

## read in data


oj_conv <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_convergent_genes_up_down.mstrg")
pj_conv <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/pj_convergent_genes_up_down.mstrg")

lv_oj <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/LV_OJ_herbcarn_DESeq_fullresults_0.05.ourid.mstrg")
lm_oj <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/LM_OJ_herbcarn_DESeq_fullresults_0.05.ourid.mstrg")
lt_oj <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/LT_OJ_herbcarn_DESeq_fullresults_0.05.ourid.mstrg")

splice <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/spliceosome_genes.ensembl", header=F) 
splice_ass <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/spliceosome_associated_genes.ensembl", header=F) 


lt_ojpj <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/LT_PJOJ_DESeq_fullresults_0.05.ourid", header=F) 
lm_ojpj <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/LM_PJOJ_DESeq_fullresults_0.05.ourid", header=F) 
lv_ojpj <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/LV_PJOJ_DESeq_fullresults_0.05.ourid", header=F) 
ojpj <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/PJOJ_convergent_genes_up_down.ourid", header=F) 
alx <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_vs_pj/alx.ourid", header=F) 



oj_conv1 <- merge(oj_conv, oj_gene, by.x="V1", by.y="gene_name")
pj_conv1 <- merge(pj_conv, pj_gene, by.x="V1", by.y="gene_name")
lv_oj1 <- merge(lv_oj, oj_gene, by.x="V1", by.y="gene_name")
lm_oj1 <- merge(lm_oj, oj_gene, by.x="V1", by.y="gene_name")
lt_oj1 <- merge(lt_oj, oj_gene, by.x="V1", by.y="gene_name")
splice1 <- merge(splice, oj_gene, by.x="V1", by.y="gene_name")
splice_ass1 <- merge(splice_ass, oj_gene, by.x="V1", by.y="gene_name")
ojpj1 <- merge(ojpj, oj_gene, by.x="V1", by.y="gene_name")
lt_ojpj1 <- merge(lt_ojpj, oj_gene, by.x="V1", by.y="gene_name")
lm_ojpj1 <- merge(lm_ojpj, oj_gene, by.x="V1", by.y="gene_name")
lv_ojpj1 <- merge(lv_ojpj, oj_gene, by.x="V1", by.y="gene_name")

alx1oj <- merge(alx, oj_gene, by.x="V1", by.y="gene_name")
alx1pj <- merge(alx, pj_gene, by.x="V1", by.y="gene_name")

boxplot(log10(oj_gene$var), log10(oj_transcript$var), log10(pj_gene$var), log10(pj_transcript$var), log10(oj_conv1$var), log10(pj_conv1$var), log10(lv_oj1$var), log10(lm_oj1$var),log10(lt_oj1$var),outline=T)

boxplot(oj_gene$var, oj_transcript$var, pj_gene$var, pj_transcript$var, oj_conv1$var, pj_conv1$var, lv_oj1$var, lm_oj1$var,lt_oj1$var,outline=F)

boxplot(log10(oj_gene$var), log10(pj_gene$var),log10(oj_conv1$var), log10(pj_conv1$var), log10(lv_oj1$var), log10(lm_oj1$var),log10(lt_oj1$var), log10(splice1$var), outline=T)

boxplot(oj_gene$var, pj_gene$var, oj_conv1$var, pj_conv1$var, lv_oj1$var, lm_oj1$var,lt_oj1$var, splice1$var, splice_ass1$var, outline=F)


boxplot(lt_ojpj1$var, lm_ojpj1$var, lv_ojpj1$var, ojpj1$var, alx1oj$var, outline=F)
wilcox.test(ojpj1$var, alx1oj$var)

