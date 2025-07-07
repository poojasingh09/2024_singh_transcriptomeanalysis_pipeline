
## upset plots for differential gene expression
## pooja.singh09@gmail.com

# load libraries
library(plyr)
library(SuperExactTest)

## set wd

setwd("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023")
## OJ

#### in shell

#cut -f 1 LT_OJ_herbcarn_DESeq_fullresults_0.05.txt | sed 1d | sort -k1n | uniq > LT_OJ_herbcarn_DESeq_fullresults_0.05.ourid
#cut -f 1 LM_OJ_herbcarn_DESeq_fullresults_0.05.txt | sed 1d | sort -k1n | uniq > LM_OJ_herbcarn_DESeq_fullresults_0.05.ourid
#cut -f 1 LV_OJ_herbcarn_DESeq_fullresults_0.05.txt | sed 1d | sort -k1n | uniq > LV_OJ_herbcarn_DESeq_fullresults_0.05.ourid

############# in R
lt <- read.table("LT_OJ_herbcarn_DESeq_fullresults_0.05.ourid", header=F)
lm <- read.table("LM_OJ_herbcarn_DESeq_fullresults_0.05.ourid", header=F)
lv <- read.table("LV_OJ_herbcarn_DESeq_fullresults_0.05.ourid", header=F)

mylist <- list(lt$V1, lm$V1, lv$V1)

#(length.gene.sets=sapply(list,length))

total = 38995


res=supertest(mylist, n=total)
res$set.names <- c("LT", "LM", "LV")

pdf("DEG_OJ_hypergeometric_overlap.pdf")
plot(res, Layout="landscape", sort.by="size", margin=c(0.5,5,1,2), ylab="Differentially expressed genes \n herbivores vs carnivores")
dev.off()

svg("DEG_OJ_hypergeometric_overlap.svg")
plot(res, Layout="landscape", sort.by="size", ylab="Differentially expressed genes \n herbivores vs carnivores")
dev.off()


ltgenes <- lt$V1
lmgenes <- lm$V1
lvgenes <- lv$V1

shared <- data.frame(intersect(intersect(ltgenes,lmgenes),lvgenes))
colnames(shared) <- "shared_ltlmlv"
write.table(shared, "OJ_convergent_genes.txt", quote=F, row.names=F, sep="\t")
dim(shared)

############################################  direction of convergence OJ

lte <- read.table("LT_OJ_herbcarn_DESeq_fullresults_0.05.txt", header=T)
lme <- read.table("LM_OJ_herbcarn_DESeq_fullresults_0.05.txt", header=T)
lve <- read.table("LV_OJ_herbcarn_DESeq_fullresults_0.05.txt", header=T)

lte1 <- lte[2]
lte1$id <- rownames(lte)
lme1 <- lme[2]
lme1$id <- rownames(lme)
lve1 <- lve[2]
lve1$id <- rownames(lve)

a <- merge(shared, lte1,by.x="shared_ltlmlv", by.y="id")
b <- merge(a, lme1, by.x="shared_ltlmlv", by.y="id", suffixes = c("_lt_oj","_lm_oj"))
c <- merge(b, lve1, by.x="shared_ltlmlv", by.y="id")
colnames(c) = c("gene_id", "lt_pj","lm_pj", "lv_pj")

c$down <- with(c,c$lt_pj < 0 & c$lm_pj < 0 & c$lv_pj < 0)
c$up <- with(c,c$lt_pj > 0 & c$lm_pj > 0 & c$lv_pj > 0)

c$convergent <- with(c,c$up == 'TRUE' | c$down == 'TRUE')
oj_convergent_genes_same_direction <- c[c$convergent == 'TRUE',]
write.table(oj_convergent_genes_same_direction, "oj_convergent_genes_up_down.txt", quote=F, row.names=F, sep="\t")
oj_convergent_genes_same_direction$ens <- sapply(strsplit(oj_convergent_genes_same_direction$gene_id, split="|", fixed = TRUE), `[`, 2)
write.table(oj_convergent_genes_same_direction$ens, "oj_convergent_genes_up_down.genes.txt", quote=F, row.names=F, sep="\t")

out <- file("oj_convergent_genes_up_down_stats.txt", "w")
writeLines((paste0("number of pj convergent genes upregulated ", nrow(c[c$up == 'TRUE',]))), out)
writeLines((paste0("number of pj convergent genes downregulated ", nrow(c[c$down == 'TRUE',]))), out)
writeLines((paste0("number of pj convergent genes total ", nrow(c[c$convergent == 'TRUE',]))), out)
close(out)


################################################################################################################################################################################## PJ


#### in shell

#cut -f 1 LT_PJ_herbcarn_DESeq_fullresults_0.05.txt | sed 1d | sort -k1n | uniq > LT_PJ_herbcarn_DESeq_fullresults_0.05.ourid
#cut -f 1 LM_PJ_herbcarn_DESeq_fullresults_0.05.txt | sed 1d | sort -k1n | uniq > LM_PJ_herbcarn_DESeq_fullresults_0.05.ourid
#cut -f 1 LV_PJ_herbcarn_DESeq_fullresults_0.05.txt | sed 1d | sort -k1n | uniq > LV_PJ_herbcarn_DESeq_fullresults_0.05.ourid


############# in R
lt <- read.table("LT_PJ_herbcarn_DESeq_fullresults_0.05.ourid", header=F)
lm <- read.table("LM_PJ_herbcarn_DESeq_fullresults_0.05.ourid", header=F)
lv <- read.table("LV_PJ_herbcarn_DESeq_fullresults_0.05.ourid", header=F)

list <- list(lt$V1, lm$V1, lv$V1)

(length.gene.sets=sapply(list,length))

total = 38995


res=supertest(list, n=total)
res$set.names <- c("LT", "LM", "LV")

pdf("DEG_PJ_hypergeometric_overlap.pdf")
plot(res, Layout="landscape", sort.by="size", margin=c(0.5,5,1,2), ylab="Differentially expressed genes \n herbivores vs carnivores")
dev.off()

svg("DEG_PJ_hypergeometric_overlap.svg")
plot(res, Layout="landscape", sort.by="size", ylab="Differentially expressed genes \n herbivores vs carnivores")
dev.off()


ltgenes <- lt$V1
lmgenes <- lm$V1
lvgenes <- lv$V1

shared <- data.frame(intersect(intersect(ltgenes,lmgenes),lvgenes))
colnames(shared) <- "shared_ltlmlv_pj"
write.table(shared, "PJ_convergent_genes.txt", quote=F, row.names=F, sep="\t")
dim(shared)


############################################  direction of convergence PJ

lte <- read.table("LT_PJ_herbcarn_DESeq_fullresults_0.05.txt", header=T)
lme <- read.table("LM_PJ_herbcarn_DESeq_fullresults_0.05.txt", header=T)
lve <- read.table("LV_PJ_herbcarn_DESeq_fullresults_0.05.txt", header=T)

lte1 <- lte[2]
lte1$id <- rownames(lte)
lme1 <- lme[2]
lme1$id <- rownames(lme)
lve1 <- lve[2]
lve1$id <- rownames(lve)

a <- merge(shared, lte1, by.x="shared_ltlmlv_pj", by.y="id")
b <- merge(a, lme1, by.x="shared_ltlmlv_pj", by.y="id", suffixes = c("_lt_pj","_lm_pj"))
c <- merge(b, lve1, by.x="shared_ltlmlv_pj", by.y="id")
colnames(c) = c("gene_id", "lt_pj","lm_pj", "lv_pj")

c$down <- with(c,c$lt_pj < 0 & c$lm_pj < 0 & c$lv_pj < 0)
c$up <- with(c,c$lt_pj > 0 & c$lm_pj > 0 & c$lv_pj > 0)

c$convergent <- with(c,c$up == 'TRUE' | c$down == 'TRUE')
pj_convergent_genes_same_direction <- c[c$convergent == 'TRUE',]
write.table(pj_convergent_genes_same_direction, "pj_convergent_genes_up_down.txt", quote=F, row.names=F, sep="\t")
pj_convergent_genes_same_direction$ens <- sapply(strsplit(pj_convergent_genes_same_direction$gene_id, split="|", fixed = TRUE), `[`, 2)
write.table(pj_convergent_genes_same_direction$ens, "pj_convergent_genes_up_down.genes.txt", quote=F, row.names=F, sep="\t")

out <- file("pj_convergent_genes_up_down_stats.txt", "w")
writeLines((paste0("number of pj convergent genes upregulated ", nrow(c[c$up == 'TRUE',]))), out)
writeLines((paste0("number of pj convergent genes downregulated ", nrow(c[c$down == 'TRUE',]))), out)
writeLines((paste0("number of pj convergent genes total ", nrow(c[c$convergent == 'TRUE',]))), out)
close(out)






