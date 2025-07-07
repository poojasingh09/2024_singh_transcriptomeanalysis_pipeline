
## upset plots

# load libraries
library(plyr)
library(SuperExactTest)

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/suppa')
## OJ

lt <- read.table("oj_lt.lakes.dpsi.sig.genes", header=F)
lm <- read.table("oj_lm.lakes.dpsi.sig.genes", header=F)
lv <- read.table("oj_lv.lakes.dpsi.sig.genes", header=F)

lt1 <- unique(data.frame(do.call('rbind', strsplit(as.character(lt$V1),';',fixed=TRUE)))[1])
lm1 <- unique(data.frame(do.call('rbind', strsplit(as.character(lm$V1),';',fixed=TRUE)))[1])
lv1 <- unique(data.frame(do.call('rbind', strsplit(as.character(lv$V1),';',fixed=TRUE)))[1])

list <- list(lt1$X1, lm1$X1, lv1$X1)

(length.gene.sets=sapply(list,length))

total = 36296


res=supertest(list, n=total)
res$set.names <- c("LT OJ", "LM OJ", "LV OJ")

svg("OJ_hypergeometric_overlap_2025_v1.svg")
plot(res, Layout="landscape", sort.by="size", margin=c(0.5,5,1,2), ylab="Differentially Spliced Genes \nherbivores vs carnivores")
dev.off()

ltgenes <- lt1$X1
lmgenes <- lm1$X1
lvgenes <- lv1$X1

shared <- data.frame(intersect(intersect(ltgenes,lmgenes),lvgenes))
write.table(shared, "OJ_convergent_isoforms.txt", quote=F, row.names=F, sep="\t")

## PJ


lt <- read.table("pj_lt.lakes.dpsi.sig.genes", header=F)
lm <- read.table("pj_lm.lakes.dpsi.sig.genes", header=F)
lv <- read.table("pj_lv.lakes.dpsi.sig.genes", header=F)


lt1 <- unique(data.frame(do.call('rbind', strsplit(as.character(lt$V1),';',fixed=TRUE)))[1])
lm1 <- unique(data.frame(do.call('rbind', strsplit(as.character(lm$V1),';',fixed=TRUE)))[1])
lv1 <- unique(data.frame(do.call('rbind', strsplit(as.character(lv$V1),';',fixed=TRUE)))[1])

list <- list(lt1$X1, lm1$X1, lv1$X1)

(length.gene.sets=sapply(list,length))

total = 36296


res=supertest(list, n=total)
res$set.names <- c("LT PJ", "LM PJ", "LV PJ")

svg("PJ_hypergeometric_overlap_2025_v1.svg")
plot(res, Layout="landscape", sort.by="size", margin=c(0.5,5,1,2), ylab="Differentially Spliced Genes \nherbivores vs carnivores")
dev.off()


ltgenes <- lt1$X1
lmgenes <- lm1$X1
lvgenes <- lv1$X1

shared <- data.frame(intersect(intersect(ltgenes,lmgenes),lvgenes))
write.table(shared, "PJ_convergent_PSI.txt", quote=F, row.names=F, sep="\t")
