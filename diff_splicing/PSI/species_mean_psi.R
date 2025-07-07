

## libraries
library(stringr)

## wd

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/suppa')

### species names
species <- rep(c("Aa", "Ab", "Ah" , "Ch" , "Gh" , "Gp" , "Ht" , "Lt"  ,"Ml" , "No" , "Pf" , "Pp" , "Prp" ,"Ps" , "Ptb", 
                  "Py" , "Sd"  ,"Sf" , "Tm" , "Tt"), each=5)

#species <- species[-58] # drop one pp

a <- read.table("filtered.gffcmp.all.annotated_isoform.psi", header=T, sep="\t", row.names=1)
dim(a)
head(a)

oj <- a[,grepl("OJ", colnames(a))]
dim(oj)
colnames(oj) <- species

pj <- a[,grepl("PJ", colnames(a))]
dim(pj)
colnames(pj) <- species


# calculate mean exp of biological replicates and manipulate gene name
oj_out <- t(apply(oj, 1, function(x) tapply(x, colnames(oj), mean)))
head(oj)
head(oj_out)

pj_out <- t(apply(pj, 1, function(x) tapply(x, colnames(pj), mean)))
head(pj)
head(pj_out)

oj_out_t <- t(oj_out)
pj_out_t <- t(pj_out)

head(oj_out)
head(pj_out)

write.table(pj_out, "filtered.gffcmp.all.annotated_isoform.psi_speciesmean_pj.csv", quote=F, sep="\t")
write.table(oj_out, "iltered.gffcmp.all.annotated_isoform.psi_speciesmean_oj.csv", quote=F, sep="\t")


