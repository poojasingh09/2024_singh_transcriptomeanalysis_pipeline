

## libraries
library(stringr)

## wd

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/')

### species names
species <- rep(c("Aa", "Ab", "Ah" , "Ch" , "Gh" , "Gp" , "Ht" , "Lt"  ,"Ml" , "No" , "Pf" , "Pp" , "Prp" ,"Ps" , "Ptb", 
                  "Py" , "Sd"  ,"Sf" , "Tm" , "Tt"), each=5)

#species <- species[-58] # drop one pp

a <- read.table("gene_count_matrix_prepDE.csv", header=T, sep=",", row.names=1)
dim(a)

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

write.table(pj_out, "gene_count_matrix_prepDE_speciesmean_pj.csv", quote=F, sep="\t")
write.table(oj_out, "gene_count_matrix_prepDE_speciesmean_oj.csv", quote=F, sep="\t")

### species means plots

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/')
oj <- read.table("gene_count_matrix_prepDE_speciesmean_oj.csv", header=T)

## fzd2

fzd2 <- oj[grep(pattern = "MSTRG.34478", x = rownames(oj)),][1,]
a <- t(fzd2)
colnames(a) <- "fzd2"

# Barplot
ggplot(a, aes(x=fzd2, y=rownames(a))) + 
  geom_bar(stat = "identity") +
  coord_flip()



