#Pooja Singh 
#march 2023
#EVEE tools, selection on gene expression https://evee-tools.github.io

library("stringr") 

##########

setwd("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/RNAseq_newgrant/allst26/analysis/new_ref/05_selection_geneexpression/selection/EVEEtools")



################# GENE EXPRESSION


#read in data
countData <- read.csv("gene_tpm_all_samples_onil.tsv", sep="\t")

## after analysing PCAs, drop outlier samples"
drop <- c("Pp.L.OJ.3")
countData <- countData[ , !(colnames(countData) %in% drop)]


##subset
# add orthogroup, i dont have orthogrouops as i used the same ref, so i just code each each as a unique OG

countData1 <- countData[,grepl("OJ", colnames(countData))]
On <- countData[,grepl("On20d3", colnames(countData))]
GENE_NAME <- countData[,grepl("Gene_ID", colnames(countData))]
countData1a <- cbind(On, countData1)
countData1b <- cbind(GENE_NAME, countData1a)
rownames(countData1b) <- paste0("OG", rownames(countData1b))

countData2 <- countData[,grepl("PJ", colnames(countData))]
countData2a <- cbind(On, countData2)
countData2b <- cbind(GENE_NAME, countData2a)
rownames(countData2b) <- paste0("OG", rownames(countData2b))

head(countData2b)
head(countData1b)

#rename header as eveetools wants it

#header <- str_remove_all(colnames(countData1a), ".L.*")
#colnames(countData1a) <- header

#header <- str_remove_all(colnames(countData2a), ".L.*")
#colnames(countData2a) <- header


write.table(countData1b, "gene_tpm_all_samples_onil_OJevee.tsv", quote=F, row.names=T, sep="\t")
write.table(countData2b, "gene_tpm_all_samples_onil_PJevee.tsv", quote=F, row.names=T, sep="\t")

## build index needed for eveee

c1 <- colnames(countData1b)
c2 <- str_remove_all(c1, ".L.*")
c3 <- cbind(data.frame(c1),data.frame(c2))
write.table(c3, "gene_tpm_all_samples_onil_OJevee.index.txt", quote=F, row.names=F, sep="\t")


c1 <- colnames(countData2b)
c2 <- str_remove_all(c1, ".L.*")
c3 <- cbind(data.frame(c1),data.frame(c2))
write.table(c3, "gene_tpm_all_samples_onil_PJevee.index.txt", quote=F, row.names=F, sep="\t")

## after the code remove the top 2 rows of the indices



################# ISOFORM EXPRESSION


#read in data
countData <- read.csv("transcript_tpm_all_samples_onil.tsv", sep="\t")

## after analysing PCAs, drop outlier samples"
drop <- c("Pp.L.OJ.3")
countData <- countData[ , !(colnames(countData) %in% drop)]


##subset
# add orthogroup, i dont have orthogrouops as i used the same ref, so i just code each each as a unique OG

countData1 <- countData[,grepl("OJ", colnames(countData))]
On <- countData[,grepl("On20d3", colnames(countData))]
transcript_NAME <- countData[,grepl("Transcript_ID", colnames(countData))]
countData1a <- cbind(On, countData1)
countData1b <- cbind(transcript_NAME, countData1a)
rownames(countData1b) <- paste0("OG", rownames(countData1b))

countData2 <- countData[,grepl("PJ", colnames(countData))]
countData2a <- cbind(On, countData2)
countData2b <- cbind(transcript_NAME, countData2a)
rownames(countData2b) <- paste0("OG", rownames(countData2b))

head(countData2b)
head(countData1b)

#rename header as eveetools wants it

#header <- str_remove_all(colnames(countData1a), ".L.*")
#colnames(countData1a) <- header

#header <- str_remove_all(colnames(countData2a), ".L.*")
#colnames(countData2a) <- header


write.table(countData1b, "transcript_tpm_all_samples_onil_OJevee.tsv", quote=F, row.names=T, sep="\t")
write.table(countData2b, "transcript_tpm_all_samples_onil_PJevee.tsv", quote=F, row.names=T, sep="\t")

## build index needed for eveee

c1 <- colnames(countData1b)
c2 <- str_remove_all(c1, ".L.*")
c3 <- cbind(data.frame(c1),data.frame(c2))
write.table(c3, "transcript_tpm_all_samples_onil_OJevee.index.txt", quote=F, row.names=F, sep="\t")


c1 <- colnames(countData2b)
c2 <- str_remove_all(c1, ".L.*")
c3 <- cbind(data.frame(c1),data.frame(c2))
write.table(c3, "transcript_tpm_all_samples_onil_PJevee.index.txt", quote=F, row.names=F, sep="\t")

## after the code remove the top row of the indices and switch col1 and col2


################### PSI

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/selection/EVEEtools/')
#read in data
countData <- read.csv("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/suppa/with_onil/filtered.gffcmp.all.annotated_onil_isoform.psi", sep="\t")

## after analysing PCAs, drop outlier samples"
drop <- c("Pp.L.OJ.3")
countData <- countData[ , !(colnames(countData) %in% drop)]


##subset
# add orthogroup, i dont have orthogrouops as i used the same ref, so i just code each each as a unique OG

countData1 <- countData[,grepl("OJ", colnames(countData))]
On <- countData[,grepl("On20d3", colnames(countData))]
transcript_NAME <- countData[,grepl("Transcript_ID", colnames(countData))]
countData1a <- cbind(On, countData1)
countData1b <- cbind(transcript_NAME, countData1a)
rownames(countData1b) <- paste0("OG", rownames(countData1b))

countData2 <- countData[,grepl("PJ", colnames(countData))]
countData2a <- cbind(On, countData2)
countData2b <- cbind(transcript_NAME, countData2a)
rownames(countData2b) <- paste0("OG", rownames(countData2b))

head(countData2b)
head(countData1b)

#rename header as eveetools wants it

#header <- str_remove_all(colnames(countData1a), ".L.*")
#colnames(countData1a) <- header

#header <- str_remove_all(colnames(countData2a), ".L.*")
#colnames(countData2a) <- header



write.table(countData1b, "isoform_psi_all_samples_onil_OJevee.tsv", quote=F, row.names=T, sep="\t")
write.table(countData2b, "isoform_psi_all_samples_onil_PJevee.tsv", quote=F, row.names=T, sep="\t")

## build index needed for eveee

c1 <- colnames(countData1b)
c2 <- str_remove_all(c1, ".L.*")
c3 <- cbind(data.frame(c1),data.frame(c2))
write.table(c3, "isoform_psi_all_samples_onil_OJevee.index.txt", quote=F, row.names=F, sep="\t")


c1 <- colnames(countData2b)
c2 <- str_remove_all(c1, ".L.*")
c3 <- cbind(data.frame(c1),data.frame(c2))
write.table(c3, "isoform_psi_all_samples_onil_PJevee.index.txt", quote=F, row.names=F, sep="\t")

## after the code remove the top row of the indices and switch col1 and col2


