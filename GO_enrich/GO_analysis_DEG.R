## ----1. Variables load or remove unneeded-----------------------------------------------
library("GSEABase")
library("Category")
library("topGO")
library("biomaRt")
library("goseq")
library(dplyr)
library(tidyverse)
library(Rgraphviz)
library(rtracklayer)

citation("topGo")
citation("biomaRt")
citation("topGO")


supportedOrganisms()


#[topGO](https://bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.pdf)

FILE.DIR        <- "/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/"
DATE            <- format(Sys.time(), "%m%d%Y")

FILE.DIR.BASE        <- "/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/DEG/"
PATH.OUT        <- FILE.DIR.BASE

# set directory 
setwd(FILE.DIR.BASE)
getwd()


#########################################################################
###  2. bioMart to get gene names and GO term
#########################################################################

# You only need to do this once 

# 16        oniloticus_gene_ensembl      Oreochromis niloticus  ENSONIG
 listMarts()
 ensembl <- useMart("ensembl")                   
 datasets <- listDatasets(ensembl)
 head(datasets)                   
 grep("G*", datasets)
 searchDatasets(mart = ensembl, pattern = "aculeatus")
 oniloticus_gene_ensembl
 listAttributes(ensembl)
 searchAttributes(mart=ensembl, pattern = "evidence")
 search for filter patterns 
 searchFilters(mart = ensembl, pattern = "gene")

myMart = useMart("ensembl", dataset="oniloticus_gene_ensembl", host = "https://www.ensembl.org")

## retrieve GO annotation from ENSEMBL
annot<-getBM(attributes = c("go_id","go_linkage_type","ensembl_gene_id"),
              filters = 'biotype',values = "protein_coding",    # only focus on protein coding genes
              mart = myMart)

annot<-getBM(attributes = c("go_id","go_linkage_type","ensembl_gene_id", "ensembl_transcript_id"),
             filters = 'biotype',values = "protein_coding",    # only focus on protein coding genes
             mart = myMart)



## retrieve GO annotations from ENSEMBL
annot<-getBM(attributes = c("go_id","go_linkage_type","ensembl_gene_id"),    # I am not sure if I should filter for something
             mart = myMart)

goframeData <- annot
goframeData <- goframeData[which(goframeData[,2]!=""),]  # remove missing values for go_linkage_type
goframeData <- goframeData[which(goframeData[,3]!="NA"),]  # remove NA in geneID
head(goframeData)

# select GOID and ensemblID
GOdata <- goframeData %>% 
  dplyr::select(go_id, ensembl_gene_id)

# save this
write_tsv(GOdata, "goframeData.txt", na = "NA", col_names = FALSE)

###########################################################################
### GO  Enrichment analysis 
###########################################################################

# build your GO2geneID or geneID2GO dataframe

GO2geneID <-readMappings(file = paste0(PATH.OUT, "goframeData.txt"), sep = "")
geneID2GO <- inverseList(GO2geneID)
geneNames <- names(geneID2GO)
str(head(geneID2GO))
length(geneNames)

#import ensemble gene universe and reduce to rna experiment universe

gene.universe.exp <-read.table("gene_universe.txt",header=F)
dim(gene.universe.exp)
sub <- as.character(gene.universe.exp$V1) %in% geneNames
geneNames.universe <- na.omit(geneNames[sub])
length(geneNames.universe)



######## Biological process enrichment 
###################################################################################### DEGSs

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/')

files = list.files(path=".", pattern="ourid.ensembl")
files

for(i in files){
  
# import significant DEGs
sig_genes <-read.table(i,header=F)
print(dim(sig_genes))
head(sig_genes)
sig_genes <- na.omit(sig_genes)
genesOfInterest <- as.character(sig_genes$V1)

# only my universe from the RNAseq experiment
geneList <- factor(as.integer(geneNames.universe %in% as.character(sig_genes$V1)))
names(geneList) <- geneNames.universe
str(geneList) 

######## Biological process enrichment 

GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO) 

# have a look at your topGOdata object: available genes (all genes in your geneList), feasable genes (genes that can be used for further analysis)
GOdata_BP
numGenes(GOdata_BP)
sg <- sigGenes(GOdata_BP)
str(sg)

# fisher classic

myGOdata <- GOdata_BP

# run the Fisher's exact tests
resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")

# see how many results we get where weight01 gives a P-value <= 0.05:
#mysummary <- summary(attributes(resultTopgo)$score <= 0.01)
#numsignif <- as.integer(mysummary[[3]])
#dim(numsignif)



# print out the top 'numsignif' results:
allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", numChar=1000)
allRes$topgoFisher.padjust <- p.adjust(allRes$topgoFisher, method="fdr")

write.table(allRes[allRes$topgoFisher.padjust < 0.05,], file = paste(i, "topgo_fisher_weight_tests_BP_q0.05.txt", sep="_"), quote = FALSE, row.names = FALSE, sep = "\t")


### genes associated with go enriched terms


## read in annotation

annos <- readGFF("O_niloticus_UMD_NMBU.99.gff3.genes")
anno <- annos[,c(10,12,13)]
head(anno)

my_list <- list()

myterms = (allRes[allRes$topgoFisher.padjust < 0.05,])$GO.ID

mygenes <- genesInTerm(myGOdata, myterms)
for (j in 1:length(myterms))
{
  myterm <- myterms[j]
  mygenesforterm <- mygenes[myterm][[1]]
  myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
  mygenesforterm2 <- mygenesforterm[myfactor == TRUE] 
  mygenesforterm3 <- data.frame(mygenesforterm2)
  mygenesforterm3$GO_id <- myterms[j]
  colnames(mygenesforterm3) <- c("gene_id", "GO_id")
  full <- merge(mygenesforterm3, anno, by="gene_id")
  my_list[[j]] <- full
}

# Write out single data frames
lst2 <- lapply(my_list,function(x) cbind(rowname=rownames(x),x))
df1 <- Reduce(function(x,y) merge(x,y,all=T),lst2)
out <- df1[,c(2:5)]
out <- out[order(out$GO_id, out$gene_id),]

details = (allRes[allRes$topgoFisher.padjust < 0.05,])[,c(1,2)]
details

out1 <- merge(out, details, by.x="GO_id", by.y="GO.ID")

write.table(out1, file = paste(i, "topgo_fisher_weight_tests_BP_q0.05sig_genes.txt", sep="_"), quote = FALSE, row.names = FALSE, sep = "\t")

}


########################################### convergent genes GO enrichment BP
###################################################################################### 


setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/')

files = list.files(path=".", pattern="_convergent_genes_up_down.txt")
files

for(i in files){
  
  # import significant DEGs
  sig_genes <-read.table(i,header=T)
  print(dim(sig_genes))
  head(sig_genes)
  sig_genes1 <- data.frame(do.call('rbind', strsplit(as.character(sig_genes$gene_id),'|',fixed=TRUE)))[2]
  genesOfInterest <- as.character(sig_genes1$X2)
  
  # only my universe from the RNAseq experiment
  geneList <- factor(as.integer(geneNames.universe %in%  as.character(sig_genes1$X2)))
  names(geneList) <- geneNames.universe
  str(geneList) 
  
  ######## Biological process enrichment 
  
  GOdata_BP <- new("topGOdata", ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO) 
  
  # have a look at your topGOdata object: available genes (all genes in your geneList), feasable genes (genes that can be used for further analysis)
  GOdata_BP
  numGenes(GOdata_BP)
  sg <- sigGenes(GOdata_BP)
  str(sg)
  
  # fisher classic
  
  myGOdata <- GOdata_BP
  
  # run the Fisher's exact tests
  resultClassic <- runTest(myGOdata, algorithm="classic", statistic="fisher")
  resultElim <- runTest(myGOdata, algorithm="elim", statistic="fisher")
  resultTopgo <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
  resultParentchild <- runTest(myGOdata, algorithm="parentchild", statistic="fisher")
  
  # see how many results we get where weight01 gives a P-value <= 0.05:
  #mysummary <- summary(attributes(resultTopgo)$score <= 0.01)
  #numsignif <- as.integer(mysummary[[3]])
  #dim(numsignif)
  
  
  
  # print out the top 'numsignif' results:
  allRes <- GenTable(myGOdata, classicFisher = resultClassic, elimFisher = resultElim, topgoFisher = resultTopgo, parentchildFisher = resultParentchild, orderBy = "topgoFisher", ranksOf = "classicFisher", numChar=1000)
  allRes$topgoFisher.padjust <- p.adjust(allRes$topgoFisher, method="fdr")
  
  write.table(allRes[allRes$topgoFisher.padjust < 0.05,], file = paste(i, "topgo_fisher_weight_tests_BP_q0.05.txt", sep="_"), quote = FALSE, row.names = FALSE, sep = "\t")
  


  ### genes associated with go enriched terms
  
  
  ## read in annotation
  
  annos <- readGFF("O_niloticus_UMD_NMBU.99.gff3.genes")
  anno <- annos[,c(10,12,13)]
  head(anno)
  
  sig_genes$genes <- sig_genes1$X2
  big <- merge(sig_genes, anno, by.x="genes", by.y="gene_id")
  write.table(big, file = paste(i, "gene_name.txt", sep="_"), quote = FALSE, row.names = FALSE, sep = "\t")
  
  
  my_list <- list()
  
  myterms = (allRes[allRes$topgoFisher.padjust < 0.05,])$GO.ID
  
  mygenes <- genesInTerm(myGOdata, myterms)
  for (j in 1:length(myterms))
  {
    myterm <- myterms[j]
    mygenesforterm <- mygenes[myterm][[1]]
    myfactor <- mygenesforterm %in% genesOfInterest # find the genes that are in the list of genes of interest
    mygenesforterm2 <- mygenesforterm[myfactor == TRUE] 
    mygenesforterm3 <- data.frame(mygenesforterm2)
    mygenesforterm3$GO_id <- myterms[j]
    colnames(mygenesforterm3) <- c("gene_id", "GO_id")
    full <- merge(mygenesforterm3, anno, by="gene_id")
    my_list[[j]] <- full
  }
  
  # Write out single data frames
  lst2 <- lapply(my_list,function(x) cbind(rowname=rownames(x),x))
  df1 <- Reduce(function(x,y) merge(x,y,all=T),lst2)
  out <- df1[,c(2:5)]
  out <- out[order(out$GO_id, out$gene_id),]
  
  details = (allRes[allRes$topgoFisher.padjust < 0.05,])[,c(1,2)]
  details
  
  out1 <- merge(out, details, by.x="GO_id", by.y="GO.ID")
  
  write.table(out1, file = paste(i, "topgo_fisher_weight_tests_BP_q0.05sig_genes.txt", sep="_"), quote = FALSE, row.names = FALSE, sep = "\t")
  
}


