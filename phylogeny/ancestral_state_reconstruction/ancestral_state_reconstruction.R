################################################################################################# ancestral trait reconstuctions of gene / isoform expression
## pooja.singh09@gmail.com
## june 2024
## investigateing evolution change in gene expression

library(edgeR)
library(ouch)
library(reshape2)
library(ggplot2)
library(phytools)

setwd("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/14_ancestral_expression/")


## read phylogeny from file
phylogeny <-read.tree("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/12_phylogeny/RAxML_bipartitions.all.filt.L.vcf.SNP.gz.recode.mergedOnil.fmiss.raxmlout.T10.clean.tre") 

## plot phylogeny
plotTree(phylogeny)


################### OJ gene expression
## read data
svl<-read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/transcript_tpm_all_samples_onil_OJ_mean.txt",row.names=1)

## change this into a vector
svl<-as.matrix(svl)[39148,]
svl

##  ancestral state reconstruction
fit <- fastAnc(phylogeny,svl,vars=TRUE,CI=TRUE)

## projection of the reconstruction onto the edges of the phylogeny
obj<-contMap(phylogeny,svl,plot=FALSE)


plot(obj,legend=0.7*max(nodeHeights(phylogeny)),
     fsize=c(1,1))

phenogram(phylogeny,svl,fsize=0.6,spread.costs=c(1,0))

## plot traitgram with 95% CI
x<-fastBM(phylogeny)
fancyTree(phylogeny,type="phenogram95",x=x,spread.cost=c(1,0))
fancyTree(phylogeny,type="phenogram95",x=svl,spread.cost=c(1,0))

obj<-contMap(phylogeny,x)
obj<-contMap(phylogeny,svl)

plotphylogeny.wBars(phylogeny, exp(svl),scale=0.002,
               fsize=0.7,tip.labels=TRUE)

### OJ gene expression of convergent DEGs

svl<-read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/transcript_tpm_all_samples_onil_OJ_mean.txt",row.names=1)
data <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/gene_tpm_all_samples_onil_OJevee.tsv", header=T)
OG2gene <- data[, 1:2]
oj_conv_genes <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_convergent_genes_up_down.txt_gene_name.txt", header=T, sep="\t")

oj_conv_genes <- na.omit(oj_conv_genes)
oj_conv_genes$mstrg <- sapply(strsplit(oj_conv_genes$gene_id, "\\|"), `[`, 1)
oj_conv_genes_mstrg <- oj_conv_genes$mstrg
OG2geneconv <- OG2gene[OG2gene$GENE_NAME %in% oj_conv_genes_mstrg,]
svl_conv <- svl[rownames(svl) %in% OG2geneconv$OG,]
svl <- svl_conv


#MSTRG.38016
svl<-read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/transcript_tpm_all_samples_onil_OJ_mean.txt",row.names=1)
data <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/gene_tpm_all_samples_onil_OJevee.tsv", header=T)
OG2gene <- data[, 1:2]
oj_conv_genes <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/oj_convergent_genes_up_down.txt_gene_name.txt", header=T, sep="\t")
oj_conv_genes <- na.omit(oj_conv_genes)
oj_conv_genes_mstrg <- "MSTRG.38016"
OG2geneconv <- OG2gene[OG2gene$GENE_NAME %in% oj_conv_genes_mstrg,]
svl_conv <- svl[rownames(svl) %in% OG2geneconv$OG,]
svl <- svl_conv

## change this into a vector
svl<-as.matrix(svl)[1,]
svl

##  ancestral state reconstruction
fit <- fastAnc(phylogeny,svl,vars=TRUE,CI=TRUE)

## projection of the reconstruction onto the edges of the phylogeny
obj<-contMap(phylogeny,svl,plot=FALSE)


svg("odam_gene_expression_plots_phenogram.svg", h=6, w=12)
# Set up multi-plot layout with increased margins
par(mfrow=c(1, 2), mar=c(5, 5, 5, 5) + 0.1, oma=c(4, 4, 4, 4), oma=c(6, 6, 6, 6), mai=c(1.5, 1.5, 1.5, 1.5))  # Increase outer margins

# Plot the phylogeny with ancestral state reconstruction
plot(obj, legend=0.7*max(nodeHeights(phylogeny)), fsize=c(1, 1))

phenogram(phylogeny, svl, fsize=0.6, spread.costs=c(1, 0))
mtext("Evolutinary distance", side=1, line=3)
mtext("odam gene expression", side=2, line=3)
dev.off()

# Calculate the number of pages needed
num_genes <- length(oj_conv_genes_mstrg)
genes_per_page <- 4  # 2 rows and 2 columns
num_pages <- ceiling(num_genes / genes_per_page)

# Create a PDF file to save the plots
pdf("gene_expression_plots_phenogram.pdf")

for (page in 1:num_pages) {
  # Set up multi-plot layout
  par(mfrow=c(2, 2), mar=c(1, 1, 1, 1) + 0.1)  # Increase margins
  
  # Determine the range of genes for this page
  start_index <- (page - 1) * genes_per_page + 1
  end_index <- min(page * genes_per_page, num_genes)
  
  for (i in start_index:end_index) {
    gene <- oj_conv_genes_mstrg[i]
    gene_name <- oj_conv_genes[gene == oj_conv_genes$mstrg,]$Name
    print(gene)
    
    # Load the data
    svl <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/transcript_tpm_all_samples_onil_OJ_mean.txt", row.names=1)
    OG2geneconv <- OG2gene[OG2gene$GENE_NAME == gene,]
    svl_conv <- svl[rownames(svl) %in% OG2geneconv$OG,]
    
    # Check if svl_conv is not empty
    if (nrow(svl_conv) > 0) {
      svl <- svl_conv
      
      # Change this into a vector
      svl <- as.matrix(svl)[1,]
      
      # Ancestral state reconstruction
      fit <- fastAnc(phylogeny, svl, vars=TRUE, CI=TRUE)
      
      # Projection of the reconstruction onto the edges of the phylogeny
      obj <- tryCatch({
        contMap(phylogeny, svl, plot=FALSE)
      }, error = function(e) {
        message(paste("Error in contMap for gene:", gene_name, "-", e$message))
        NULL
      })
      
      if (!is.null(obj)) {
        # Plot the phylogeny with ancestral state reconstruction
        plot(obj, legend=0.7*max(nodeHeights(phylogeny)), fsize=c(0.6, 0.8), tip.labels=TRUE)  # Reduce font size
        #title(main=paste("", gene_name), cex.main=1.2)  # Adjust title size
        
        # Add legend manually if needed
        legend("topright", legend=paste("", gene_name), bty="n")
        
        # Plot the phenogram
        phenogram(phylogeny, svl, fsize=0.6, spread.costs=c(1,0), xlab="Time", ylab="Trait Value")
        #title(main=paste("Gene:", gene_name), cex.main=1.2)  # Adjust title size
      }
    } else {
      message(paste("No data found for gene:", gene_name))
    }
  }
}

# Close the PDF file
dev.off()

# Reset the plotting layout to default
par(mfrow=c(1, 1))



## plot traitgram with 95% CI
x<-fastBM(phylogeny)
fancyTree(phylogeny,type="phenogram95",x=x,spread.cost=c(1,0))
fancyTree(phylogeny,type="phenogram95",x=svl,spread.cost=c(1,0))

obj<-contMap(phylogeny,x)
obj<-contMap(phylogeny,svl)

plotphylogeny.wBars(phylogeny, exp(svl),scale=0.002,
                    fsize=0.7,tip.labels=TRUE)


### PJ gene expression of convergent DEGs

svl<-read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/transcript_tpm_all_samples_onil_PJ_mean.txt",row.names=1)
data <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/gene_tpm_all_samples_onil_PJevee.tsv", header=T)
OG2gene <- data[, 1:2]
PJ_conv_genes <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/PJ_convergent_genes_up_down.txt_gene_name.txt", header=T, sep="\t")

PJ_conv_genes <- na.omit(PJ_conv_genes)
PJ_conv_genes$mstrg <- sapply(strsplit(PJ_conv_genes$gene_id, "\\|"), `[`, 1)
PJ_conv_genes_mstrg <- PJ_conv_genes$mstrg
OG2geneconv <- OG2gene[OG2gene$GENE_NAME %in% PJ_conv_genes_mstrg,]
svl_conv <- svl[rownames(svl) %in% OG2geneconv$OG,]
svl <- svl_conv




# Calculate the number of pages needed
num_genes <- length(PJ_conv_genes_mstrg)
genes_per_page <- 4  # 2 rows and 2 columns
num_pages <- ceiling(num_genes / genes_per_page)

# Create a PDF file to save the plots
pdf("PJ_gene_expression_plots_phenogram.pdf")

for (page in 1:num_pages) {
  # Set up multi-plot layout
  par(mfrow=c(2, 2), mar=c(1, 1, 1, 1) + 0.1)  # Increase margins
  
  # Determine the range of genes for this page
  start_index <- (page - 1) * genes_per_page + 1
  end_index <- min(page * genes_per_page, num_genes)
  
  for (i in start_index:end_index) {
    gene <- PJ_conv_genes_mstrg[i]
    gene_name <- PJ_conv_genes[gene == PJ_conv_genes$mstrg,]$Name
    print(gene)
    
    # Load the data
    svl <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/transcript_tpm_all_samples_onil_PJ_mean.txt", row.names=1)
    OG2geneconv <- OG2gene[OG2gene$GENE_NAME == gene,]
    svl_conv <- svl[rownames(svl) %in% OG2geneconv$OG,]
    
    # Check if svl_conv is not empty
    if (nrow(svl_conv) > 0) {
      svl <- svl_conv
      
      # Change this into a vector
      svl <- as.matrix(svl)[1,]
      
      # Ancestral state reconstruction
      fit <- fastAnc(phylogeny, svl, vars=TRUE, CI=TRUE)
      
      # PrPJection of the reconstruction onto the edges of the phylogeny
      obj <- tryCatch({
        contMap(phylogeny, svl, plot=FALSE)
      }, error = function(e) {
        message(paste("Error in contMap for gene:", gene_name, "-", e$message))
        NULL
      })
      
      if (!is.null(obj)) {
        # Plot the phylogeny with ancestral state reconstruction
        plot(obj, legend=0.7*max(nodeHeights(phylogeny)), fsize=c(0.6, 0.8), tip.labels=TRUE)  # Reduce font size
        #title(main=paste("", gene_name), cex.main=1.2)  # Adjust title size
        
        # Add legend manually if needed
        legend("topright", legend=paste("", gene_name), bty="n")
        
        # Plot the phenogram
        phenogram(phylogeny, svl, fsize=0.6, spread.costs=c(1,0), xlab="Time", ylab="Trait Value")
        #title(main=paste("Gene:", gene_name), cex.main=1.2)  # Adjust title size
      }
    } else {
      message(paste("No data found for gene:", gene_name))
    }
  }
}

# Close the PDF file
dev.off()

# Reset the plotting layout to default
par(mfrow=c(1, 1))


## OJ PSI 
index <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/gene_tpm_all_samples_onil_OJevee.index.txt", header=F, comment.char="")
phylogeny <- read.delim("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/05_selection_geneexpression/evolution/st26_phylogeny_ouch.txt")
exp1<-read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/suppa/with_onil/filtered.gffcmp.all.annotated_onil_isoform.psi.txt",header=T, comment.char="", row.names=1, stringsAsFactors = F)
exp <- exp1[ , !(colnames(exp1) %like% ".PJ.")] ## OJ only

##get average of species
data=exp[1:ncol(exp)]
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
rownames(data_mean) <- rownames(exp)
write.table(data_mean,"filtered.gffcmp.all.annotated_onil_isoform.psi_mean.txt", quote=F, sep="\t")
svl <- data_mean
svl1 <- na.omit(svl)

## change this into a vector
svl<-as.matrix(svl1)[9040,]
svl

##  ancestral state reconstruction
fit <- fastAnc(phylogeny,svl,vars=TRUE,CI=TRUE)

## projection of the reconstruction onto the edges of the phylogeny
obj1<-contMap(phylogeny,svl,plot=FALSE)


plot(obj1,legend=0.7*max(nodeHeights(phylogeny)),
     fsize=c(0.7,0.9))

phenogram(phylogeny,svl1,fsize=0.6,spread.costs=c(1,0))

## plot traitgram with 95% CI
x<-fastBM(phylogeny)
fancyTree(phylogeny,type="phenogram95",x=x,spread.cost=c(1,0))
fancyTree(phylogeny,type="phenogram95",x=svl,spread.cost=c(1,0))

obj<-contMap(phylogeny,x)
obj<-contMap(phylogeny,svl)

plotphylogeny.wBars(phylogeny, exp(svl),scale=0.002,
                    fsize=0.7,tip.labels=TRUE)


