

## pooja.singh09@gmail.com
## july 2024

library(ggplot2)
library(cowplot)
library(ggpubr)
library(ComplexHeatmap)
if (!require("devtools")) {
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
}
install_github("raivokolde/pheatmap")
library(pheatmap)
library(RColorBrewer)

## setwd
setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/10_GOenrichment/DEG')
#make figure 5 for paper

oj_lt <- read.table("LT_OJ_herbcarn_DESeq_fullresults_0.05.ourid.ensembl_topgo_fisher_weight_tests_BP_q0.05.txt", header=T,sep="\t")
oj_lt$trait <- "oj_lt"
oj_lm <- read.table("LM_OJ_herbcarn_DESeq_fullresults_0.05.ourid.ensembl_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
oj_lm $trait <- "oj_lm "
oj_lv <- read.table("LV_OJ_herbcarn_DESeq_fullresults_0.05.ourid.ensembl_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
oj_lv$trait <- "oj_lv"

pj_lt <- read.table("LT_PJ_herbcarn_DESeq_fullresults_0.05.ourid.ensembl_topgo_fisher_weight_tests_BP_q0.05.txt", header=T,sep="\t")
pj_lt$trait <- "pj_lt"
pj_lm <- read.table("LM_PJ_herbcarn_DESeq_fullresults_0.05.ourid.ensembl_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
pj_lm $trait <- "pj_lm "
pj_lv <- read.table("LV_PJ_herbcarn_DESeq_fullresults_0.05.ourid.ensembl_topgo_fisher_weight_tests_BP_q0.05.txt", header=T, sep="\t")
pj_lv$trait <- "pj_lv"

## just shortening some terms for plotting
oj_lt$Term[3] <- "exonucleolytic trimming"
oj_lm$Term[10] <- "positive regulation of cortical granule exocytosis"
pj_lv$Term[6]<- "positive regulation of cortical granule exocytosis"
pj_lv$Term[9] <- "adenylate cyclase-modulating GPCR signaling pathway"


## reformat data for heatmap
traits <- list(oj_lt,oj_lm,oj_lv, pj_lt,pj_lm,pj_lv)
names(traits) <- 

traits_new <- list()
counter=0
for (trait in traits){
  counter = counter + 1
  trait1 <- trait[,c(2,11)]
  trait1$topgoFisher.padjust <- -log10(trait1$topgoFisher.padjust)
  colnames(trait1) <- c("term", trait[1,12])
  traits_new[[counter]] <- trait1
  
}

input <- Reduce(function(x, y) merge(x, y, all=T, by=c("term")), traits_new, accumulate=F)


input1 <- input[,c(2:7)]
rownames(input1) <- input$term


svg("figure_s9_go_enrichment_heatmap_full.svg", height=7.5, width=8)
pheatmap(as.matrix(input1), labels_col= c("LT OJ", "LM OJ","LV OJ","LT PJ", "LM PJ", "LV PJ"), cluster_rows=FALSE, 
         cluster_cols=FALSE, na_col="white", cellwidth = 12, fontsize_row = 10,
         color = colorRampPalette((brewer.pal(n = 7, name ="Greens")))(10))

dev.off()

input2 <- input1[c(1,7,8,26,29,31,47),]

svg("figure_5d_go_enrichment_heatmap_subset.svg", height=4, width=4)
pheatmap(as.matrix(input2), labels_col= c("LT OJ", "LM OJ","LV OJ","LT PJ", "LM PJ", "LV PJ"), cluster_rows=FALSE, 
         cluster_cols=FALSE, na_col="white", cellwidth = 12, fontsize_row = 10,
         color = colorRampPalette((brewer.pal(n = 7, name ="Greens")))(10))

dev.off()


