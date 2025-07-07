
## gene expression level of DEGs and DSGs and overlap
## nov 2024


library(ggVennDiagram)

setwd("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/")

deg_lt_oj <- read.table("LT_OJ_herbcarn_DESeq_fullresults_0.05.ourid.mstrg", col.names="name")
deg_lm_oj <- read.table("LM_OJ_herbcarn_DESeq_fullresults_0.05.ourid.mstrg", col.names="name")
deg_lv_oj <- read.table("LV_OJ_herbcarn_DESeq_fullresults_0.05.ourid.mstrg", col.names="name")

deg_lt_pj <- read.table("LT_PJ_herbcarn_DESeq_fullresults_0.05.ourid.mstrg", col.names="name")
deg_lm_pj <- read.table("LM_PJ_herbcarn_DESeq_fullresults_0.05.ourid.mstrg", col.names="name")
deg_lv_pj <- read.table("LV_PJ_herbcarn_DESeq_fullresults_0.05.ourid.mstrg", col.names="name")



setwd("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/suppa/")

dsg_lt_oj <- read.table("oj_lt.lakes.dpsi.sig.genes.mstrg_id", col.names="name")
dsg_lm_oj <- read.table("oj_lm.lakes.dpsi.sig.genes.mstrg_id", col.names="name")
dsg_lv_oj <- read.table("oj_lv.lakes.dpsi.sig.genes.mstrg_id", col.names="name")

dsg_lt_pj <- read.table("pj_lt.lakes.dpsi.sig.genes.mstrg_id", col.names="name")
dsg_lm_pj <- read.table("pj_lm.lakes.dpsi.sig.genes.mstrg_id", col.names="name")
dsg_lv_pj <- read.table("pj_lv.lakes.dpsi.sig.genes.mstrg_id", col.names="name")


# Chart
x <- list(deg = unique(deg_lt_oj$name), dsg = unique(dsg_lt_oj$name))
ltoj <- ggVennDiagram(x)

x <- list(deg = unique(deg_lm_oj$name), dsg = unique(dsg_lm_oj$name))
lmoj <- ggVennDiagram(x)


x <- list(deg = unique(deg_lv_oj$name), dsg = unique(dsg_lv_oj$name))
lvoj <- ggVennDiagram(x)


# Chart
x <- list(deg = unique(deg_lt_pj$name), dsg = unique(dsg_lt_pj$name))
ltpj <- ggVennDiagram(x)

x <- list(deg = unique(deg_lm_pj$name), dsg = unique(dsg_lm_pj$name))
lmpj <- ggVennDiagram(x)


x <- list(deg =unique(deg_lv_pj$name), dsg = unique(dsg_lv_pj$name))
lvpj <- ggVennDiagram(x)

## gene express

meanoj <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/gene_count_matrix_prepDE_speciesmean_oj.csv")
meanoj$mean <- rowMeans(meanoj)
name <- data.frame(do.call('rbind', strsplit(as.character(rownames(meanoj)),'|',fixed=TRUE)))[1]
meanoj1 <- cbind(meanoj, name)
meanoj1 <- meanoj1[,c(21,22)]


deg_lt_oj_exp <- meanoj1[meanoj1$X1 %in% deg_lt_oj$name, ]
dim(deg_lt_oj_exp)
deg_lm_oj_exp <- meanoj1[meanoj1$X1 %in% deg_lm_oj$name, ]
dim(deg_lm_oj_exp)
deg_lv_oj_exp <- meanoj1[meanoj1$X1 %in% deg_lv_oj$name, ]
dim(deg_lv_oj_exp)

dsg_lt_oj_exp <- meanoj1[meanoj1$X1 %in% dsg_lt_oj$name, ]
dim(dsg_lt_oj_exp)
dsg_lm_oj_exp <- meanoj1[meanoj1$X1 %in% dsg_lm_oj$name, ]
dim(dsg_lm_oj_exp)
dsg_lv_oj_exp <- meanoj1[meanoj1$X1 %in% dsg_lv_oj$name, ]
dim(dsg_lv_oj_exp)

meanpj <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/gene_count_matrix_prepDE_speciesmean_pj.csv")
meanpj$mean <- rowMeans(meanpj)
name <- data.frame(do.call('rbind', strsplit(as.character(rownames(meanpj)),'|',fixed=TRUE)))[1]
meanpj1 <- cbind(meanpj, name)
meanpj1 <- meanpj1[,c(21,22)]


deg_lt_pj_exp <- meanpj1[meanpj1$X1 %in% deg_lt_pj$name, ]
dim(deg_lt_pj_exp)
deg_lm_pj_exp <- meanpj1[meanpj1$X1 %in% deg_lm_pj$name, ]
dim(deg_lm_pj_exp)
deg_lv_pj_exp <- meanpj1[meanpj1$X1 %in% deg_lv_pj$name, ]
dim(deg_lv_pj_exp)

dsg_lt_pj_exp <- meanpj1[meanpj1$X1 %in% dsg_lt_pj$name, ]
dim(dsg_lt_pj_exp)
dsg_lm_pj_exp <- meanpj1[meanpj1$X1 %in% dsg_lm_pj$name, ]
dim(dsg_lm_pj_exp)
dsg_lv_pj_exp <- meanpj1[meanpj1$X1 %in% dsg_lv_pj$name, ]
dim(dsg_lv_pj_exp)


boxplot(c(meanoj1[1], deg_lt_oj_exp[1], deg_lm_oj_exp[1], deg_lv_oj_exp[1],
          dsg_lt_oj_exp[1], dsg_lm_oj_exp[1], dsg_lv_oj_exp[1]),outline=F)


boxplot(c(meanpj1[1], deg_lt_pj_exp[1], deg_lm_pj_exp[1], deg_lv_pj_exp[1],
          dsg_lt_pj_exp[1], dsg_lm_pj_exp[1], dsg_lv_pj_exp[1]),outline=F)


##suppmat S9
library(patchwork)
library(ggplot2)

# Arrange the plots in a 2x3 grid
combined_plot <- (ltoj | lmoj | lvoj) /
  (ltpj | lmpj | lvpj)

# Display the combined plot

ggsave("overlap_dsg_deg_figure_s9.svg", plot = combined_plot, width = 12, height = 8, units = "in")


