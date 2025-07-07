## PSI plots of species HvsC
## May 2024
# pooja.singh09@gmail.com
# errors arise from lib conflicsts so beware! :)

# libraries
library(tidyverse)
library(ggpubr)
library(gdata)
library(patchwork)
library(ggrepel)
library(pheatmap)
library(data.table)
library(cowplot)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(tidyr)


DIR.SUPPA   <- "/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/SUPPA/"
PATH.OUT    <- "/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/14_ancestral_expression/"
JAW         <- c("oj", "pj")
LAKES         <- c("lv", "lm","lt")
DATE        <- format(Sys.time(), "%Y%m%d")


setwd(PATH.OUT)

# load psivec file and rename single colums
psi_values <- read.delim("filtered.gffcmp.all.annotated_isoform.psi.txt") 
head(psi_values)
psi_values_oj <- psi_values[,grepl("OJ", colnames(psi_values))]
psi_values_oj$id <- rownames(psi_values_oj)

psi_values_pj <- psi_values[,grepl("PJ", colnames(psi_values))]
psi_values_pj$id <- rownames(psi_values_pj)


# color palette
# colorBlindGrey  <- c("#C5C1C1" ,"#464343", "#E69F00", "#56B4E9", "#009E73",
#                      "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

COL <- c(LT="goldenrod2", LM="#009E73", LV="#E49EC2", RV="#56B4E9")

res <- list()

########## LV OJ
# load dpsi file with only significant comparison to select genes of interest
#sig <-  read_delim(file.path(paste(DIR.SUPPA, JAW, "_", LAKES, ".lakes.dpsi.sig", sep = "")), col_names =  T, delim = "\t") %>% dplyr::select(filtered-filtered_dPSI, filtered-filtered_p-val)

sig <- read.delim(file.path(paste(DIR.SUPPA, "oj", "_", "lv", ".lakes.dpsi.sig", sep = "")))
sig$id <-rownames(sig)
nrow(sig)


# calculatemedian PSI of biological replicates and manipulate gene name
sig_df <- psi_values_oj %>%
  dplyr::filter(psi_values_oj$id %in% sig$id) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    RV =median(c(Aa.L.OJ.1,Aa.L.OJ.2,Aa.L.OJ.3,Aa.L.OJ.4,Aa.L.OJ.5,Ab.L.OJ.1,Ab.L.OJ.2,Ab.L.OJ.3,Ab.L.OJ.4,Ab.L.OJ.5), na.rm = T),
    LT =median(c(Ch.L.OJ.1,Ch.L.OJ.2,Ch.L.OJ.3,Ch.L.OJ.4,Ch.L.OJ.5,Gp.L.OJ.1,Gp.L.OJ.2,Gp.L.OJ.3,Gp.L.OJ.4,Gp.L.OJ.5,Pf.L.OJ.1,Pf.L.OJ.2,Pf.L.OJ.3,Pf.L.OJ.4,Pf.L.OJ.5,Pp.L.OJ.1,Pp.L.OJ.2,Pp.L.OJ.3,Pp.L.OJ.4,Pp.L.OJ.5,Sd.L.OJ.1,Sd.L.OJ.2,Sd.L.OJ.3,Sd.L.OJ.4,Sd.L.OJ.5,Tm.L.OJ.1,Tm.L.OJ.2,Tm.L.OJ.3,Tm.L.OJ.4,Tm.L.OJ.5), na.rm = T),
    LM =median(c(Ah.L.OJ.1,Ah.L.OJ.2,Ah.L.OJ.3,Ah.L.OJ.4,Ah.L.OJ.5,Lt.L.OJ.1,Lt.L.OJ.2,Lt.L.OJ.3,Lt.L.OJ.4,Lt.L.OJ.5,Ptb.L.OJ.1,Ptb.L.OJ.2,Ptb.L.OJ.3,Ptb.L.OJ.4,Ptb.L.OJ.5,Py.L.OJ.1,Py.L.OJ.2,Py.L.OJ.3,Py.L.OJ.4,Py.L.OJ.5,Sf.L.OJ.1,Sf.L.OJ.2,Sf.L.OJ.3,Sf.L.OJ.4,Sf.L.OJ.5,Tt.L.OJ.1,Tt.L.OJ.2,Tt.L.OJ.3,Tt.L.OJ.4,Tt.L.OJ.5), na.rm = T),
    LV =median(c(Gh.L.OJ.1,Gh.L.OJ.2,Gh.L.OJ.3,Gh.L.OJ.4,Gh.L.OJ.5,Ht.L.OJ.1,Ht.L.OJ.2,Ht.L.OJ.3,Ht.L.OJ.4,Ht.L.OJ.5,Ml.L.OJ.1,Ml.L.OJ.2,Ml.L.OJ.3,Ml.L.OJ.4,Ml.L.OJ.5,No.L.OJ.1,No.L.OJ.2,No.L.OJ.3,No.L.OJ.4,No.L.OJ.5,Prp.L.OJ.1,Prp.L.OJ.2,Prp.L.OJ.3,Prp.L.OJ.4,Prp.L.OJ.5,Ps.L.OJ.1,Ps.L.OJ.2,Ps.L.OJ.3,Ps.L.OJ.4,Ps.L.OJ.5), na.rm = T),
  ) %>% 
  dplyr::left_join(sig, by = "id")

plot_rvlv <- sig_df[,c("id", "RV", "LV")]
plot_rvlv$direction <- ifelse(plot_rvlv$LV > plot_rvlv$RV &
                                plot_rvlv$LV - plot_rvlv$RV > 0.2, TRUE,FALSE)
plot_rvlv$novel <- ifelse(plot_rvlv$RV == 0 & plot_rvlv$LV > plot_rvlv$RV &
                                plot_rvlv$LV - plot_rvlv$RV > 0, TRUE,FALSE)

nrow(na.omit(plot_rvlv[plot_rvlv$direction==TRUE,]))
nrow(na.omit(plot_rvlv[plot_rvlv$direction==FALSE,]))

nrow(na.omit(plot_rvlv[plot_rvlv$novel==TRUE,]))
nrow(na.omit(plot_rvlv[plot_rvlv$novel==FALSE,]))

res$lv_oj <- plot_rvlv

plotdata <- plot_rvlv %>% dplyr::select(RV, LV, id, direction, novel) %>%
  pivot_longer(RV:LV, names_to = "Species", values_to = "PSI")



plotdata$Species <- factor(plotdata$Species, levels=unique(plotdata$Species))
                              
p1 <- ggplot(plotdata, aes(x = Species, y = PSI, group = id, colour = direction, linetype = novel)) +
  geom_point(show.legend = F)+
  geom_line(show.legend = F)+
  scale_color_manual(values=c("lightgrey", "#E49EC2"))+
  scale_x_discrete(expand = c(0.04,0))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



########## LM OJ
# load dpsi file with only significant comparison to select genes of interest
#sig <-  read_delim(file.path(paste(DIR.SUPPA, JAW, "_", LAKES, ".lakes.dpsi.sig", sep = "")), col_names =  T, delim = "\t") %>% dplyr::select(filtered-filtered_dPSI, filtered-filtered_p-val)

sig <- read.delim(file.path(paste(DIR.SUPPA, "oj", "_", "lm", ".lakes.dpsi.sig", sep = "")))
sig$id <-rownames(sig)
nrow(sig)


# calculatemedian PSI of biological replicates and manipulate gene name
sig_df <- psi_values_oj %>%
  
  dplyr::filter(psi_values_oj$id %in% sig$id) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    RV =median(c(Aa.L.OJ.1,Aa.L.OJ.2,Aa.L.OJ.3,Aa.L.OJ.4,Aa.L.OJ.5,Ab.L.OJ.1,Ab.L.OJ.2,Ab.L.OJ.3,Ab.L.OJ.4,Ab.L.OJ.5), na.rm = T),
    LT =median(c(Ch.L.OJ.1,Ch.L.OJ.2,Ch.L.OJ.3,Ch.L.OJ.4,Ch.L.OJ.5,Gp.L.OJ.1,Gp.L.OJ.2,Gp.L.OJ.3,Gp.L.OJ.4,Gp.L.OJ.5,Pf.L.OJ.1,Pf.L.OJ.2,Pf.L.OJ.3,Pf.L.OJ.4,Pf.L.OJ.5,Pp.L.OJ.1,Pp.L.OJ.2,Pp.L.OJ.3,Pp.L.OJ.4,Pp.L.OJ.5,Sd.L.OJ.1,Sd.L.OJ.2,Sd.L.OJ.3,Sd.L.OJ.4,Sd.L.OJ.5,Tm.L.OJ.1,Tm.L.OJ.2,Tm.L.OJ.3,Tm.L.OJ.4,Tm.L.OJ.5), na.rm = T),
    LM =median(c(Ah.L.OJ.1,Ah.L.OJ.2,Ah.L.OJ.3,Ah.L.OJ.4,Ah.L.OJ.5,Lt.L.OJ.1,Lt.L.OJ.2,Lt.L.OJ.3,Lt.L.OJ.4,Lt.L.OJ.5,Ptb.L.OJ.1,Ptb.L.OJ.2,Ptb.L.OJ.3,Ptb.L.OJ.4,Ptb.L.OJ.5,Py.L.OJ.1,Py.L.OJ.2,Py.L.OJ.3,Py.L.OJ.4,Py.L.OJ.5,Sf.L.OJ.1,Sf.L.OJ.2,Sf.L.OJ.3,Sf.L.OJ.4,Sf.L.OJ.5,Tt.L.OJ.1,Tt.L.OJ.2,Tt.L.OJ.3,Tt.L.OJ.4,Tt.L.OJ.5), na.rm = T),
    LM =median(c(Gh.L.OJ.1,Gh.L.OJ.2,Gh.L.OJ.3,Gh.L.OJ.4,Gh.L.OJ.5,Ht.L.OJ.1,Ht.L.OJ.2,Ht.L.OJ.3,Ht.L.OJ.4,Ht.L.OJ.5,Ml.L.OJ.1,Ml.L.OJ.2,Ml.L.OJ.3,Ml.L.OJ.4,Ml.L.OJ.5,No.L.OJ.1,No.L.OJ.2,No.L.OJ.3,No.L.OJ.4,No.L.OJ.5,Prp.L.OJ.1,Prp.L.OJ.2,Prp.L.OJ.3,Prp.L.OJ.4,Prp.L.OJ.5,Ps.L.OJ.1,Ps.L.OJ.2,Ps.L.OJ.3,Ps.L.OJ.4,Ps.L.OJ.5), na.rm = T),
  ) %>% 
  dplyr::left_join(sig, by = "id")



plot_rvLM <- sig_df[,c("id", "RV", "LM")]
plot_rvLM$direction <- ifelse(plot_rvLM$LM > plot_rvLM$RV &
                                plot_rvLM$LM - plot_rvLM$RV > 0.2, TRUE,FALSE)
plot_rvLM$novel <- ifelse(plot_rvLM$RV == 0 & plot_rvLM$LM > plot_rvLM$RV &
                            plot_rvLM$LM - plot_rvLM$RV > 0, TRUE,FALSE)

nrow(na.omit(plot_rvLM[plot_rvLM$direction==TRUE,]))
nrow(na.omit(plot_rvLM[plot_rvLM$direction==FALSE,]))

nrow(na.omit(plot_rvLM[plot_rvLM$novel==TRUE,]))
nrow(na.omit(plot_rvLM[plot_rvLM$novel==FALSE,]))

res$lm_oj <- plot_rvLM

plotdata <- plot_rvLM %>% dplyr::select(RV, LM, id, direction, novel) %>%
  pivot_longer(RV:LM, names_to = "Species", values_to = "PSI")

plotdata$Species <- factor(plotdata$Species, levels=unique(plotdata$Species))

p2 <-ggplot(plotdata, aes(x = Species, y = PSI, group = id, colour = direction, linetype = novel)) +
  geom_point(show.legend = F)+
  geom_line(show.legend = F)+
  scale_color_manual(values=c("lightgrey", "#009E73"))+
  scale_x_discrete(expand = c(0.04,0))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


########## LT OJ
# load dpsi file with only significant comparison to select genes of interest
#sig <-  read_delim(file.path(paste(DIR.SUPPA, JAW, "_", LAKES, ".lakes.dpsi.sig", sep = "")), col_names =  T, delim = "\t") %>% dplyr::select(filtered-filtered_dPSI, filtered-filtered_p-val)

sig <- read.delim(file.path(paste(DIR.SUPPA, "oj", "_", "lt", ".lakes.dpsi.sig", sep = "")))
sig$id <-rownames(sig)
nrow(sig)


# calculatemedian PSI of biological replicates and manipulate gene name
sig_df <- psi_values_oj %>%
  dplyr::filter(rownames(psi_values_oj) %in% sig$id) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    RV =median(c(Aa.L.OJ.1,Aa.L.OJ.2,Aa.L.OJ.3,Aa.L.OJ.4,Aa.L.OJ.5,Ab.L.OJ.1,Ab.L.OJ.2,Ab.L.OJ.3,Ab.L.OJ.4,Ab.L.OJ.5), na.rm = T),
    LT =median(c(Ch.L.OJ.1,Ch.L.OJ.2,Ch.L.OJ.3,Ch.L.OJ.4,Ch.L.OJ.5,Gp.L.OJ.1,Gp.L.OJ.2,Gp.L.OJ.3,Gp.L.OJ.4,Gp.L.OJ.5,Pf.L.OJ.1,Pf.L.OJ.2,Pf.L.OJ.3,Pf.L.OJ.4,Pf.L.OJ.5,Pp.L.OJ.1,Pp.L.OJ.2,Pp.L.OJ.3,Pp.L.OJ.4,Pp.L.OJ.5,Sd.L.OJ.1,Sd.L.OJ.2,Sd.L.OJ.3,Sd.L.OJ.4,Sd.L.OJ.5,Tm.L.OJ.1,Tm.L.OJ.2,Tm.L.OJ.3,Tm.L.OJ.4,Tm.L.OJ.5), na.rm = T),
    LM =median(c(Ah.L.OJ.1,Ah.L.OJ.2,Ah.L.OJ.3,Ah.L.OJ.4,Ah.L.OJ.5,Lt.L.OJ.1,Lt.L.OJ.2,Lt.L.OJ.3,Lt.L.OJ.4,Lt.L.OJ.5,Ptb.L.OJ.1,Ptb.L.OJ.2,Ptb.L.OJ.3,Ptb.L.OJ.4,Ptb.L.OJ.5,Py.L.OJ.1,Py.L.OJ.2,Py.L.OJ.3,Py.L.OJ.4,Py.L.OJ.5,Sf.L.OJ.1,Sf.L.OJ.2,Sf.L.OJ.3,Sf.L.OJ.4,Sf.L.OJ.5,Tt.L.OJ.1,Tt.L.OJ.2,Tt.L.OJ.3,Tt.L.OJ.4,Tt.L.OJ.5), na.rm = T),
    LT =median(c(Gh.L.OJ.1,Gh.L.OJ.2,Gh.L.OJ.3,Gh.L.OJ.4,Gh.L.OJ.5,Ht.L.OJ.1,Ht.L.OJ.2,Ht.L.OJ.3,Ht.L.OJ.4,Ht.L.OJ.5,Ml.L.OJ.1,Ml.L.OJ.2,Ml.L.OJ.3,Ml.L.OJ.4,Ml.L.OJ.5,No.L.OJ.1,No.L.OJ.2,No.L.OJ.3,No.L.OJ.4,No.L.OJ.5,Prp.L.OJ.1,Prp.L.OJ.2,Prp.L.OJ.3,Prp.L.OJ.4,Prp.L.OJ.5,Ps.L.OJ.1,Ps.L.OJ.2,Ps.L.OJ.3,Ps.L.OJ.4,Ps.L.OJ.5), na.rm = T),
  ) %>% 
  dplyr::left_join(sig, by = "id")



plot_rvLT <- sig_df[,c("id", "RV", "LT")]
plot_rvLT$direction <- ifelse(plot_rvLT$LT > plot_rvLT$RV &
                                plot_rvLT$LT - plot_rvLT$RV > 0.2, TRUE,FALSE)
plot_rvLT$novel <- ifelse(plot_rvLT$RV == 0 & plot_rvLT$LT > plot_rvLT$RV &
                            plot_rvLT$LT - plot_rvLT$RV > 0, TRUE,FALSE)

nrow(na.omit(plot_rvLT[plot_rvLT$direction==TRUE,]))
nrow(na.omit(plot_rvLT[plot_rvLT$direction==FALSE,]))

nrow(na.omit(plot_rvLT[plot_rvLT$novel==TRUE,]))
nrow(na.omit(plot_rvLT[plot_rvLT$novel==FALSE,]))

res$lt_oj <- plot_rvLT

plotdata <- plot_rvLT %>% dplyr::select(RV, LT, id, direction, novel) %>%
  pivot_longer(RV:LT, names_to = "Species", values_to = "PSI")

plotdata$Species <- factor(plotdata$Species, levels=unique(plotdata$Species))

p3 <-ggplot(plotdata, aes(x = Species, y = PSI, group = id, colour = direction, linetype = novel)) +
  geom_point(show.legend = F)+
  geom_line(show.legend = F)+
  scale_color_manual(values=c("lightgrey", "goldenrod2"))+
  scale_x_discrete(expand = c(0.04,0))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))




########## LV PJ
# load dpsi file with only significant comparison to select genes of interest
#sig <-  read_delim(file.path(paste(DIR.SUPPA, JAW, "_", LAKES, ".lakes.dpsi.sig", sep = "")), col_names =  T, delim = "\t") %>% dplyr::select(filtered-filtered_dPSI, filtered-filtered_p-val)

sig <- read.delim(file.path(paste(DIR.SUPPA, "pj", "_", "lv", ".lakes.dpsi.sig", sep = "")))
sig$id <-rownames(sig)
nrow(sig)


# calculatemedian PSI of biological replicates and manipulate gene name
sig_df <- psi_values_pj %>%
  dplyr::filter(rownames(psi_values_pj) %in% sig$id) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    RV =median(c(Aa.L.PJ.1,Aa.L.PJ.2,Aa.L.PJ.3,Aa.L.PJ.4,Aa.L.PJ.5,Ab.L.PJ.1,Ab.L.PJ.2,Ab.L.PJ.3,Ab.L.PJ.4,Ab.L.PJ.5), na.rm = T),
    LT =median(c(Ch.L.PJ.1,Ch.L.PJ.2,Ch.L.PJ.3,Ch.L.PJ.4,Ch.L.PJ.5,Gp.L.PJ.1,Gp.L.PJ.2,Gp.L.PJ.3,Gp.L.PJ.4,Gp.L.PJ.5,Pf.L.PJ.1,Pf.L.PJ.2,Pf.L.PJ.3,Pf.L.PJ.4,Pf.L.PJ.5,Pp.L.PJ.1,Pp.L.PJ.2,Pp.L.PJ.3,Pp.L.PJ.4,Pp.L.PJ.5,Sd.L.PJ.1,Sd.L.PJ.2,Sd.L.PJ.3,Sd.L.PJ.4,Sd.L.PJ.5,Tm.L.PJ.1,Tm.L.PJ.2,Tm.L.PJ.3,Tm.L.PJ.4,Tm.L.PJ.5), na.rm = T),
    LM =median(c(Ah.L.PJ.1,Ah.L.PJ.2,Ah.L.PJ.3,Ah.L.PJ.4,Ah.L.PJ.5,Lt.L.PJ.1,Lt.L.PJ.2,Lt.L.PJ.3,Lt.L.PJ.4,Lt.L.PJ.5,Ptb.L.PJ.1,Ptb.L.PJ.2,Ptb.L.PJ.3,Ptb.L.PJ.4,Ptb.L.PJ.5,Py.L.PJ.1,Py.L.PJ.2,Py.L.PJ.3,Py.L.PJ.4,Py.L.PJ.5,Sf.L.PJ.1,Sf.L.PJ.2,Sf.L.PJ.3,Sf.L.PJ.4,Sf.L.PJ.5,Tt.L.PJ.1,Tt.L.PJ.2,Tt.L.PJ.3,Tt.L.PJ.4,Tt.L.PJ.5), na.rm = T),
    LV =median(c(Gh.L.PJ.1,Gh.L.PJ.2,Gh.L.PJ.3,Gh.L.PJ.4,Gh.L.PJ.5,Ht.L.PJ.1,Ht.L.PJ.2,Ht.L.PJ.3,Ht.L.PJ.4,Ht.L.PJ.5,Ml.L.PJ.1,Ml.L.PJ.2,Ml.L.PJ.3,Ml.L.PJ.4,Ml.L.PJ.5,No.L.PJ.1,No.L.PJ.2,No.L.PJ.3,No.L.PJ.4,No.L.PJ.5,Prp.L.PJ.1,Prp.L.PJ.2,Prp.L.PJ.3,Prp.L.PJ.4,Prp.L.PJ.5,Ps.L.PJ.1,Ps.L.PJ.2,Ps.L.PJ.3,Ps.L.PJ.4,Ps.L.PJ.5), na.rm = T),
  ) %>% 
  dplyr::left_join(sig, by = "id")



plot_rvlv <- sig_df[,c("id", "RV", "LV")]
plot_rvlv$direction <- ifelse(plot_rvlv$LV > plot_rvlv$RV &
                                plot_rvlv$LV - plot_rvlv$RV > 0.2, TRUE,FALSE)
plot_rvlv$novel <- ifelse(plot_rvlv$RV == 0 & plot_rvlv$LV > plot_rvlv$RV &
                            plot_rvlv$LV - plot_rvlv$RV > 0, TRUE,FALSE)

nrow(na.omit(plot_rvlv[plot_rvlv$direction==TRUE,]))
nrow(na.omit(plot_rvlv[plot_rvlv$direction==FALSE,]))

nrow(na.omit(plot_rvlv[plot_rvlv$novel==TRUE,]))
nrow(na.omit(plot_rvlv[plot_rvlv$novel==FALSE,]))
res$lv_pj <- plot_rvlv

plotdata <- plot_rvlv %>% dplyr::select(RV, LV, id, direction, novel) %>%
  pivot_longer(RV:LV, names_to = "Species", values_to = "PSI")

plotdata$Species <- factor(plotdata$Species, levels=unique(plotdata$Species))

p4 <-ggplot(plotdata, aes(x = Species, y = PSI, group = id, colour = direction, linetype = novel)) +
  geom_point(show.legend = F)+
  geom_line(show.legend = F)+
  scale_color_manual(values=c("lightgrey", "#E49EC2"))+
  scale_x_discrete(expand = c(0.04,0))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



########## LM PJ
# load dpsi file with only significant comparison to select genes of interest
#sig <-  read_delim(file.path(paste(DIR.SUPPA, JAW, "_", LAKES, ".lakes.dpsi.sig", sep = "")), col_names =  T, delim = "\t") %>% dplyr::select(filtered-filtered_dPSI, filtered-filtered_p-val)

sig <- read.delim(file.path(paste(DIR.SUPPA, "pj", "_", "lm", ".lakes.dpsi.sig", sep = "")))
sig$id <-rownames(sig)
nrow(sig)


# calculatemedian PSI of biological replicates and manipulate gene name
sig_df <- psi_values_pj %>%
  dplyr::filter(rownames(psi_values_pj) %in% sig$id) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    RV =median(c(Aa.L.PJ.1,Aa.L.PJ.2,Aa.L.PJ.3,Aa.L.PJ.4,Aa.L.PJ.5,Ab.L.PJ.1,Ab.L.PJ.2,Ab.L.PJ.3,Ab.L.PJ.4,Ab.L.PJ.5), na.rm = T),
    LT =median(c(Ch.L.PJ.1,Ch.L.PJ.2,Ch.L.PJ.3,Ch.L.PJ.4,Ch.L.PJ.5,Gp.L.PJ.1,Gp.L.PJ.2,Gp.L.PJ.3,Gp.L.PJ.4,Gp.L.PJ.5,Pf.L.PJ.1,Pf.L.PJ.2,Pf.L.PJ.3,Pf.L.PJ.4,Pf.L.PJ.5,Pp.L.PJ.1,Pp.L.PJ.2,Pp.L.PJ.3,Pp.L.PJ.4,Pp.L.PJ.5,Sd.L.PJ.1,Sd.L.PJ.2,Sd.L.PJ.3,Sd.L.PJ.4,Sd.L.PJ.5,Tm.L.PJ.1,Tm.L.PJ.2,Tm.L.PJ.3,Tm.L.PJ.4,Tm.L.PJ.5), na.rm = T),
    LM =median(c(Ah.L.PJ.1,Ah.L.PJ.2,Ah.L.PJ.3,Ah.L.PJ.4,Ah.L.PJ.5,Lt.L.PJ.1,Lt.L.PJ.2,Lt.L.PJ.3,Lt.L.PJ.4,Lt.L.PJ.5,Ptb.L.PJ.1,Ptb.L.PJ.2,Ptb.L.PJ.3,Ptb.L.PJ.4,Ptb.L.PJ.5,Py.L.PJ.1,Py.L.PJ.2,Py.L.PJ.3,Py.L.PJ.4,Py.L.PJ.5,Sf.L.PJ.1,Sf.L.PJ.2,Sf.L.PJ.3,Sf.L.PJ.4,Sf.L.PJ.5,Tt.L.PJ.1,Tt.L.PJ.2,Tt.L.PJ.3,Tt.L.PJ.4,Tt.L.PJ.5), na.rm = T),
    LM =median(c(Gh.L.PJ.1,Gh.L.PJ.2,Gh.L.PJ.3,Gh.L.PJ.4,Gh.L.PJ.5,Ht.L.PJ.1,Ht.L.PJ.2,Ht.L.PJ.3,Ht.L.PJ.4,Ht.L.PJ.5,Ml.L.PJ.1,Ml.L.PJ.2,Ml.L.PJ.3,Ml.L.PJ.4,Ml.L.PJ.5,No.L.PJ.1,No.L.PJ.2,No.L.PJ.3,No.L.PJ.4,No.L.PJ.5,Prp.L.PJ.1,Prp.L.PJ.2,Prp.L.PJ.3,Prp.L.PJ.4,Prp.L.PJ.5,Ps.L.PJ.1,Ps.L.PJ.2,Ps.L.PJ.3,Ps.L.PJ.4,Ps.L.PJ.5), na.rm = T),
  ) %>% 
  dplyr::left_join(sig, by = "id")



plot_rvLM <- sig_df[,c("id", "RV", "LM")]
plot_rvLM$direction <- ifelse(plot_rvLM$LM > plot_rvLM$RV &
                                plot_rvLM$LM - plot_rvLM$RV > 0.2, TRUE,FALSE)
plot_rvLM$novel <- ifelse(plot_rvLM$RV == 0 & plot_rvLM$LM > plot_rvLM$RV &
                            plot_rvLM$LM - plot_rvLM$RV > 0, TRUE,FALSE)

nrow(na.omit(plot_rvLM[plot_rvLM$direction==TRUE,]))
nrow(na.omit(plot_rvLM[plot_rvLM$direction==FALSE,]))

nrow(na.omit(plot_rvLM[plot_rvLM$novel==TRUE,]))
nrow(na.omit(plot_rvLM[plot_rvLM$novel==FALSE,]))
res$lm_pj <- plot_rvLM

plotdata <- plot_rvLM %>% dplyr::select(RV, LM, id, direction, novel) %>%
  pivot_longer(RV:LM, names_to = "Species", values_to = "PSI")

plotdata$Species <- factor(plotdata$Species, levels=unique(plotdata$Species))

p5 <-ggplot(plotdata, aes(x = Species, y = PSI, group = id, colour = direction, linetype = novel)) +
  geom_point(show.legend = F)+
  geom_line(show.legend = F)+
  scale_color_manual(values=c("lightgrey", "#009E73"))+
  scale_x_discrete(expand = c(0.04,0))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


########## LT PJ
# load dpsi file with only significant comparison to select genes of interest
#sig <-  read_delim(file.path(paste(DIR.SUPPA, JAW, "_", LAKES, ".lakes.dpsi.sig", sep = "")), col_names =  T, delim = "\t") %>% dplyr::select(filtered-filtered_dPSI, filtered-filtered_p-val)

sig <- read.delim(file.path(paste(DIR.SUPPA, "pj", "_", "lt", ".lakes.dpsi.sig", sep = "")))
sig$id <-rownames(sig)
nrow(sig)


# calculatemedian PSI of biological replicates and manipulate gene name
sig_df <- psi_values_pj %>%
  dplyr::filter(rownames(psi_values_pj) %in% sig$id) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    RV =median(c(Aa.L.PJ.1,Aa.L.PJ.2,Aa.L.PJ.3,Aa.L.PJ.4,Aa.L.PJ.5,Ab.L.PJ.1,Ab.L.PJ.2,Ab.L.PJ.3,Ab.L.PJ.4,Ab.L.PJ.5), na.rm = T),
    LT =median(c(Ch.L.PJ.1,Ch.L.PJ.2,Ch.L.PJ.3,Ch.L.PJ.4,Ch.L.PJ.5,Gp.L.PJ.1,Gp.L.PJ.2,Gp.L.PJ.3,Gp.L.PJ.4,Gp.L.PJ.5,Pf.L.PJ.1,Pf.L.PJ.2,Pf.L.PJ.3,Pf.L.PJ.4,Pf.L.PJ.5,Pp.L.PJ.1,Pp.L.PJ.2,Pp.L.PJ.3,Pp.L.PJ.4,Pp.L.PJ.5,Sd.L.PJ.1,Sd.L.PJ.2,Sd.L.PJ.3,Sd.L.PJ.4,Sd.L.PJ.5,Tm.L.PJ.1,Tm.L.PJ.2,Tm.L.PJ.3,Tm.L.PJ.4,Tm.L.PJ.5), na.rm = T),
    LM =median(c(Ah.L.PJ.1,Ah.L.PJ.2,Ah.L.PJ.3,Ah.L.PJ.4,Ah.L.PJ.5,Lt.L.PJ.1,Lt.L.PJ.2,Lt.L.PJ.3,Lt.L.PJ.4,Lt.L.PJ.5,Ptb.L.PJ.1,Ptb.L.PJ.2,Ptb.L.PJ.3,Ptb.L.PJ.4,Ptb.L.PJ.5,Py.L.PJ.1,Py.L.PJ.2,Py.L.PJ.3,Py.L.PJ.4,Py.L.PJ.5,Sf.L.PJ.1,Sf.L.PJ.2,Sf.L.PJ.3,Sf.L.PJ.4,Sf.L.PJ.5,Tt.L.PJ.1,Tt.L.PJ.2,Tt.L.PJ.3,Tt.L.PJ.4,Tt.L.PJ.5), na.rm = T),
    LT =median(c(Gh.L.PJ.1,Gh.L.PJ.2,Gh.L.PJ.3,Gh.L.PJ.4,Gh.L.PJ.5,Ht.L.PJ.1,Ht.L.PJ.2,Ht.L.PJ.3,Ht.L.PJ.4,Ht.L.PJ.5,Ml.L.PJ.1,Ml.L.PJ.2,Ml.L.PJ.3,Ml.L.PJ.4,Ml.L.PJ.5,No.L.PJ.1,No.L.PJ.2,No.L.PJ.3,No.L.PJ.4,No.L.PJ.5,Prp.L.PJ.1,Prp.L.PJ.2,Prp.L.PJ.3,Prp.L.PJ.4,Prp.L.PJ.5,Ps.L.PJ.1,Ps.L.PJ.2,Ps.L.PJ.3,Ps.L.PJ.4,Ps.L.PJ.5), na.rm = T),
  ) %>% 
  left_join(sig, by = "id")



plot_rvLT <- sig_df[,c("id", "RV", "LT")]
plot_rvLT$direction <- ifelse(plot_rvLT$LT > plot_rvLT$RV &
                                plot_rvLT$LT - plot_rvLT$RV > 0.2, TRUE,FALSE)
plot_rvLT$novel <- ifelse(plot_rvLT$RV == 0 & plot_rvLT$LT > plot_rvLT$RV &
                            plot_rvLT$LT - plot_rvLT$RV > 0, TRUE,FALSE)

nrow(na.omit(plot_rvLT[plot_rvLT$direction==TRUE,]))
nrow(na.omit(plot_rvLT[plot_rvLT$direction==FALSE,]))

nrow(na.omit(plot_rvLT[plot_rvLT$novel==TRUE,]))
nrow(na.omit(plot_rvLT[plot_rvLT$novel==FALSE,]))

res$lt_pj <- plot_rvLT
plotdata <- plot_rvLT %>% dplyr::select(RV, LT, id, direction, novel) %>%
  pivot_longer(RV:LT, names_to = "Species", values_to = "PSI")

plotdata$Species <- factor(plotdata$Species, levels=unique(plotdata$Species))

p6 <-ggplot(plotdata, aes(x = Species, y = PSI, group = id, colour = direction, linetype = novel)) +
  geom_point(show.legend = F)+
  geom_line(show.legend = F)+
  scale_color_manual(values=c("lightgrey", "goldenrod2"))+
  scale_x_discrete(expand = c(0.04,0))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


svg("oj_psi_rvlake_plotsv2.svg", width=5, height=3)
plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12, nrow=1)
dev.off()

svg("pj_psi_rvlake_plotsv2.svg")
plot_grid(p4, p5, p6, labels = c('A', 'B', 'C'), label_size = 12, nrow=1)
dev.off()

svg("ojpj_psi_rvlake_plotsv2.svg")
plot_grid(p1, p2, p3, p4, p5, p6, labels = c('A', 'B', 'C', 'D', 'E','F'), label_size = 12, nrow=2)
dev.off()


res_out <- rbindlist(res, fill=TRUE, idcol = "comparison")
res_out$gene_id <-  sapply(strsplit(res_out$id, split=";", fixed = TRUE), `[`, 1)
res_out$transcript_id <-  sapply(strsplit(res_out$id, split=";", fixed = TRUE), `[`, 2)

#################read anno so plots can be annotated
ourgff <- readGFF("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/filtered.gffread.gffcmp.all.annotated.gtf")
ourgff1 <- ourgff[ourgff$type == "transcript",]
head(ourgff1)
res_out_anno <- merge(res_out, ourgff1, by=c("transcript_id"))
head(res_out_anno)

gff <- readGFF("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/O_niloticus_UMD_NMBU.99.gff3.genes")
gff$name_full <-  sapply(strsplit(gff$description, split=" [", fixed = TRUE), `[`, 1)
head(gff)


res_out_anno1 <- merge(res_out_anno, gff, by.x="gene_name", by.y="gene_id")
head(res_out_anno1)
#write.table(res_out_anno1, "psi_novel_standing_results.txt", quote=F, row.names=F)

res_out_anno1[res_out_anno1$comparison == "lv_oj" & res_out_anno1$novel == TRUE,]
res_out_anno1[res_out_anno1$comparison == "lm_oj" & res_out_anno1$novel == TRUE,]
res_out_anno1[res_out_anno1$comparison == "lt_oj" & res_out_anno1$novel == TRUE,]
res_out_anno1[res_out_anno1$comparison == "lv_pj" & res_out_anno1$novel == TRUE,]
res_out_anno1[res_out_anno1$comparison == "lm_pj" & res_out_anno1$novel == TRUE,]
res_out_anno1[res_out_anno1$comparison == "lt_pj" & res_out_anno1$novel == TRUE,]


res_out_anno1[res_out_anno1$comparison == "lv_oj" & res_out_anno1$direction == TRUE,]
res_out_anno1[res_out_anno1$comparison == "lm_oj" & res_out_anno1$direction == TRUE,]
res_out_anno1[res_out_anno1$comparison == "lt_oj" & res_out_anno1$direction == TRUE,]
res_out_anno1[res_out_anno1$comparison == "lv_pj" & res_out_anno1$direction == TRUE,]
res_out_anno1[res_out_anno1$comparison == "lm_pj" & res_out_anno1$direction == TRUE,]
res_out_anno1[res_out_anno1$comparison == "lt_pj" & res_out_anno1$direction== TRUE,]

#### highlight gene models of novel 

### ggtanscript

library(magrittr)
library(dplyr)
library(ggtranscript)
library(ggplot2)

gtf <- rtracklayer::import("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/filtered.gffread.gffcmp.all.annotated.gtf")


a_df=as.data.frame(gtf)
a_df %>% head()

a_exons <- a_df %>% dplyr::filter(type == "exon" & gene_id == "MSTRG.8460")

svg("mstrg.8460.genemodel.svg",h=4,w=6)
a_exons %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range(
    aes(fill = transcript_id)
  ) +
  geom_intron(
    data = to_intron(a_exons, "transcript_id"),
    aes(strand = strand)
  ) + scale_fill_manual(values = c("white", "lightgrey"))+ theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + xlab("position (bp)")

dev.off()

### plot exp

oj_mean <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/suppa/filtered.gffcmp.all.annotated_isoform.psi_speciesmean_oj.csv")
pj_mean <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/suppa/filtered.gffcmp.all.annotated_isoform.psi_speciesmean_pj.csv")

oj_mean1 <-  oj_mean[, c("Tm","Sd","Pp","Pf","Ch","Gp","Tt","Lt","Ptb","Py","Sf","Ah", 
                                                "Ps","Ml","No","Gh","Ht","Prp",
                                                "Ab","Aa")]

pj_mean1 <-  pj_mean[, c("Tm","Sd","Pp","Pf","Ch","Gp","Tt","Lt","Ptb","Py","Sf","Ah", 
                         "Ps","Ml","No","Gh","Ht","Prp",
                         "Ab","Aa")]

novel_oj_psi <- oj_mean1 %>% filter(rownames(oj_mean1) == "MSTRG.8460;MSTRG.8460.2")
d = melt(novel_oj_psi)

svg("mstrg.8460.2.expression.svg",h=4,w=8)
ggplot(data = d, aes(x = variable, y = value, fill=variable)) + 
  geom_col() + scale_fill_manual(values=c("goldenrod2","goldenrod2","goldenrod2","goldenrod2","goldenrod2","goldenrod2",
                              "#009E73","#009E73","#009E73","#009E73","#009E73","#009E73",
                              "#E49EC2","#E49EC2","#E49EC2","#E49EC2","#E49EC2","#E49EC2",
                              "#56B4E9","#56B4E9")) + theme_bw() + 
                              theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                              axis.text.x = element_text(size = 12),  # Increase x-axis text size
                              axis.text.y = element_text(size = 12)) +
                              xlab("species") + ylab("mean PSI") + ylim(0,1)
dev.off()                              
                              
                              
novel_pj_psi <- pj_mean1 %>% filter(rownames(pj_mean1) == "MSTRG.13205;MSTRG.13205.2")
d = melt(novel_pj_psi)
ggplot(data = d, aes(x = variable, y = value, fill=variable)) + 
  geom_col() + scale_fill_manual(values=c("goldenrod2","goldenrod2","goldenrod2","goldenrod2","goldenrod2","goldenrod2",
                                          "#009E73","#009E73","#009E73","#009E73","#009E73","#009E73",
                                          "#E49EC2","#E49EC2","#E49EC2","#E49EC2","#E49EC2","#E49EC2",
                                          "#56B4E9","#56B4E9")) + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("species") + ylab("mean PSI") + ylim(0,1)




#####
COL <- c(Tm="goldenrod2", Sd="goldenrod2", Pp="goldenrod2", Pf="goldenrod2", Ch="goldenrod2", Gp="goldenrod2",
         Tt="#009E73", Lt="#009E73", Ptb="#009E73", Py="#009E73",Sf="#009E73", Ah="#009E73", 
         Ps="#E49EC2",Ml="#E49EC2",No="#E49EC2",Gh="#E49EC2",Ht="#E49EC2",Prp="#E49EC2",
         Ab="#56B4E9", Ab="#56B4E9")
