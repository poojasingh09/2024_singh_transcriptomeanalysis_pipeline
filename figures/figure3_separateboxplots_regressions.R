library(reshape2)
library(ggpubr)
library(cowplot)


## transcript boxplots

a <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/01_assembly_annotation/transcripts_persample_plot.txt", header=T)

a1 <- melt(a,id.vars='lake', measure.vars=c('oj','pj'))

my_comparisons <- list( c("LT", "LM"), c("LT", "LV"), c("LM", "LV"), c("LT", "RV"), c("LM", "RV"), c("LV", "RV"))

pj <- ggboxplot(a1[a1$variable == "pj",], x="lake", y="value", color="lake", palette = c(LT="goldenrod2", LM="#009E73", LV="#E49EC2", RV="#56B4E9"), add="jitter",  order = c("LT", "LM", "LV", "RV"), ylab="PJ isoforms")  + theme(legend.position = "none")
oj <- ggboxplot(a1[a1$variable == "oj",], x="lake", y="value", color="lake", palette = c(LT="goldenrod2", LM="#009E73", LV="#E49EC2", RV="#56B4E9"), add="jitter",  order = c("LT", "LM", "LV", "RV"), ylab="OJ isoforms")  + theme(legend.position = "none")

plot1 <- oj + stat_compare_means(comparisons = my_comparisons, method="wilcox.test")
plot2 <- pj + stat_compare_means(comparisons = my_comparisons, method="wilcox.test")

#svg("/Users/pooja/Desktop/Papers/2018_Postdoc_papers/RNAseq_newgrant/allst26/manuscript/figure3_transcriptcount.svg")
plot_grid(plot1, plot2, labels = c('A', 'B'), label_size = 12)
#dev.off()



## regression of transcripts and age of radiaiton

b <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/01_assembly_annotation/transcripts_persample_plotage.txt", header=T)

svg("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/manuscript/figure_S3_B.svg", h=6, w=12)
par(mfrow=c(1,2)) 

#rsq <- function (x, y) cor(x, y) ^ 2


#oj

mod1 = lm(b$age ~ b$oj, data = b)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
modsum$coefficients
my.p = modsum$coefficients[2,4]
plot(b$age, b$oj, pch=19, xlab="divergence time of radiation (years)", ylab="no. of isoforms expressed in OJ")
#r <- rsq(b$age, b$oj)
abline(lm(b$oj ~ b$age), lwd=2)
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 2000000, y = 30500, labels = mylabel)
mylabel2 = bquote(italic(p)== .(format(my.p, digits = 3)))
text(x = 2000000, y = 30200, labels = mylabel2)


#pj
mod1 = lm(b$age ~ b$pj, data = b)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
modsum$coefficients
my.p = modsum$coefficients[2,4]
plot(b$age, b$pj, pch=19, col="grey", xlab="divergence time of radiation (years)", ylab="no. of isoforms expressed in PJ")
#r <- rsq(b$age, b$pj)
abline(lm(b$pj ~ b$age), lwd=2, col="darkgrey")
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 2000000, y = 30000, labels = mylabel)
mylabel2 = bquote(italic(p)== .(format(my.p, digits = 3)))
text(x = 2000000, y = 29700, labels = mylabel2)

dev.off()





################ GENES



## gene boxplots




a <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/01_assembly_annotation/genes_persample_plot.txt", header=T)

a1 <- melt(a,id.vars='lake', measure.vars=c('oj','pj'))

my_comparisons <- list( c("LT", "LM"), c("LT", "LV"), c("LM", "LV"), c("LT", "RV"), c("LM", "RV"), c("LV", "RV"))

pj <- ggboxplot(a1[a1$variable == "pj",], x="lake", y="value", color="lake", palette = c(LT="goldenrod2", LM="#009E73", LV="#E49EC2", RV="#56B4E9"), add="jitter",  order = c("LT", "LM", "LV", "RV"), ylab="PJ genes") + theme(legend.position = "none")
oj <- ggboxplot(a1[a1$variable == "oj",], x="lake", y="value", color="lake", palette = c(LT="goldenrod2", LM="#009E73", LV="#E49EC2", RV="#56B4E9"), add="jitter",  order = c("LT", "LM", "LV", "RV"), ylab="OJ genes")  + theme(legend.position = "none")

plot5 <- oj + stat_compare_means(comparisons = my_comparisons, method="wilcox.test")
plot6 <- pj + stat_compare_means(comparisons = my_comparisons, method="wilcox.test")

#svg("/Users/pooja/Desktop/Papers/2018_Postdoc_papers/RNAseq_newgrant/allst26/manuscript/figure3_genecount.svg")
#plot_grid(plot5, plot6, labels = c('A', 'B'), label_size = 12)
#dev.off()


## regression of genes and age of radiaiton

b <- read.table("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/01_assembly_annotation/genes_persample_plotage.txt", header=T)

svg("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/manuscript/figure_S3_A.svg", h=6, w=12)
par(mfrow=c(1,2)) 

#rsq <- function (x, y) cor(x, y) ^ 2


#oj

mod1 = lm(b$age ~ b$oj, data = b)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
modsum$coefficients
my.p = modsum$coefficients[2,4]
plot(b$age, b$oj, pch=19, xlab="divergence time of radiation (years)", ylab="no. of genes expressed in OJ")
#r <- rsq(b$age, b$oj)
abline(lm(b$oj ~ b$age), lwd=2)
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 1500000, y = 28000, labels = mylabel)
mylabel2 = bquote(italic(p)== .(format(my.p, digits = 3)))
text(x = 1500000, y = 27700, labels = mylabel2)



#rsq <- function (x, y) cor(x, y) ^ 2

#pj
mod1 = lm(b$age ~ b$pj, data = b)
modsum = summary(mod1)
r2 = modsum$adj.r.squared
modsum$coefficients
my.p = modsum$coefficients[2,4]
plot(b$age, b$pj, pch=19, col="grey", xlab="divergence time of radiation (years)", ylab="no. of genes expressed in PJ")
#r <- rsq(b$age, b$pj)
abline(lm(b$pj ~ b$age), lwd=2, col="darkgrey")
mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
text(x = 1500000, y = 28000, labels = mylabel)
mylabel2 = bquote(italic(p)== .(format(my.p, digits = 3)))
text(x = 1500000, y = 27700, labels = mylabel2)

dev.off()

#### extra crap

oj + stat_compare_means(method="kruskal.test", label.y = 32000)+ stat_compare_means(label = "p.signif", method = "wilcox.test",ref.group = "lt")

pj + stat_compare_means(method="kruskal.test", label.y = 32000)+ stat_compare_means(label = "p.signif", method = "wilcox.test",ref.group = "lt")

#both <- ggboxplot(a1, x="lake", y="value", color="variable", palette = c(LT="goldenrod2", LM="#009E73", LV="#E49EC2", RV="#56B4E9"), add="jitter", order = c("lt", "lm", "lv", "rv"), ylab="genes", xlab="lake")


### plots 3 barplots of expressed and uniquely expressed genes and isoforms
library(ggplot2)
library(ggpattern)
library(RColorBrewer)



setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/01_assembly_annotation')

input <- read.table("figure3_input.txt", header=T, sep="\t")

### with pattern colouring


#svg("figure3_barplot_v1.svg")
plot9 <- ggplot(input) + 
  aes(x = type2, 
      fill = type1, 
      weight = percentage, 
      pattern = type2) +  # Differentiate patterns by type2
  scale_x_discrete(drop = TRUE) + 
  geom_bar_pattern(position = "dodge", 
                   pattern_density = 0.5,  
                   pattern_fill = "black",  # Keeps patterns visible
                   pattern_color = NA) +  # Removes pattern outlines
  scale_pattern_manual(values = c("none", "stripe", "none", "stripe")) +  # Define patterns
  scale_fill_brewer(palette = "Pastel3") +  # Fill colors from Pastel1 palette
  facet_wrap(~type1, strip.position = "bottom", 
             scales = "free_x", 
             nrow = 1) +  # Keep facets in a single row
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.grid = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        strip.placement = "inside")

#dev.off()

## figure 3 final 2025

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/manuscript')

svg("figure3_2025.svg", h=8, w=8)
row1 <- plot_grid(plot1, plot2, ncol = 2, labels = c('A', 'B'), label_size = 14)
plot_grid(row1, plot9, ncol = 1, rel_heights = c(1.4, 2), labels = c('A', 'C'), label_size = 14)
dev.off()

svg("figure_S2_2025.svg", h=4, w=6)
plot_grid(plot5, plot6, ncol = 2,labels = c('A', 'B'), label_size = 14)

dev.off()



