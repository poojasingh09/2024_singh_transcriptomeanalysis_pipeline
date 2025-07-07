library(scales)
library(ggplot2)
library(ggrepel)

# install.packages("scales")
# install.packages("ggplot2")
# install.packages("ggrepel")

setwd('/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/03_DTU/suppa')
#Load the dspi file (output from SUPPA diffSplice)
dpsi <- read.table(file="oj_lt.lakes.dpsi",sep="\t")
head(dpsi)
colnames(dpsi) <- c("dPSI","p_value")

#Load the psi file (output from SUPPA generateEvents)
psi_events <- read.table(file="filtered.gffcmp.all.annotated_isoform.psi_speciesmean_oj.csv", sep="\t", header=T)
##colnames(psi_events) <- c("CTRL1","CTRL2","CTRL3","KD1","KD2","KD3")


#Load the tpm of the events (output from SUPPA diffSplice, activating flag --save_tpm_events)
##event_TPM <- read.table(file="~/TRA2_diffSplice_avglogtpm.tab",sep="\t",header = TRUE)
##colnames(event_TPM) <- c("event","mean_TPM")


#Merge dpsi and psi
merge1 <- merge(dpsi,psi_events,by="row.names")
#merge2 <- merge(merge1,event_TPM,by.x="Row.names",by.y="event")
rownames(merge1) <- merge1$Row.names
final_table <- merge1
final_table <- final_table[,-1]
final_table <- final_table[!is.nan(final_table$dPSI),]
final_table$cpval <- p.adjust(final_table$p_value, method = "bonferroni")
final_table$log10pval <- -log10(final_table$p_value)
final_table$sig <- "not sig"
final_table[final_table$p_value < 0.05,]$sig <- "sig"
# final_table[, c(2:13)] <- sapply(final_table[, c(2:13)], function(x) as.numeric(gsub(",", ".", x)))
##final_table$logRNAc <-final_table$mean_TPM

dim(final_table)
final_table$gene_mstrg <- sapply(strsplit(rownames(final_table),";"), `[`, 1)
final_table$isoform_mstrg <- sapply(strsplit(rownames(final_table),";"), `[`, 2)
final_table$isoform <- sapply(strsplit(final_table$isoform_mstrg, ".", fixed=TRUE), tail, 1)

#annotate genes 
dpsi_anno <- read.table(file="oj_lt.lakes.dpsi.sig.genes.mstrg_id.ensembl",sep="\t")
final_table <- merge(dpsi_anno, final_table, by.x="V1", by.y="gene_mstrg", all=T)

dim(final_table)

gff <- readGFF("/Users/singhpoo/Desktop/Papers/2019_postdoc_papers/2019_RNAseq_newgrant/2023_allst26/analysis/new_ref/02_DEG/deseq_new/2023/O_niloticus_UMD_NMBU.99.gff3.genes")
final_table <- merge(gff, final_table, by.x="gene_id", by.y="V2", all=T)
final_table$anno <- paste(final_table$Name, final_table$isoform, sep="_")
head(final_table)
dim(final_table)

#Volcano plot

#jpeg(file = "~/volcano_TRA2.jpeg")

p <- ggplot(final_table, aes(x=dPSI, y=log10pval, color=sig))
p + geom_point() + geom_text_repel(aes(x=dPSI, y=log10pval), label = ifelse((final_table$dPSI < -0.5 | final_table$dPSI > 0.5), final_table$anno, ""))+
  geom_vline(xintercept=c(-0.5,0.5), linetype="solid", size=1) +
  geom_hline(yintercept=1.3, size=1) +
  xlab(expression(~Delta~PSI)) + ylab("-log10(p-value)") + 
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    legend.position="none",
    legend.text=element_blank(),
    legend.title=element_blank(),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20)) +
  scale_color_manual(values=c("sig" = "red", "not sig" = "gray40", "nan" = "gray80")) + 
  scale_x_continuous(breaks=pretty_breaks(n=5)) +
  scale_y_continuous(breaks=pretty_breaks(n=5)) +
  guides(fill = guide_legend(reverse = FALSE))

#dev.off()




