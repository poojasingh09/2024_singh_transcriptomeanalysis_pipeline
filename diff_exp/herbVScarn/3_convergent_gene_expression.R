##pooja.singh09@gmail.com
## finding convergent gene expression


### in shell #####

fgrep -f OJ_convergent_genes.txt LM_OJ_herbcarn_DESeq_fullresults_0.05.txt.sig > OJ_convergent_genes.LM
fgrep -f OJ_convergent_genes.txt LV_OJ_herbcarn_DESeq_fullresults_0.05.txt.sig > OJ_convergent_genes.LV
fgrep -f OJ_convergent_genes.txt LT_OJ_herbcarn_DESeq_fullresults_0.05.txt.sig > OJ_convergent_genes.LT


fgrep -f PJ_convergent_genes.txt LM_PJ_herbcarn_DESeq_fullresults_0.05.txt.sig > PJ_convergent_genes.LM
fgrep -f PJ_convergent_genes.txt LV_PJ_herbcarn_DESeq_fullresults_0.05.txt.sig > PJ_convergent_genes.LV
fgrep -f PJ_convergent_genes.txt LT_PJ_herbcarn_DESeq_fullresults_0.05.txt.sig > PJ_convergent_genes.LT


##### in R ####

###OJ herb vs carn

lt <- read.table("OJ_convergent_genes.LT", header=F)
lm <- read.table("OJ_convergent_genes.LM", header=F)
lv <- read.table("OJ_convergent_genes.LV", header=F)


lt <- lt[,c(1,3)]
lm <- lm[,c(1,3)]
lv <- lv[,c(1,3)]

one <- merge(lt, lm, by="V1")
two <- merge(one, lv, by="V1")

colnames(two) <- c("gene", "lt", "lm", "lv")

carn_up <- two[(two$lt > 0) & (two$lm > 0) & (two$lv > 0),]
herb_up <- two[(two$lt < 0) & (two$lm < 0) & (two$lv < 0),]

write.table(carn_up, "OJ_convergent_genes_carn_up,txt", sep="\t", quote=F, row.names=F)
write.table(herb_up, "OJ_convergent_genes_herb_up,txt", sep="\t", quote=F, row.names=F)

###PJ herb vs carn

lt <- read.table("PJ_convergent_genes.LT", header=F)
lm <- read.table("PJ_convergent_genes.LM", header=F)
lv <- read.table("PJ_convergent_genes.LV", header=F)


lt <- lt[,c(1,3)]
lm <- lm[,c(1,3)]
lv <- lv[,c(1,3)]

one_pj <- merge(lt, lm, by="V1")
two_pj <- merge(one_pj, lv, by="V1")

colnames(two_pj) <- c("gene", "lt", "lm", "lv")

carn_up <- two_pj[(two_pj$lt > 0) & (two_pj$lm > 0) & (two_pj$lv > 0),]
herb_up <- two_pj[(two_pj$lt < 0) & (two_pj$lm < 0) & (two_pj$lv < 0),]

write.table(carn_up, "PJ_convergent_genes_carn_up,txt", sep="\t", quote=F, row.names=F)
write.table(herb_up, "PJ_convergent_genes_herb_up,txt", sep="\t", quote=F, row.names=F)


