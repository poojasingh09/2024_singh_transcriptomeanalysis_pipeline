#filtering annotation.gtf files for u and p class codes
#PS: ensure that you only remove monoexonic GENES and not monoexonic transcripts from the annotation

#copy files over to gffread folder coz i dont want to mess with the annotated gtfs in the gffcomapr efolder
#cp /cl_tmp/singh_duenser/all_st26/gffcompare/gffcmp*annotated.gtf /cl_tmp/singh_duenser/all_st26/gffcompare/gffread


FILE.DIR = "/cl_tmp/singh_duenser/all_st26/gffcompare/gffread" #!change this
setwd(FILE.DIR) 

#install.packages("tidyverse", repos = 'https://cran.wu.ac.at/')
library(tidyverse)
#install.packages("rtracklayer", repos = 'https://cran.wu.ac.at/')
library(rtracklayer)
library(dplyr)
library("BiocGenerics")


DATE <- format(Sys.time(), "%m%d%Y")

all.files <- list.files(FILE.DIR, full.names = FALSE)
gffcmp_files <- as.vector(grep("*gtf", all.files, value = TRUE))

my_log <- file(paste("log_filterscritp_R_", DATE, ".txt", sep = ""))
sink(my_log, append = TRUE, type = "output", split = TRUE)

filter_function <- function(filename_1){
  print("Loaded Files:")
  print(filename_1)
  l1 <- length(filename_1)
  print("Number of Files:")
  print(l1)
  
  for (i in filename_1){
    print("Filtering:")
    print(i)

    # import annotation file
    # import gffread file with long introns filtered out
    all <- import(i)
    all <- as_tibble(all)
    print("Imported Rows:")
    print(nrow(all))
    # filter all single exons not in the reference (class code "u")
    # filter class p"
    #  do not mess with the filtering steps. n() is needed to keep grouping alive a bit longer!
    
    # group by gene_id, to exclude genes that also include monoexonic isoforms and only 
    # filter monoexonic transcripts that ar not in the reference
    all_without_singleUexons <- all %>% group_by(gene_id) %>% filter(!(n() == 2 && class_code == "u"))
    print("Removed single exon with class_code 'u'")
    print(nrow(all_without_singleUexons)-nrow(all))
    
    # filter_o <- filtered_exons_introns %>% group_by(transcript_id) %>% filter(!(n() >= 2 && class_code == "o"))
    filter_p <- all_without_singleUexons %>% group_by(transcript_id) %>% filter(!(n() >= 2 && class_code == "p"))
    print("Removed transcripts with class_code 'p'")
    print(nrow(filter_p)-nrow(all))
    
    filename <- paste("filtered.", i, sep = "")
    print("Saving")
    print(filename)
    export(filter_p, filename, format="GTF")
    print(Sys.time())
  }
  return(TRUE)
}


filter_function(gffcmp_files)

sink()


##########################################################################################
quit(save="no")
