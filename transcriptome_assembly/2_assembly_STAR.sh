#!/bin/bash

#pooja.singh09@gmail.com
#RNA star assembles reads into transcriptsusing a reference genome guide
########### Change '-N' jobname, '-l h_vmem=' (storage), '-pe smp #' (number of cores)
########### mapping to reference O. niloticus reference genome with 5 bioreprs per species per jaw


for i in `cat inputfiles`;

do

qsub -q all.q -pe smp 8 -l h_vmem=4G -cwd -V -N $i".STAR" -e $i".STAR.log" -o $i".STAR.err" -b y STAR --runThreadN 8 \
 --genomeDir /cl_tmp/singh_duenser/reference \
 --readFilesCommand gunzip -c \
 --readFilesIn $i".1.paired.fq.gz" $i".2.paired.fq.gz" \
 --outSAMtype BAM SortedByCoordinate \
 --twopassMode None \
 --outFileNamePrefix ./star_assembly/$i"." \
 --quantMode GeneCounts \
 --outSAMstrandField intronMotif \
 --outSAMattrIHstart 1 \
 --outSAMattributes NH HI AS nM NM MD MC jM jI ch XS \
 --outSAMprimaryFlag OneBestScore \
 --outSAMmapqUnique 60 \
 --outSAMunmapped Within  \
 --outFilterIntronStrands RemoveInconsistentStrands \
 --outBAMsortingBinsN 50  \
 --limitBAMsortRAM 4000000000 &

done
