#!/bin/bash
# multiqc
#### input: fastqc from initial .fq files, fastqc files after trimmomatic, mapped .bam files
FILE_DIR="/cl_tmp/singh_duenser/all_st26"

qsub -q all.q -b y -cwd -N multiqc -l h_vmem=4G -pe smp 8 -o ./multiqc/"log.multiqc.out" -e ./multiqc/"log.multiqc.err" /usr/people/EDVZ/duenser/.local/bin/multiqc \
--interactive $FILE_DIR/fastqc_raw/ $FILE_DIR/trimmomatic/fastqc/ $FILE_DIR/star_assembly/
(base) singh@IT010043:/cl_tmp/singh_duenser/all_st26$ cat 6_rename_bam.sh 

#!bin/bash
# April 2020
# rename larvae and put a ".L." in the name

FILE_DIR="/cl_tmp/singh_duenser/all_st26/star_assembly"

# loop through .bam files
for i in `ls $FILE_DIR/*.bam | cut -d "/" -f 6`;

do
# cut only first part off
FIRST="$(echo $i | cut -d '.' -f 1)"
# cut the rest of needed name of
LAST="$(echo $i | cut -d '.' -f 2-3)"
# create new name
NEW_NAME=$FIRST".L."$LAST".bam"

mv $FILE_DIR/$i $FILE_DIR/$NEW_NAME

done
