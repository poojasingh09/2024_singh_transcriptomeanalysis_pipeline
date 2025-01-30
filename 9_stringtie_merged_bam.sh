#!bin/bash
# picard
## run stringtie on single files

FILE_DIR="/cl_tmp/singh_duenser/all_st26/star_assembly"


qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N stringtie.all -e ./stringtie_single/stringtie.all.err -o ./stringtie_single/stringtie.all.out \
/usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/stringtie \
$FILE_DIR/*.bam --rf -f 0.15 -m 200 -a 10 -j 1 -c 2 -g 50 -M 0.95 -o ./stringtie_single/all.bam.gtf &

sleep 2

(base) singh@IT010043:/cl_tmp/singh_duenser/all_st26$ cat 9_stringtie_merged_bam.sh 
#!bin/bash

# anna.duenser@gmail.com
# March 2020
# stringtie
## run stringtie on merged bam files = MergedBam

FILE_DIR="/cl_tmp/singh_duenser/all_st26"


for i in `ls $FILE_DIR/merged_bam_files/*.bam | cut -d "/" -f 6`;

do

S_NAME="$(echo $i | cut -d'.' -f1-3)"


qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N "stringtieMB."$S_NAME -e ./stringtie_merged_bam/"stringtie_MergedBam."$S_NAME".err" \
-o ./stringtie_merged_bam/"stringtie_MergedBam."$S_NAME".out" \
/usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/stringtie \
$FILE_DIR/merged_bam_files/$i --rf -f 0.15 -m 200 -a 10 -j 1 -c 2 -g 50 -M 0.95 -o ./stringtie_merged_bam/$S_NAME".MergedBam.gtf" &

sleep 1

done
