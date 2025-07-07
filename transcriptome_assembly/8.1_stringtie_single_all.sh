#!bin/bash

# picard
## run stringtie on single bam files

FILE_DIR="/cl_tmp/singh_duenser/all_st26/star_assembly"


for i in `ls $FILE_DIR/*.bam | cut -d "/" -f 6`;

do

S_NAME="$(echo $i | cut -d'.' -f1-4)"

qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N "stringtie."$S_NAME -e ./stringtie_single/$S_NAME".stringtie.single.err" -o ./stringtie_single/$S_NAME".stringtie.single.out" \
/usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/stringtie \
$FILE_DIR/$i --rf -f 0.15 -m 200 -a 10 -j 1 -c 2 -g 50 -M 0.95 -o ./stringtie_single/$S_NAME".gtf" &

sleep 2

done
