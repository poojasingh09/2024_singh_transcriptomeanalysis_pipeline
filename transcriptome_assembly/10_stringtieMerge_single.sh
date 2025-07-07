#!bin/bash
# stringtie --merge
## merge single .gtf files for each species

FILE_DIR="/cl_tmp/singh_duenser/all_st26/stringtie_single"


for i in `ls $FILE_DIR/*.1.gtf | cut -d "/" -f 6 | cut -d'.' -f1-3`;

do

qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N "stringtieMmerge."$i -e $FILE_DIR/stringtieMerge_single/"stringtieM."$i".err" \
-o $FILE_DIR/stringtieMerge_single/"stringtieM."$i".out" \
/usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/stringtie --merge \
-F 1.0 -T 1.0 $FILE_DIR/$i"."*".gtf" \
-o $FILE_DIR/stringtieMerge_single/"stringtieMerge_single."$i".gtf" &

sleep 1
echo "yippie!"

done
