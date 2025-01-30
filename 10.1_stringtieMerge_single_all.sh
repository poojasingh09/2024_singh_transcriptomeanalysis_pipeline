#!bin/bash
# stringtie --merge
## merge single .gtf files for each species

FILE_DIR="/cl_tmp/singh_duenser/all_st26/stringtie_single"


qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N stringtieMmerge.allgtf -e $FILE_DIR/stringtieMerge_single/stringtieM.allgtf.err \
-o $FILE_DIR/stringtieMerge_single/stringtieM.allgtf.out \
/usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/stringtie --merge \
-F 1.0 -T 1.0 $FILE_DIR/*.gtf \
-o $FILE_DIR/stringtieMerge_single/stringtieMerge.allgtf.gtf &

sleep 1
echo "yippie!"
