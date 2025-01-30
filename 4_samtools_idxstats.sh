#!bin/bash
# samtool idxstats - .bam mapping statistics
# March 2020
########### Change '-N' jobname, '-l h_vmem=' (storage), '-pe smp #' (number of cores)
########### samtools idxstats

FILE_DIR="/cl_tmp/singh_duenser/all_st26/star_assembly"

for i in `cat samtools_idx.txt`;


do

qsub -q all.q -pe smp 8 -l h_vmem=4G -cwd -V -N $i".idxstat" -o $FILE_DIR/$i".idxstat.out" -e $FILE_DIR/$i".idxstats.err" -b y samtools idxstats --threads 12 \
 $FILE_DIR/$i".Aligned.sortedByCoord.out.bam" &

done
