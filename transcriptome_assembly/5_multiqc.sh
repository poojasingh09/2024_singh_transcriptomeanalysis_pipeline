#!bin/bash

# samtool idxstats - .bam mapping statistics
########### Change '-N' jobname, '-l h_vmem=' (storage), '-pe smp #' (number of cores)
########### samtools idxstats

FILE_DIR="/cl_tmp/singh_duenser/all_st26/star_assembly"

for i in `cat samtools_idx.txt`;


do

qsub -q all.q -pe smp 8 -l h_vmem=4G -cwd -V -N $i".idxstat" -o $FILE_DIR/$i".idxstat.out" -e $FILE_DIR/$i".idxstats.err" -b y samtools idxstats --threads 12 \
 $FILE_DIR/$i".Aligned.sortedByCoord.out.bam" &

done

(base) singh@IT010043:/cl_tmp/singh_duenser/all_st26$ cat 5_multiqc.sh 
#!/bin/bash

# anna.duenser@gmail.com
# March 2020
# multiqc
#### input: fastqc from initial .fq files, fastqc files after trimmomatic, mapped .bam files
FILE_DIR="/cl_tmp/singh_duenser/all_st26"

qsub -q all.q -b y -cwd -N multiqc -l h_vmem=4G -pe smp 8 -o ./multiqc/"log.multiqc.out" -e ./multiqc/"log.multiqc.err" /usr/people/EDVZ/duenser/.local/bin/multiqc \
--interactive $FILE_DIR/fastqc_raw/ $FILE_DIR/trimmomatic/fastqc/ $FILE_DIR/star_assembly/
