#!bin/bash
# fastqc on trimmed reads (this is a trial on sequential job sumission)
# March 2020
########### Change '-N' jobname, '-l h_vmem=' (storage), '-pe smp #' (number of cores)
########### fastqc of trimmed reads


for i in `ls *.paired.fq.gz | cut -d "." -f 1-4`;

do

qsub -q all.q -pe smp 8 -l h_vmem=4G -cwd -V -N $i".fastqc" -o "log."$i".fastqc.out" -e "log."$i".fastqc.err" -b y fastqc -t 12 \
-o ./fastqc /cl_tmp/singh_duenser/all_st26/trimmomatic/$i".paired.fq.gz" &

done
