# psingh
# run prepDE.py which generates a count matrix for DESeq

DIR="/cl_tmp/singh_duenser/all_st26/count_matrices"

qsub -l h_vmem=10G -cwd -V -N prepDE_st26 -o $DIR/prepDE.log -e $DIR/prepDE.err -b y /usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/prepDE.py -i $DIR/prepDE_samples -g $DIR/gene_count_matrix_prepDE.csv -t $DIR/transcript_count_matrix_prepDE.csv  -l 125


(base) singh@IT010043:/cl_tmp/singh_duenser/all_st26$ cat 17_prepDE.sh 
Display all 483 possibilities? (y or n)
(base) singh@IT010043:/cl_tmp/singh_duenser/all_st26$ cat 17_prepDE.sh 
# psingh
#aug2020
# run prepDE.py which generates a count matrix for DESeq

DIR="/cl_tmp/singh_duenser/all_st26/count_matrices"

qsub -l h_vmem=10G -cwd -V -N prepDE_st26 -o $DIR/prepDE.log -e $DIR/prepDE.err -b y /usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/prepDE.py -i $DIR/prepDE_samples -g $DIR/gene_count_matrix_prepDE.csv -t $DIR/transcript_count_matrix_prepDE.csv  -l 125

