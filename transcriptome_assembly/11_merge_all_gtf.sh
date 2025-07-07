#!bin/bash
#poojasingh
#aug2020
#merging gtfs originating from single .bam files (and then merged) and the ones from merged .bam files
#merged for OJ and PJ separately for comparison and then OJ and PJ will be merged for final output
#the reason for this is that I want to compare the number of transcripts assmbled for OJ and PJ too
# I also merged all samples at once versus first merging OJ and PJ, just to see if there was a difference
# I find more trancripts when I merged OJ and PJ separately than when I merge all (77k versus 84k)




FILE_DIR="/cl_tmp/singh_duenser/all_st26"


#OJ (69519 transcripts)
qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N all_gtf_mergeOJ \
-e $FILE_DIR/all_gtf_merged/all_gtf.OJerr -o $FILE_DIR/all_gtf_merged/all_gtf.OJ.out \
/usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/stringtie --merge -F 1.0 -T 1.0 \
$FILE_DIR/stringtie_single/stringtieMerge_single/stringtieMerge_single.*.OJ.gtf \
$FILE_DIR/stringtie_merged_bam/*.OJ.MergedBam.gtf \
-p 10 -o $FILE_DIR/all_gtf_merged/all_gtf_merged.OJ.gtf 

sleep 1

#PJ (67057 transcripts)
qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N all_gtf_mergePJ \
-e $FILE_DIR/all_gtf_merged/all_gtf.PJerr -o $FILE_DIR/all_gtf_merged/all_gtf.PJout \
/usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/stringtie --merge -F 1.0 -T 1.0 \
$FILE_DIR/stringtie_single/stringtieMerge_single/stringtieMerge_single.*.PJ.gtf \
$FILE_DIR/stringtie_merged_bam/*.PJ.MergedBam.gtf \
-p 10 -o $FILE_DIR/all_gtf_merged/all_gtf_merged.PJ.gtf

sleep 1

#both (84k)

qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N all_gtf_merge \
-e $FILE_DIR/all_gtf_merged/all_gtf.err -o $FILE_DIR/all_gtf_merged/all_gtf.out \
/usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/stringtie --merge -F 1.0 -T 1.0 \
$FILE_DIR/all_gtf_merged/all_gtf_merged.OJ.gtf \
$FILE_DIR/all_gtf_merged/all_gtf_merged.PJ.gtf \
-p 10 -o $FILE_DIR/all_gtf_merged/all_gtf_merged.OJPJ.gtf

#all (77k)

qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N all_gtf_merge \
-e $FILE_DIR/all_gtf_merged/all_gtf.err -o $FILE_DIR/all_gtf_merged/all_gtf.out \
/usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/stringtie --merge -F 1.0 -T 1.0 \
$FILE_DIR/stringtie_single/stringtieMerge_single/stringtieMerge_single.*.gtf \
$FILE_DIR/stringtie_merged_bam/*.MergedBam.gtf \
-p 10 -o $FILE_DIR/all_gtf_merged/all_gtf_merged.all.gtf
