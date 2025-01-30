#!bin/bash
#poojasingh
#comparing assembled gtfs to ref Onil and annotating, i compare OJ and PJ and both merged

# -r comparing with a reference annotation - will produce class codes: https://ccb.jhu.edu/software/stringtie/gffcompare.shtml
# -e Maximum distance (range) allowed from free ends of terminal exons of reference transcripts when assessing exon accuracy. By default, this is 100.
# -d Maximum distance (range) for grouping transcript start sites, by default 100.

FILE_DIR="/cl_tmp/singh_duenser/all_st26"


#OJ (69519 transcripts)


FILE_DIR="/cl_tmp/singh_duenser/all_st26"


#OJ 
qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N gffcmp.OJ \
-e $FILE_DIR/gffcompare/gffcmp.OJ.err -o $FILE_DIR/gffcompare/gffcmp.OJ.out \
/usr/people/EDVZ/duenser/tools/gffcompare-0.11.2.Linux_x86_64/gffcompare -r /cl_tmp/singh_duenser/reference/O_niloticus_UMD_NMBU.99.gff3 -e 100 -d 100 \
$FILE_DIR/all_gtf_merged/all_gtf_merged.OJ.gtf -o $FILE_DIR/gffcompare/gffcmp.OJ


#PJ 
qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N gffcmp.PJ \
-e $FILE_DIR/gffcompare/gffcmp.PJ.err -o $FILE_DIR/gffcompare/gffcmp.PJ.out \
/usr/people/EDVZ/duenser/tools/gffcompare-0.11.2.Linux_x86_64/gffcompare -r /cl_tmp/singh_duenser/reference/O_niloticus_UMD_NMBU.99.gff3 -e 100 -d 100 \
$FILE_DIR/all_gtf_merged/all_gtf_merged.PJ.gtf -o $FILE_DIR/gffcompare/gffcmp.PJ


#OJPJ 
qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N gffcmp.OJPJ \
-e $FILE_DIR/gffcompare/gffcmp.OJPJ.err -o $FILE_DIR/gffcompare/gffcmp.OJPJ.out \
/usr/people/EDVZ/duenser/tools/gffcompare-0.11.2.Linux_x86_64/gffcompare -r /cl_tmp/singh_duenser/reference/O_niloticus_UMD_NMBU.99.gff3 -e 100 -d 100 \
$FILE_DIR/all_gtf_merged/all_gtf_merged.OJPJ.gtf -o $FILE_DIR/gffcompare/gffcmp.OJPJ


#all

qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N gffcmp.all \
-e $FILE_DIR/gffcompare/gffcmp.all.ere -o $FILE_DIR/gffcompare/gffcmp.all.out \
/usr/people/EDVZ/duenser/tools/gffcompare-0.11.2.Linux_x86_64/gffcompare -r /cl_tmp/singh_duenser/reference/O_niloticus_UMD_NMBU.99.gff3 -e 100 -d 100 \
$FILE_DIR/all_gtf_merged/all_gtf_merged.all.gtf -o $FILE_DIR/gffcompare/gffcmp.all


done
