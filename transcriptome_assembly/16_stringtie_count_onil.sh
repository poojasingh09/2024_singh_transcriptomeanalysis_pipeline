# pooja.singh09@gmail.com
# stringtie get counts with final annotation
# -M to reduce multimapping, -e to only alow transcripts in the hopefully soon filtered annotation file, --rf strand info, -B ballgown files
#wget https://raw.githubusercontent.com/griffithlab/rnabio.org/master/assets/scripts/stringtie_expression_matrix.pl
#chmod +x stringtie_expression_matrix.pl


DIR_IN="/cl_tmp/singh_duenser/all_st26/star_assembly"
DIR_OUT="/cl_tmp/singh_duenser/all_st26/count_matrices/"

'''
#for i in `ls $DIR_IN/*.bam | cut -f 6 -d"/" | cut -f -4 -d"."`;

for i in `ls $DIR_IN/On*.bam | cut -f 6 -d"/" | cut -f 1 -d"."`;

do

mkdir $DIR_OUT/$i


qsub -pe mpi 2 -l h_vmem=5G -cwd -V -N $i"_stringtie_out_onil" \
-o $DIR_OUT/$i"_stringtie_count_onil.out" -e $DIR_OUT/$i"_stringtie_count_onil.err" -b y \
/usr/people/EDVZ/duenser/tools/stringtie-2.0.6.Linux_x86_64/stringtie -p 2 -M 0.0 -e -B -G \
/cl_tmp/singh_duenser/all_st26/gffcompare/gffread/filtered.gffread.gffcmp.all.annotated.gtf \
--rf  -o $DIR_OUT/$i/transcripts_onil.gtf \
-A $DIR_OUT/$i/gene_abundances_onil.tsv \
$DIR_IN/$i.bam




done
'''

cd $DIR_OUT

/cl_tmp/singh_duenser/all_st26/count_matrices/stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='On20d3,Aa.L.OJ.1,Aa.L.OJ.2,Aa.L.OJ.3,Aa.L.OJ.4,Aa.L.OJ.5,Aa.L.PJ.1,Aa.L.PJ.2,Aa.L.PJ.3,Aa.L.PJ.4,Aa.L.PJ.5,Ab.L.OJ.1,Ab.L.OJ.2,Ab.L.OJ.3,Ab.L.OJ.4,Ab.L.OJ.5,Ab.L.PJ.1,Ab.L.PJ.2,Ab.L.PJ.3,Ab.L.PJ.4,Ab.L.PJ.5,Ah.L.OJ.1,Ah.L.OJ.2,Ah.L.OJ.3,Ah.L.OJ.4,Ah.L.OJ.5,Ah.L.PJ.1,Ah.L.PJ.2,Ah.L.PJ.3,Ah.L.PJ.4,Ah.L.PJ.5,Ch.L.OJ.1,Ch.L.OJ.2,Ch.L.OJ.3,Ch.L.OJ.4,Ch.L.OJ.5,Ch.L.PJ.1,Ch.L.PJ.2,Ch.L.PJ.3,Ch.L.PJ.4,Ch.L.PJ.5,Gh.L.OJ.1,Gh.L.OJ.2,Gh.L.OJ.3,Gh.L.OJ.4,Gh.L.OJ.5,Gh.L.PJ.1,Gh.L.PJ.2,Gh.L.PJ.3,Gh.L.PJ.4,Gh.L.PJ.5,Gp.L.OJ.1,Gp.L.OJ.2,Gp.L.OJ.3,Gp.L.OJ.4,Gp.L.OJ.5,Gp.L.PJ.1,Gp.L.PJ.2,Gp.L.PJ.3,Gp.L.PJ.4,Gp.L.PJ.5,Ht.L.OJ.1,Ht.L.OJ.2,Ht.L.OJ.3,Ht.L.OJ.4,Ht.L.OJ.5,Ht.L.PJ.1,Ht.L.PJ.2,Ht.L.PJ.3,Ht.L.PJ.4,Ht.L.PJ.5,Lt.L.OJ.1,Lt.L.OJ.2,Lt.L.OJ.3,Lt.L.OJ.4,Lt.L.OJ.5,Lt.L.PJ.1,Lt.L.PJ.2,Lt.L.PJ.3,Lt.L.PJ.4,Lt.L.PJ.5,Ml.L.OJ.1,Ml.L.OJ.2,Ml.L.OJ.3,Ml.L.OJ.4,Ml.L.OJ.5,Ml.L.PJ.1,Ml.L.PJ.2,Ml.L.PJ.3,Ml.L.PJ.4,Ml.L.PJ.5,No.L.OJ.1,No.L.OJ.2,No.L.OJ.3,No.L.OJ.4,No.L.OJ.5,No.L.PJ.1,No.L.PJ.2,No.L.PJ.3,No.L.PJ.4,No.L.PJ.5,Pf.L.OJ.1,Pf.L.OJ.2,Pf.L.OJ.3,Pf.L.OJ.4,Pf.L.OJ.5,Pf.L.PJ.1,Pf.L.PJ.2,Pf.L.PJ.3,Pf.L.PJ.4,Pf.L.PJ.5,Pp.L.OJ.1,Pp.L.OJ.2,Pp.L.OJ.3,Pp.L.OJ.4,Pp.L.OJ.5,Pp.L.PJ.1,Pp.L.PJ.2,Pp.L.PJ.3,Pp.L.PJ.4,Pp.L.PJ.5,Prp.L.OJ.1,Prp.L.OJ.2,Prp.L.OJ.3,Prp.L.OJ.4,Prp.L.OJ.5,Prp.L.PJ.1,Prp.L.PJ.2,Prp.L.PJ.3,Prp.L.PJ.4,Prp.L.PJ.5,Ps.L.OJ.1,Ps.L.OJ.2,Ps.L.OJ.3,Ps.L.OJ.4,Ps.L.OJ.5,Ps.L.PJ.1,Ps.L.PJ.2,Ps.L.PJ.3,Ps.L.PJ.4,Ps.L.PJ.5,Ptb.L.OJ.1,Ptb.L.OJ.2,Ptb.L.OJ.3,Ptb.L.OJ.4,Ptb.L.OJ.5,Ptb.L.PJ.1,Ptb.L.PJ.2,Ptb.L.PJ.3,Ptb.L.PJ.4,Ptb.L.PJ.5,Py.L.OJ.1,Py.L.OJ.2,Py.L.OJ.3,Py.L.OJ.4,Py.L.OJ.5,Py.L.PJ.1,Py.L.PJ.2,Py.L.PJ.3,Py.L.PJ.4,Py.L.PJ.5,Sd.L.OJ.1,Sd.L.OJ.2,Sd.L.OJ.3,Sd.L.OJ.4,Sd.L.OJ.5,Sd.L.PJ.1,Sd.L.PJ.2,Sd.L.PJ.3,Sd.L.PJ.4,Sd.L.PJ.5,Sf.L.OJ.1,Sf.L.OJ.2,Sf.L.OJ.3,Sf.L.OJ.4,Sf.L.OJ.5,Sf.L.PJ.1,Sf.L.PJ.2,Sf.L.PJ.3,Sf.L.PJ.4,Sf.L.PJ.5,Tm.L.OJ.1,Tm.L.OJ.2,Tm.L.OJ.3,Tm.L.OJ.4,Tm.L.OJ.5,Tm.L.PJ.1,Tm.L.PJ.2,Tm.L.PJ.3,Tm.L.PJ.4,Tm.L.PJ.5,Tt.L.OJ.1,Tt.L.OJ.2,Tt.L.OJ.3,Tt.L.OJ.4,Tt.L.OJ.5,Tt.L.PJ.1,Tt.L.PJ.2,Tt.L.PJ.3,Tt.L.PJ.4,Tt.L.PJ.5' --transcript_matrix_file=transcript_tpm_all_samples_onil.tsv --gene_matrix_file=gene_tpm_all_samples_onil.tsv &

/cl_tmp/singh_duenser/all_st26/count_matrices/stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs='On20d3,Aa.L.OJ.1,Aa.L.OJ.2,Aa.L.OJ.3,Aa.L.OJ.4,Aa.L.OJ.5,Aa.L.PJ.1,Aa.L.PJ.2,Aa.L.PJ.3,Aa.L.PJ.4,Aa.L.PJ.5,Ab.L.OJ.1,Ab.L.OJ.2,Ab.L.OJ.3,Ab.L.OJ.4,Ab.L.OJ.5,Ab.L.PJ.1,Ab.L.PJ.2,Ab.L.PJ.3,Ab.L.PJ.4,Ab.L.PJ.5,Ah.L.OJ.1,Ah.L.OJ.2,Ah.L.OJ.3,Ah.L.OJ.4,Ah.L.OJ.5,Ah.L.PJ.1,Ah.L.PJ.2,Ah.L.PJ.3,Ah.L.PJ.4,Ah.L.PJ.5,Ch.L.OJ.1,Ch.L.OJ.2,Ch.L.OJ.3,Ch.L.OJ.4,Ch.L.OJ.5,Ch.L.PJ.1,Ch.L.PJ.2,Ch.L.PJ.3,Ch.L.PJ.4,Ch.L.PJ.5,Gh.L.OJ.1,Gh.L.OJ.2,Gh.L.OJ.3,Gh.L.OJ.4,Gh.L.OJ.5,Gh.L.PJ.1,Gh.L.PJ.2,Gh.L.PJ.3,Gh.L.PJ.4,Gh.L.PJ.5,Gp.L.OJ.1,Gp.L.OJ.2,Gp.L.OJ.3,Gp.L.OJ.4,Gp.L.OJ.5,Gp.L.PJ.1,Gp.L.PJ.2,Gp.L.PJ.3,Gp.L.PJ.4,Gp.L.PJ.5,Ht.L.OJ.1,Ht.L.OJ.2,Ht.L.OJ.3,Ht.L.OJ.4,Ht.L.OJ.5,Ht.L.PJ.1,Ht.L.PJ.2,Ht.L.PJ.3,Ht.L.PJ.4,Ht.L.PJ.5,Lt.L.OJ.1,Lt.L.OJ.2,Lt.L.OJ.3,Lt.L.OJ.4,Lt.L.OJ.5,Lt.L.PJ.1,Lt.L.PJ.2,Lt.L.PJ.3,Lt.L.PJ.4,Lt.L.PJ.5,Ml.L.OJ.1,Ml.L.OJ.2,Ml.L.OJ.3,Ml.L.OJ.4,Ml.L.OJ.5,Ml.L.PJ.1,Ml.L.PJ.2,Ml.L.PJ.3,Ml.L.PJ.4,Ml.L.PJ.5,No.L.OJ.1,No.L.OJ.2,No.L.OJ.3,No.L.OJ.4,No.L.OJ.5,No.L.PJ.1,No.L.PJ.2,No.L.PJ.3,No.L.PJ.4,No.L.PJ.5,Pf.L.OJ.1,Pf.L.OJ.2,Pf.L.OJ.3,Pf.L.OJ.4,Pf.L.OJ.5,Pf.L.PJ.1,Pf.L.PJ.2,Pf.L.PJ.3,Pf.L.PJ.4,Pf.L.PJ.5,Pp.L.OJ.1,Pp.L.OJ.2,Pp.L.OJ.3,Pp.L.OJ.4,Pp.L.OJ.5,Pp.L.PJ.1,Pp.L.PJ.2,Pp.L.PJ.3,Pp.L.PJ.4,Pp.L.PJ.5,Prp.L.OJ.1,Prp.L.OJ.2,Prp.L.OJ.3,Prp.L.OJ.4,Prp.L.OJ.5,Prp.L.PJ.1,Prp.L.PJ.2,Prp.L.PJ.3,Prp.L.PJ.4,Prp.L.PJ.5,Ps.L.OJ.1,Ps.L.OJ.2,Ps.L.OJ.3,Ps.L.OJ.4,Ps.L.OJ.5,Ps.L.PJ.1,Ps.L.PJ.2,Ps.L.PJ.3,Ps.L.PJ.4,Ps.L.PJ.5,Ptb.L.OJ.1,Ptb.L.OJ.2,Ptb.L.OJ.3,Ptb.L.OJ.4,Ptb.L.OJ.5,Ptb.L.PJ.1,Ptb.L.PJ.2,Ptb.L.PJ.3,Ptb.L.PJ.4,Ptb.L.PJ.5,Py.L.OJ.1,Py.L.OJ.2,Py.L.OJ.3,Py.L.OJ.4,Py.L.OJ.5,Py.L.PJ.1,Py.L.PJ.2,Py.L.PJ.3,Py.L.PJ.4,Py.L.PJ.5,Sd.L.OJ.1,Sd.L.OJ.2,Sd.L.OJ.3,Sd.L.OJ.4,Sd.L.OJ.5,Sd.L.PJ.1,Sd.L.PJ.2,Sd.L.PJ.3,Sd.L.PJ.4,Sd.L.PJ.5,Sf.L.OJ.1,Sf.L.OJ.2,Sf.L.OJ.3,Sf.L.OJ.4,Sf.L.OJ.5,Sf.L.PJ.1,Sf.L.PJ.2,Sf.L.PJ.3,Sf.L.PJ.4,Sf.L.PJ.5,Tm.L.OJ.1,Tm.L.OJ.2,Tm.L.OJ.3,Tm.L.OJ.4,Tm.L.OJ.5,Tm.L.PJ.1,Tm.L.PJ.2,Tm.L.PJ.3,Tm.L.PJ.4,Tm.L.PJ.5,Tt.L.OJ.1,Tt.L.OJ.2,Tt.L.OJ.3,Tt.L.OJ.4,Tt.L.OJ.5,Tt.L.PJ.1,Tt.L.PJ.2,Tt.L.PJ.3,Tt.L.PJ.4,Tt.L.PJ.5' --transcript_matrix_file=transcript_fpkm_all_samples_onil.tsv --gene_matrix_file=gene_fpkm_all_samples_onil.tsv &

/cl_tmp/singh_duenser/all_st26/count_matrices/stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs='On20d3,Aa.L.OJ.1,Aa.L.OJ.2,Aa.L.OJ.3,Aa.L.OJ.4,Aa.L.OJ.5,Aa.L.PJ.1,Aa.L.PJ.2,Aa.L.PJ.3,Aa.L.PJ.4,Aa.L.PJ.5,Ab.L.OJ.1,Ab.L.OJ.2,Ab.L.OJ.3,Ab.L.OJ.4,Ab.L.OJ.5,Ab.L.PJ.1,Ab.L.PJ.2,Ab.L.PJ.3,Ab.L.PJ.4,Ab.L.PJ.5,Ah.L.OJ.1,Ah.L.OJ.2,Ah.L.OJ.3,Ah.L.OJ.4,Ah.L.OJ.5,Ah.L.PJ.1,Ah.L.PJ.2,Ah.L.PJ.3,Ah.L.PJ.4,Ah.L.PJ.5,Ch.L.OJ.1,Ch.L.OJ.2,Ch.L.OJ.3,Ch.L.OJ.4,Ch.L.OJ.5,Ch.L.PJ.1,Ch.L.PJ.2,Ch.L.PJ.3,Ch.L.PJ.4,Ch.L.PJ.5,Gh.L.OJ.1,Gh.L.OJ.2,Gh.L.OJ.3,Gh.L.OJ.4,Gh.L.OJ.5,Gh.L.PJ.1,Gh.L.PJ.2,Gh.L.PJ.3,Gh.L.PJ.4,Gh.L.PJ.5,Gp.L.OJ.1,Gp.L.OJ.2,Gp.L.OJ.3,Gp.L.OJ.4,Gp.L.OJ.5,Gp.L.PJ.1,Gp.L.PJ.2,Gp.L.PJ.3,Gp.L.PJ.4,Gp.L.PJ.5,Ht.L.OJ.1,Ht.L.OJ.2,Ht.L.OJ.3,Ht.L.OJ.4,Ht.L.OJ.5,Ht.L.PJ.1,Ht.L.PJ.2,Ht.L.PJ.3,Ht.L.PJ.4,Ht.L.PJ.5,Lt.L.OJ.1,Lt.L.OJ.2,Lt.L.OJ.3,Lt.L.OJ.4,Lt.L.OJ.5,Lt.L.PJ.1,Lt.L.PJ.2,Lt.L.PJ.3,Lt.L.PJ.4,Lt.L.PJ.5,Ml.L.OJ.1,Ml.L.OJ.2,Ml.L.OJ.3,Ml.L.OJ.4,Ml.L.OJ.5,Ml.L.PJ.1,Ml.L.PJ.2,Ml.L.PJ.3,Ml.L.PJ.4,Ml.L.PJ.5,No.L.OJ.1,No.L.OJ.2,No.L.OJ.3,No.L.OJ.4,No.L.OJ.5,No.L.PJ.1,No.L.PJ.2,No.L.PJ.3,No.L.PJ.4,No.L.PJ.5,Pf.L.OJ.1,Pf.L.OJ.2,Pf.L.OJ.3,Pf.L.OJ.4,Pf.L.OJ.5,Pf.L.PJ.1,Pf.L.PJ.2,Pf.L.PJ.3,Pf.L.PJ.4,Pf.L.PJ.5,Pp.L.OJ.1,Pp.L.OJ.2,Pp.L.OJ.3,Pp.L.OJ.4,Pp.L.OJ.5,Pp.L.PJ.1,Pp.L.PJ.2,Pp.L.PJ.3,Pp.L.PJ.4,Pp.L.PJ.5,Prp.L.OJ.1,Prp.L.OJ.2,Prp.L.OJ.3,Prp.L.OJ.4,Prp.L.OJ.5,Prp.L.PJ.1,Prp.L.PJ.2,Prp.L.PJ.3,Prp.L.PJ.4,Prp.L.PJ.5,Ps.L.OJ.1,Ps.L.OJ.2,Ps.L.OJ.3,Ps.L.OJ.4,Ps.L.OJ.5,Ps.L.PJ.1,Ps.L.PJ.2,Ps.L.PJ.3,Ps.L.PJ.4,Ps.L.PJ.5,Ptb.L.OJ.1,Ptb.L.OJ.2,Ptb.L.OJ.3,Ptb.L.OJ.4,Ptb.L.OJ.5,Ptb.L.PJ.1,Ptb.L.PJ.2,Ptb.L.PJ.3,Ptb.L.PJ.4,Ptb.L.PJ.5,Py.L.OJ.1,Py.L.OJ.2,Py.L.OJ.3,Py.L.OJ.4,Py.L.OJ.5,Py.L.PJ.1,Py.L.PJ.2,Py.L.PJ.3,Py.L.PJ.4,Py.L.PJ.5,Sd.L.OJ.1,Sd.L.OJ.2,Sd.L.OJ.3,Sd.L.OJ.4,Sd.L.OJ.5,Sd.L.PJ.1,Sd.L.PJ.2,Sd.L.PJ.3,Sd.L.PJ.4,Sd.L.PJ.5,Sf.L.OJ.1,Sf.L.OJ.2,Sf.L.OJ.3,Sf.L.OJ.4,Sf.L.OJ.5,Sf.L.PJ.1,Sf.L.PJ.2,Sf.L.PJ.3,Sf.L.PJ.4,Sf.L.PJ.5,Tm.L.OJ.1,Tm.L.OJ.2,Tm.L.OJ.3,Tm.L.OJ.4,Tm.L.OJ.5,Tm.L.PJ.1,Tm.L.PJ.2,Tm.L.PJ.3,Tm.L.PJ.4,Tm.L.PJ.5,Tt.L.OJ.1,Tt.L.OJ.2,Tt.L.OJ.3,Tt.L.OJ.4,Tt.L.OJ.5,Tt.L.PJ.1,Tt.L.PJ.2,Tt.L.PJ.3,Tt.L.PJ.4,Tt.L.PJ.5' --transcript_matrix_file=transcript_coverage_all_samples_onil.tsv --gene_matrix_file=gene_coverage_all_samples_onil.tsv &
