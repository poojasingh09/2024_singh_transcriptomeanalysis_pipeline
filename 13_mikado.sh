#!bin/bash
#poojasingh
#mikado: get gtf statistics
# i will run this for each species (for tables) as well as the merged OJPJ gtf
# this is necessary to achieve an understanding of the exons per transcripts and intron lengths
#this was run locally because the server hates mikado




for i in `ls *.gtf`;

do

mikado util stats  $i > $i".stats"


done




echo "sample" > sample.txt
echo "genes_total" > genes_total.txt
echo "transcripts_total" > transcripts_total.txt
echo "transcripts_avrg" > transcripts_avrg.txt
echo "monoexonic_transcripts" > transcripts_monoexonic.txt
echo "intron_length95" > intron_length95.txt

for i in `ls *.stats | cut -f -2 -d"."`;

do

echo $i >> sample.txt
grep -w "Number of genes" $i".gtf.stats" | grep -v "coding" | cut -f 2 >> genes_total.txt
grep -w "Transcripts per gene" $i".gtf.stats" | cut -f 2 >> transcripts_total.txt
grep -w "Transcripts per gene" $i".gtf.stats" | cut -f 3 >> transcripts_avrg.txt
grep -w "Monoexonic transcripts" $i".gtf.stats" | cut -f 2 >> transcripts_monoexonic.txt
grep -w "Intron lengths" $i".gtf.stats" | grep -v CDS | grep -v mRNA | cut -f 13 >> intron_length95.txt

done

paste sample.txt genes_total.txt transcripts_total.txt transcripts_avrg.txt transcripts_monoexonic.txt intron_length95.txt | sed 's/,//g' > stats.all
