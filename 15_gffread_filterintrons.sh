#!bin/bash
#poojasingh
#intron filter: the longest intron in the O.niloticus reference is 200kb long so remove introns longer than that


for i in `ls /cl_tmp/singh_duenser/all_st26/gffcompare/gffread/filtered*.annotated.gtf | cut -d. -f 2- `;

do
echo $i
gffread /cl_tmp/singh_duenser/all_st26/gffcompare/gffread/filtered.$i -o /cl_tmp/singh_duenser/all_st26/gffcompare/gffread/filtered.gffread.$i -i 200000 -T 

done
