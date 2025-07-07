#!bin/bash

## merge single .bam files for bio replicates per species per dev stage

FILE_DIR="/cl_tmp/singh_duenser/all_st26/star_assembly"
# define temp files
COUNTER=0
LAST=""
FULL=""

# loop whole folder and cut path
for i in `ls $FILE_DIR/*.bam | cut -d "/" -f 6`;
do

# for sanity check we need to remove numbering
SANITY="$(echo $i | cut -d'.' -f1-3)"

# counter is only needed for first loop, because $LAST is not defined
COUNTER=$[COUNTER + 1]

if [[ "$SANITY" != "$LAST" && "$COUNTER" != 1 ]]; 
then
echo "######"
echo $LAST
echo $FULL

qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N picard -e ./merged_bam_files/"picard."$LAST".err" -o ./merged_bam_files/"picard."$LAST".out" \
java -Dpicard.useLegacyParser=false -jar -Xmx2g /usr/people/EDVZ/duenser/tools/picard.jar MergeSamFiles \
$FULL \
-O ./merged_bam_files/$LAST".merged_files.bam" -MSD=true \
-VALIDATION_STRINGENCY=LENIENT

# reset string
FULL=""

# wait before submitting next job

sleep 10

fi

LAST=$SANITY
# build Input string
TMP="-I "$FILE_DIR/$i
FULL=$FULL" "$TMP

done

echo "#######"
echo $LAST
echo $FULL
echo "END"

# we need to submit the last loop result here

qsub -pe mpi 10 -l h_vmem=10G -b y -cwd -N picard -e ./merged_bam_files/"picard."$LAST".err" -o ./merged_bam_files/"picard."$LAST".out" \
java -Dpicard.useLegacyParser=false -jar /usr/people/EDVZ/duenser/tools/picard.jar MergeSamFiles \
$FULL \
-O ./merged_bam_files/$LAST".merged_files.bam" -MSD=true -VALIDATION_STRINGENCY=LENIENT
