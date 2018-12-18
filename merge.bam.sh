#!/bin/bash
for i in /hdd/TropInst/P3_sRNA/mapping/trimgalore/lane1/bam/*.Aligned.sortedByCoord.out.bam
do

bam=$(basename $i ) # get folder name
lane2=/hdd/TropInst/P3_sRNA/mapping/trimgalore/lane2/bam
merge=/hdd/TropInst/P3_sRNA/mapping/trimgalore/merge

echo "bam is $bam"
echo "lane2 is $lane2"
echo "merge is $merge"
echo "$lane2/${bam}"
echo "${merge}/${bam}"



samtools merge -l 9 --output-fmt BAM --threads 12 $merge/${bam} $i ${lane2}/${bam} 


done



for i in /hdd/TropInst/P3_sRNA/mapping/trimgalore/lane1/bam/*.Aligned.sortedByCoord.out.bam
do

bam=$(basename $i ) # get folder name
lane2=/hdd/TropInst/P3_sRNA/mapping/trimgalore/lane2/bam
merge=/hdd/TropInst/P3_sRNA/mapping/trimgalore/merge

echo "bam is $bam"
echo "lane2 is $lane2"
echo "merge is $merge"
echo "$lane2/${bam}"
echo "${merge}/${bam}"



samtools merge -l 9 --output-fmt BAM --threads 12 $merge/${bam} $i ${lane2}/${bam} 


done
