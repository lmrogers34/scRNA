#!/bin/bash
#@author : Lisa Rogers
#@modified : 16-07-2018 @ 09:27; 18-07-2018 @8:56
#2. Mapping
#Load genome
#STAR --runThreadN 14 \
#	--genomeDir /home/lisa/TropInst/P1_RNA_TB/mapping/index \
#	--genomeLoad LoadAndExit

# Mapping
#Load genome
STAR --runThreadN 4 \
	--genomeDir /home/lisa/TropInst/P1_RNA_TB/mapping/index \
	--genomeLoad LoadAndExit

for i in /home/lisa/TropInst/P3_sRNA/reads/trimmed/trimgalore/lane1/*_1_val_1.fq
do
echo $i


#R1 : rnauvb1_61_ba142q_0_S13_L004_R1_1.fq.gz
#R2 : rnauvb1_61_ba142q_0_S13_L004_R2_2.fq.gz
read1=$i
read2=${read1/_1_val_1.fq/_2_val_2.fq} #names the second read
echo "read1 is $read1"
echo "read2 is $read2"
outname=${read1/_1_val_1.fq/} #names the second read
outname=$(basename "$outname" ) # get folder name
echo "outname is $outname"

if [ -f /home/lisa/TropInst/P3_sRNA/mapping/trimgalore/lane1/bam/${outname}Aligned.sortedByCoord.out.bam ]
then
echo "STAR bam file exists for $outname"
else

echo "Beginning STAR for $outname"
STAR --runThreadN 4 \
	--genomeDir /home/lisa/TropInst/P1_RNA_TB/mapping/index \
	--genomeLoad LoadAndKeep \
	--readFilesIn $read1 $read2 \
	--outFileNamePrefix /home/lisa/TropInst/P3_sRNA/mapping/trimgalore/lane1/bam/${outname}. \
	--outSAMtype BAM SortedByCoordinate \
	--limitBAMsortRAM 10000000000


fi
done
echo "STAR done for all"
## End of storage; unload genome
## End of storage; unload genome
STAR --runThreadN 14 \
	--genomeDir /home/lisa/TropInst/P1_RNA_TB/mapping/index \
	--genomeLoad Remove 

cd /home/lisa/TropInst/P3_sRNA/mapping/trimgalore/lane2

# Mapping
#Load genome
STAR --runThreadN 4 \
	--genomeDir /home/lisa/TropInst/P1_RNA_TB/mapping/index \
	--genomeLoad LoadAndExit

for i in /home/lisa/TropInst/P3_sRNA/reads/trimmed/trimgalore/lane2/*_1_val_1.fq
do
echo $i


#R1 : rnauvb1_61_ba142q_0_S13_L004_R1_1.fq.gz
#R2 : rnauvb1_61_ba142q_0_S13_L004_R2_2.fq.gz
read1=$i
read2=${read1/_1_val_1.fq/_2_val_2.fq} #names the second read
echo "read1 is $read1"
echo "read2 is $read2"
outname=${read1/_1_val_1.fq/} #names the second read
outname=$(basename "$outname" ) # get folder name
echo "outname is $outname"

if [ -f /home/lisa/TropInst/P3_sRNA/mapping/trimgalore/lane2/bam/${outname}Aligned.sortedByCoord.out.bam ]
then
echo "STAR bam file exists for $outname"
else

echo "Beginning STAR for $outname"
STAR --runThreadN 4 \
	--genomeDir /home/lisa/TropInst/P1_RNA_TB/mapping/index \
	--genomeLoad LoadAndKeep \
	--readFilesIn $read1 $read2 \
	--outFileNamePrefix /home/lisa/TropInst/P3_sRNA/mapping/trimgalore/lane2/bam/${outname}. \
	--outSAMtype BAM SortedByCoordinate \
	--limitBAMsortRAM 10000000000


fi
done
echo "STAR done for all"

## End of storage; unload genome
STAR --runThreadN 14 \
	--genomeDir /home/lisa/TropInst/P1_RNA_TB/mapping/index \
	--genomeLoad Remove 


#Mapping commands
#	--genomeLoad LoadAndKeep \

#--runThreadN #no of threads used, in this case 14/28 #-- 1 cpu, 14 cores, each with 2 threads.  1x14x2=28
#--genomeDir \ #The dir where the indices should be store.  Should be made before STAR RUN.  Large files. Ref genome
#--readFilesIn $read1 $read2 \ # path to sequences to be mapped
#--readFilesCommand zcat \ # Used for when reads are compressed
#--genomeLoad LoandAndRemove \ #Used to load the Star Reference in shared memory until all the jobs are finished.
