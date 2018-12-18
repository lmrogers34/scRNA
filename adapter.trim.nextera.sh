#!/bin/bash
for i in *_1.fastq
do
read1=$(basename "$i" ) # get folder name
read2=${read1/_1/_2} #names the second read


echo "read 1 is $read1"
echo "read 2 is $read2"


trim_galore --quality 0 --paired --nextera $read1 $read2


done
