#!/bin/bash
# To be run ~/TropInst/P3_sRNA/reads/raw/lane1 and ~/TropInst/P3_sRNA/reads/raw/lane2
#raw_fastq - contains untouvhed raw reads as given to me
# fastq - will contain syb links of original fastq reads
#Modified on 21st Aug 2018


#for i in /run/user/1001/gvfs/smb-share:server=172.24.134.212,share=mtb_tcells_ngs/rawdata/fastq/raw/KHeld_lane1/* # all the dirs
for i in /home/lisa/TropInst/P3_sRNA/nas/rawdata/fastq/raw/KHeld_lane1/* # all the dirs
do
if [ -d $i ]
then
echo $i
folder=$(basename "$i" ) # get folder name
cd $i # mv into the folder dir

echo "folder is $folder"
for j in *.fastq
do
echo $j

# sym link ln -s source link

ln -s ${i}/${j} /home/lisa/TropInst/P3_sRNA/reads/raw/lane1/${j}
done

cd ../ # mv out of the folder dir
fi
done
##Lane2

#for i in /run/user/1001/gvfs/smb-share:server=172.24.134.212,share=mtb_tcells_ngs/rawdata/fastq/raw/KHeld_lane2/* # all the dirs
for i in /home/lisa/TropInst/P3_sRNA/nas/rawdata/fastq/raw/KHeld_lane2/* # all the dirs
do
if [ -d $i ]
then
echo $i
folder=$(basename "$i" ) # get folder name
cd $i # mv into the folder dir

echo "folder is $folder"
for j in *.fastq
do
echo $j

# sym link ln -s source link

ln -s ${i}/${j} /home/lisa/TropInst/P3_sRNA/reads/raw/lane2/${j}
done

cd ../ # mv out of the folder dir
fi
done


#


#
#
