#!/bin/bash
#@author : Lisa Rogers
#@created : 23-07-2018 @ 9:07 
#@modified : 10-09-18 : 11:09 
#Script to annotate the merged sorted bam files
#dir : /home/lisa/TropInst/P3_sRNA/featureCounts



featureCounts -p \
	-a /home/lisa/TropInst/P1_RNA_TB/mapping/ref/hg38.p12.ref.gtf \
	-o all.count \
	-F GTF  \
	-t exon \
	-g gene_name \
	-T 12 \
	--primary \
	-p \
	-s 0 \
	--extraAttributes gene_type \
	/home/lisa/TropInst/P3_sRNA/mapping/merge/*.bam

#Stranded -s1
featureCounts -p \
	-a /home/lisa/TropInst/P1_RNA_TB/mapping/ref/hg38.p12.ref.gtf \
	-o all.stranded.S1.count \
	-F GTF  \
	-t exon \
	-g gene_name \
	-T 12 \
	--primary \
	-p \
	-s 1 \
	--extraAttributes gene_type \
	/home/lisa/TropInst/P3_sRNA/mapping/merge/*.bam

#Stranded -s2
featureCounts -p \
	-a /home/lisa/TropInst/P1_RNA_TB/mapping/ref/hg38.p12.ref.gtf \
	-o all.stranded.S2.count \
	-F GTF  \
	-t exon \
	-g gene_name \
	-T 12 \
	--primary \
	-p \
	-s 2 \
	--extraAttributes gene_type \
	/home/lisa/TropInst/P3_sRNA/mapping/merge/*.bam



: <<'COMMENT'
	-a gtf  	# location of the GTF Annotation file
	-o NAME.count 	# name of the output file ; NAME.count.summary is also included
	-F GTF  	# Format of annotation file. GTF by default
	-t exon 	# Specific feature type in GTF file. exon by default
	-g gene_id 	# Specify the attribute type in GTF annotation.  gene_id by default
	-T 4		# no. of threads, max here 28
	-f		# counts reads at feature level rather # eg count reads at exons rather than genes
	--primary	# Counts primary alignments only (0X100 in SAM/BAM flag)
	-O		# Assigns reads to all overlapping meta-features 
	--minOverlap No	# Min number of overlapping bases in a read that is required for read assignment.  1 by default
	-s [012]	# Strand specific read counting. 0 - unstranded; 1 - stranded; 2- reversely stranded
	-p		# Fragments counted instead of reads, used for paired end reads
	-B		# only count read pairs that have both ends aligned
	inputFILE 	# a list of sam/bam file formats


for i in /home/lisa/TropInst/P1_RNA_TB/mapping/bam/*.merge.bam 
do 
echo $i
folder=$(basename "$i" ) # get folder name
ln -s $i $folder
done


COMMENT


