#!/bin/bash
#qc-check

fastqc --threads 8 --noextract --outdir /home/lisa/TropInst/P3_sRNA/QC/raw/lane1/ *.fastq


fastqc --threads 8 --noextract --outdir /home/lisa/TropInst/P3_sRNA/QC/raw/lane2/ *.fastq




fastqc --threads 12 --noextract --outdir /home/lisa/TropInst/P3_sRNA/QC/fastQC/trimmed/lane1/ *.fastq


fastqc --threads 12 --noextract --outdir /home/lisa/TropInst/P3_sRNA/QC/fastQC/trimmed/lane2/ *.fastq



cd /home/lisa/TropInst/P3_sRNA/reads/raw/lane1
fastqc --threads 12 --noextract --outdir /home/lisa/TropInst/P3_sRNA/QC/fastQC/full-trimmed/lane1/ *.trim.fasta 
cd /home/lisa/TropInst/P3_sRNA/reads/raw/lane2
fastqc --threads 12 --noextract --outdir /home/lisa/TropInst/P3_sRNA/QC/fastQC/full-trimmed/lane2/ *.trim.fasta 





fast#FastQC

##-- trimgalore


cd /home/lisa/TropInst/P3_sRNA/reads/trimmed/trimgalore/lane1
fastqc --threads 12 --noextract --outdir /home/lisa/TropInst/P3_sRNA/QC/fastQC/trimmed/trimgalore/lane1 *.fq
cd /home/lisa/TropInst/P3_sRNA/reads/trimmed/trimgalore/lane2
fastqc --threads 12 --noextract --outdir /home/lisa/TropInst/P3_sRNA/QC/fastQC/trimmed/trimgalore/lane2 *.fq




##-- cutadapt


cd /home/lisa/TropInst/P3_sRNA/reads/trimmed/cutadapt/lane1
fastqc --threads 12 --noextract --outdir /home/lisa/TropInst/P3_sRNA/QC/fastQC/trimmed/cutadapt/lane1 *.fasta
cd /home/lisa/TropInst/P3_sRNA/reads/trimmed/cutadapt/lane2
fastqc --threads 12 --noextract --outdir /home/lisa/TropInst/P3_sRNA/QC/fastQC/trimmed/cutadapt/lane2 *.fasta




