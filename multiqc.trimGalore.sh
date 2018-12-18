#!/bin/bash
#
##Lane 1
cd /home/lisa/TropInst/P3_sRNA/QC/fastQC/trimmed/trimgalore/lane1
## ALL
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_all_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/All ./*
## Read 1
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_R1_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/read1 ./*_1_*
## Read 2
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_R2_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/read2 ./*_2_*
#A
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_A_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/A ./A*
#B
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_B_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/B ./B*
#C
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_C_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/C ./C*
#D
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_D_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/D ./D*
#E
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_E_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/E ./E*
#F
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_F_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/F ./F*
#G
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_G_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/G ./G*
#H
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/Lane1_H_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane1/H ./H*


##Lane 2

cd /home/lisa/TropInst/P3_sRNA/QC/fastQC/trimmed/trimgalore/lane2
## ALL
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_all_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/All ./*
## Read 1
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_R1_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/read1 ./*_1_*
## Read 2
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_R2_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/read2 ./*_2_*
#A
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_A_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/A ./A*
#B
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_B_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/B ./B*
#C
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_C_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/C ./C*
#D
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_D_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/D ./D*
#E
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_E_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/E ./E*
#F
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_F_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/F ./F*
#G
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_G_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/G ./G*
#H
multiqc --filename /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/Lane2_H_multiqc_report.html --outdir /home/lisa/TropInst/P3_sRNA/QC/multiqc/trimmed/trimgalore/lane2/H ./H*


