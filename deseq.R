#!/usr/bin/Rscript
## @author : Lisa Rogers
## @date : 01-08-2018
## @ purpose : RNA-seq analysis with DeSeq2

### Step 0 : Call Libraries

library(DESeq2)
library(ggplot2)
#install.packages('Seurat')
# On Command line
#sudo apt-get install libhdf5-dev
library(Seurat)


# 1, Read in input DATA
# from feature countdata, command line : featureCounts -p -a /home/lisa/TropInst/P1_RNA_TB/mapping/ref/hg38.p12.ref.gtf -o all.count -F GTF -t exon -g gene_name -T 12, --primary -p -s 0 --extraAttributes gene_type *.merge.bam

countdata=read.csv("all.stranded.count", sep="", head=T, skip=1, row.names = "Geneid") # ignores first line
rawdata=countdata
#ignores first 6 columns (chr, start, end, strand, length, gene_type)
colnames(countdata)[7:67]
countdata <- countdata[ ,7:ncol(countdata)]

# Remove .merge.bam from file names
colnames(countdata) <- gsub("\\.merge.bam$", "", colnames(countdata))

colnames(countdata) <- gsub("X.hdd.TropInst.P1_RNA_TB.mapping.bam_2.", "", colnames(countdata))
colnames(countdata) <- gsub("X.hdd.TropInst.P1_RNA_TB.mapping.bam_12.", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)



 [1] "rnauvb1_102_ta168r_2_S6"   "rnauvb1_106_uc174n_2_S10" 
 [3] "rnauvb1_10_km506x_2_S10"   "rnauvb1_110_uc194f_2_S14" 
 [5] "rnauvb1_114_mb151n_2_S18"  "rnauvb1_118_mb158n_2_S21" 
 [7] "rnauvb1_14_km546j_2_S14"   "rnauvb1_18_km523j_2_S2"   
 [9] "rnauvb1_22_km552x_2_S6"    "rnauvb1_26_km540u_2_S2"   
[11] "rnauvb1_2_km530d_2_S2"     "rnauvb1_30_km513t_2_S14"  
[13] "rnauvb1_34_ba179f_2_S18"   "rnauvb1_38_ba119q_2_S22"  
[15] "rnauvb1_42_ba163b_2_S18"   "rnauvb1_46_ba173s_2_S22"  
[17] "rnauvb1_50_ba127s_2_S10"   "rnauvb1_54_ba181u_2_S6"   
[19] "rnauvb1_58_ba109z_2_S10"   "rnauvb1_62_ba142q_2_S14"  
[21] "rnauvb1_66_ba162v_2_S18"   "rnauvb1_6_km526c_2_S6"    
[23] "rnauvb1_70_ba135u_2_S22"   "rnauvb1_74_ba104s_2_S2"   
[25] "rnauvb1_78_km504k_2_S6"    "rnauvb1_82_ta147d_2_S10"  
[27] "rnauvb1_86_uc121e_2_S14"   "rnauvb1_90_uc189y_2_S18"  
[29] "rnauvb1_94_wt115x_2_S22"   "rnauvb1_98_ta137p_2_S2"   
[31] "rnauvb1_103_ta168r_12_S7"  "rnauvb1_107_uc174n_12_S11"
[33] "rnauvb1_111_uc194f_12_S15" "rnauvb1_115_mb151n_12_S19"
[35] "rnauvb1_119_mb158n_12_S22" "rnauvb1_11_km506x_12_S11" 
[37] "rnauvb1_15_km546j_12_S15"  "rnauvb1_19_km523j_12_S3"  
[39] "rnauvb1_23_km552x_12_SM"   "rnauvb1_27_km540u_12_S3"  
[41] "rnauvb1_31_km513t_12_S15"  "rnauvb1_35_ba179f_12_S19" 
[43] "rnauvb1_39_ba119q_12_S23"  "rnauvb1_3_km530d_12_S3"   
[45] "rnauvb1_43_ba163b_12_S19"  "rnauvb1_47_ba173s_12_S23" 
[47] "rnauvb1_51_ba127s_12_S11"  "rnauvb1_55_ba181u_12_S7"  
[49] "rnauvb1_59_ba109z_12_S11"  "rnauvb1_63_ba142q_12_S15" 
[51] "rnauvb1_67_ba162v_12_S19"  "rnauvb1_71_ba135u_12_S23" 
[53] "rnauvb1_79_km504k_12_S7"   "rnauvb1_7_km526c_12_S7"   
[55] "rnauvb1_83_ta147d_12_S11"  "rnauvb1_87_uc121e_12_S15" 
[57] "rnauvb1_91_uc189y_12_S19"  "rnauvb1_95_wt115x_12_S23" 
[59] "rnauvb1_99_ta137p_12_S3



# Assign conditions (timepoint, SampleID)
(timepoint <- c("26", "0", "26", "0", "26", "0", "26", "0", "26", "26","0", "26", "0", "26", "0", "0", "26", "0", "26", "0", "26", "0", "26", "0", "26", "0", "26", "0", "26", "0", "26", "0", "26", "26", "0", "26", "0", "0", "26", "0", "26", "0", "26", "0", "26", "0", "26", "0", "26", "0", "26", "0", "26", "0", "26", "26", "0", "26", "0", "0"))

(SampleID <- c("ta137p", "ta168r", "ta168r", "uc174n", "uc174n", "uc194f", "uc194f", "mb151n", "mb151n", "mb158n","mb108h", "km506x", "km546j", "km546j", "km523j", "km530d", "km523j", "km552x", "km552x", "km540u", "km540u", "km513t", "km513t", "ba179f", "ba179f", "ba119q", "ba119q", "ba163b", "ba163b", "ba173s", "ba173s", "ba127s", "km530d", "ba127s", "ba181u", "ba181u", "ba109z", "km526c", "ba109z", "ba142q", "ba142q", "ba162v", "ba162v", "ba135u", "ba135u", "ba104s", "ba104s", "km504k", "km504k", "ta147d", "ta147d", "uc121e", "uc121e", "uc189y", "km526c", "uc189y", "wt115x", "wt115x", "ta137p", "km506x" ))


(date <- c("28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "04.09.2017", "04.09.2017", "04.09.2017", "04.09.2017", "04.07.2017", "04.07.2017", "04.07.2017", "21.07.2017", "04.07.2017", "21.07.2017R", "21.07.2017R", "21.07.2017R", "21.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "04.07.2017", "04.07.2017", "04.07.2017", "04.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "04.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "04.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "21.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "28.07.2017", "04.07.2017"))

## (placeholder <- c("rnauvb1_100_ta137p_26_S4", "rnauvb1_101_ta168r_0_S5", "rnauvb1_104_ta168r_26_S8", "rnauvb1_105_uc174n_0_S9","rnauvb1_108_uc174n_26_S12", "rnauvb1_109_uc194f_0_S13", "rnauvb1_112_uc194f_26_S16", "rnauvb1_113_mb151n_0_S17", "rnauvb1_116_mb151n_26_S20", "rnauvb1_120_mb158n_26_S23", "rnauvb1_121_mb108h_0_S24", "rnauvb1_12_km506x_26_S12", "rnauvb1_13_km546j_0_S13", "rnauvb1_16_km546j_26_S16", "rnauvb1_17_km523j_0_S1", "rnauvb1_1_km530d_0_S1", "rnauvb1_20_km523j_26_SM", "rnauvb1_21_km552x_0_SM", "rnauvb1_24_km552x_26_SM", "rnauvb1_25_km540u_0_S1", "rnauvb1_28_km540u_26_S12", "rnauvb1_29_km513t_0_S13", "rnauvb1_32_km513t_26_S16", "rnauvb1_33_ba179f_0_S17", "rnauvb1_36_ba179f_26_S20", "rnauvb1_37_ba119q_0_S21", "rnauvb1_40_ba119q_26_S24", "rnauvb1_41_ba163b_0_S17", "rnauvb1_44_ba163b_26_S20", "rnauvb1_45_ba173s_0_S21", "rnauvb1_48_ba173s_26_S24", "rnauvb1_49_ba127s_0_S9", "rnauvb1_4_km530d_26_S4", "rnauvb1_52_ba127s_26_S4", "rnauvb1_53_ba181u_0_S5", "rnauvb1_56_ba181u_26_S8", "rnauvb1_57_ba109z_0_S9", "rnauvb1_5_km526c_0_S5", "rnauvb1_60_ba109z_26_S12", "rnauvb1_61_ba142q_0_S13", "rnauvb1_64_ba142q_26_S16", "rnauvb1_65_ba162v_0_S17", "rnauvb1_68_ba162v_26_S20", "rnauvb1_69_ba135u_0_S21", "rnauvb1_72_ba135u_26_S24", "rnauvb1_73_ba104s_0_S1", "rnauvb1_76_ba104s_26_S4", "rnauvb1_77_km504k_0_S5", "rnauvb1_80_km504k_26_S8", "rnauvb1_81_ta147d_0_S9", "rnauvb1_84_ta147d_26_S12", "rnauvb1_85_uc121e_0_S13", "rnauvb1_88_uc121e_26_S16", "rnauvb1_89_uc189y_0_S17", "rnauvb1_8_km526c_26_S8", "rnauvb1_92_uc189y_26_S20", "rnauvb1_93_wt115x_0_S21", "rnauvb1_96_wt115x_26_S24", "rnauvb1_97_ta137p_0_S1", "rnauvb1_9_km506x_0_S9"))

#(prep <- c("lib_prep_3_green", "lib_prep_3_green", "lib_prep_3_green", "lib_prep_3_whiteStar","lib_prep_3_whiteStar", "lib_prep_3_whiteStar", "lib_prep_3_whiteStar", "lib_prep_4", "lib_prep_4", "lib_prep_4", "lib_prep_4", "lib_prep_1", "lib_prep_1", "lib_prep_1", "lib_prep_2_white", "lib_prep_1", "lib_prep_2R", "lib_prep_2R", "lib_prep_2R", "lib_prep_2_green", "lib_prep_2_green", "lib_prep_2_green", "lib_prep_2_green", "lib_prep_1", "lib_prep_1", "lib_prep_1", "lib_prep_1", "lib_prep_2_yellow", "lib_prep_2_yellow", "lib_prep_2_yellow", "lib_prep_2_yellow", "lib_prep_2_blue", "lib_prep_1", "lib_prep_2_blue", "lib_prep_2_blue", "lib_prep_2_blue", "lib_prep_2_whiteStar", "lib_prep_1", "lib_prep_2_whiteStar", "lib_prep_2_whiteStar", "lib_prep_2_whiteStar", "lib_prep_2_greenStar", "lib_prep_2_greenStar", "lib_prep_2_greenStar", "lib_prep_2_greenStar", "lib_prep_3_white", "lib_prep_3_white", "lib_prep_3_white", "lib_prep_3_white", "lib_prep_3_yellow", "lib_prep_3_yellow", "lib_prep_3_yellow", "lib_prep_3_yellow", "lib_prep_3_blue", "lib_prep_1", "lib_prep_3_blue", "lib_prep_3_blue", "lib_prep_3_blue", "lib_prep_3_blue", "lib_prep_1"))


(prep <- c("lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_3","lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_4", "lib_prep_4", "lib_prep_4", "lib_prep_4", "lib_prep_1", "lib_prep_1", "lib_prep_1", "lib_prep_2", "lib_prep_1", "lib_prep_2R", "lib_prep_2R", "lib_prep_2R", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_1", "lib_prep_1", "lib_prep_1", "lib_prep_1", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_1", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_1", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_2", "lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_1", "lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_3", "lib_prep_1"))


# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), timepoint, SampleID, date, prep))

#Subset Data

#Subset Data
countdata <- subset(countdata, select=-c(rnauvb1_37_ba119q_0_S21, rnauvb1_36_ba179f_26_S20))
coldata <- subset(coldata, !(rownames(coldata) %in% c("rnauvb1_37_ba119q_0_S21", "rnauvb1_36_ba179f_26_S20")))



#countdata <- subset(countdata, select=c(rnauvb1_101_ta168r_0_S5, rnauvb1_104_ta168r_26_S8, rnauvb1_109_uc194f_0_S13, rnauvb1_112_uc194f_26_S16,  rnauvb1_65_ba162v_0_S17, rnauvb1_68_ba162v_26_S20, rnauvb1_69_ba135u_0_S21, rnauvb1_72_ba135u_26_S24, rnauvb1_81_ta147d_0_S9, rnauvb1_84_ta147d_26_S12 ))

#coldata <- subset(coldata, (rownames(coldata) %in% c("rnauvb1_101_ta168r_0_S5", "rnauvb1_104_ta168r_26_S8", "rnauvb1_109_uc194f_0_S13", "rnauvb1_112_uc194f_26_S16",  "rnauvb1_65_ba162v_0_S17", "rnauvb1_68_ba162v_26_S20", "rnauvb1_69_ba135u_0_S21", "rnauvb1_72_ba135u_26_S24", "rnauvb1_81_ta147d_0_S9", "rnauvb1_84_ta147d_26_S12" )))


## (placeholder <- c("rnauvb1_101_ta168r_0_S5", "rnauvb1_104_ta168r_26_S8", "rnauvb1_109_uc194f_0_S13", "rnauvb1_112_uc194f_26_S16",  "rnauvb1_65_ba162v_0_S17", "rnauvb1_68_ba162v_26_S20", "rnauvb1_69_ba135u_0_S21", "rnauvb1_72_ba135u_26_S24", "rnauvb1_81_ta147d_0_S9", "rnauvb1_84_ta147d_26_S12"))



#dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~timepoint)
#dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~ date + prep_code)
## Want to test for the efect of timepoint (last factor) while controlling for the effect of the sampleID (first factor)
ddsMat <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~SampleID + timepoint)

## Change Factor Level
## Can change the factor levels of variables in the design before running the DESeq2 analysis

ddsMat$timepoint <- factor(ddsMat$timepoint, levels = c("26","0")) ## F
ddsMat

# Run the DESeq pipeline
dds <- DESeq(ddsMat)
res <- results(dds)








# Plot dispersions
#A simple function that plots the per-gene dispersion estimates together with the fitted mean-dispersion relationship.

png("qc-dispersions.2.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
vst <- vst(dds)

head(assay(rld))
hist(assay(rld))


head(assay(vst))
hist(assay(vst))

##Histograms
png("rld.histo.png", 1000, 1000, pointsize=20)
hist(assay(rld))
dev.off()

png("vst.histo.png", 1000, 1000, pointsize=20)
hist(assay(vst))
dev.off()


# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(timepoint)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(timepoint))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
 col=colorpanel(100, "black", "white"),
 ColSideColors=mycols[timepoint], RowSideColors=mycols[timepoint],
 margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="timepoint")


##20

pdf("pca_timepoint.20.pdf", width = 20, height = 20)
DESeq2::plotPCA(rld, intgroup="timepoint")
dev.off()
pdf("pca_SampleID.20.pdf", width = 20, height = 20)
DESeq2::plotPCA(rld, intgroup="SampleID")
dev.off()

## 15
pdf("pca_timepoint.15.pdf", width = 15, height = 15)
DESeq2::plotPCA(rld, intgroup="timepoint")
dev.off()
pdf("pca_SampleID.15.pdf", width = 15, height = 15)
DESeq2::plotPCA(rld, intgroup="SampleID")
dev.off()

##10
pdf("pca_timepoint.10.pdf", width = 10, height = 10)
DESeq2::plotPCA(rld, intgroup="timepoint")
dev.off()
pdf("pca_SampleID.10.pdf", width = 10, height = 10)
DESeq2::plotPCA(rld, intgroup="SampleID")
dev.off()



# PCA -----------Raw counts--------------------------------

project.pca <- prcomp(t(countdata))
summary(project.pca)
#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100

##Scree plot
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))


png("scree-plot.counts.png", 1000, 1000, pointsize=20)
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))
dev.off()

## Pairs Plot

par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(project.pca$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
pairs(project.pca$x[,6:10], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)


png("bi-plot.counts.1-5.png", 1000, 1000, pointsize=20)
pairs(project.pca$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
dev.off()
png("bi-plot.counts.6-10.png", 1000, 1000, pointsize=20)
pairs(project.pca$x[,6:10], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)
dev.off()

png("bi-plot.counts.11-20.png", 1000, 1000, pointsize=20)
pairs(project.pca$x[,11:20], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)
dev.off()


## Bi-Plots
png("bi-plot.counts.triple.png", 1000, 1000, pointsize=20)
par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)
#Plots scatter plot for PC 1 and 2
plot(project.pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
points(project.pca$x, col="black", pch=16, cex=1)
#Plots scatter plot for PC 1 and 3
plot(project.pca$x[,1], project.pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"))
points(project.pca$x[,1], project.pca$x[,3], col="black", pch=16, cex=1)
#Plots scatter plot for PC 2 and 3
plot(project.pca$x[,2], project.pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"), ylab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"))
points(project.pca$x[,2], project.pca$x[,3], col="black", pch=16, cex=1)
dev.off()



# PCA -----------Normalised counts rld--------------------------------

library(genefilter)
ntop <- 500
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay(rld)[select, ] )
project.pca <- prcomp(t(mat))





summary(project.pca)
#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100

##Scree plot
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))


png("scree-plot.rld.png", 1000, 1000, pointsize=20)
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))
dev.off()

## Pairs Plot

par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(project.pca$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
pairs(project.pca$x[,6:10], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)


png("bi-plot.rld.1-5.png", 1000, 1000, pointsize=20)
pairs(project.pca$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
dev.off()
png("bi-plot.rld.6-10.png", 1000, 1000, pointsize=20)
pairs(project.pca$x[,6:10], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)
dev.off()

png("bi-plot.rld.11-20.png", 1000, 1000, pointsize=20)
pairs(project.pca$x[,11:20], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)
dev.off()


## Bi-Plots
png("bi-plot.rld.triple.png", 1000, 1000, pointsize=20)
par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)
#Plots scatter plot for PC 1 and 2
plot(project.pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
points(project.pca$x, col="black", pch=16, cex=1)
#Plots scatter plot for PC 1 and 3
plot(project.pca$x[,1], project.pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"))
points(project.pca$x[,1], project.pca$x[,3], col="black", pch=16, cex=1)
#Plots scatter plot for PC 2 and 3
plot(project.pca$x[,2], project.pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"), ylab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"))
points(project.pca$x[,2], project.pca$x[,3], col="black", pch=16, cex=1)
dev.off()


# PCA -----------Normalised counts vst--------------------------------







library(genefilter)
ntop <- 500
rv <- rowVars(assay(vst))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay(rld)[select, ] )
project.pca <- prcomp(t(mat))





summary(project.pca)
#Determine the proportion of variance of each component
#Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100

##Scree plot
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))


png("scree-plot.vst.png", 1000, 1000, pointsize=20)
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))
dev.off()

## Pairs Plot

par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(project.pca$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
pairs(project.pca$x[,6:10], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)


png("bi-plot.vst.1-5.png", 1000, 1000, pointsize=20)
pairs(project.pca$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
dev.off()
png("bi-plot.vst.6-10.png", 1000, 1000, pointsize=20)
pairs(project.pca$x[,6:10], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)
dev.off()

png("bi-plot.vst.11-20.png", 1000, 1000, pointsize=20)
pairs(project.pca$x[,11:20], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)
dev.off()
png("bi-plot.vst.1-20.png", 1000, 1000, pointsize=20)
pairs(project.pca$x[,1:20], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)
dev.off()


## Bi-Plots
png("bi-plot.vst.triple.png", 1000, 1000, pointsize=20)
par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)
#Plots scatter plot for PC 1 and 2
plot(project.pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
points(project.pca$x, col="black", pch=16, cex=1)
#Plots scatter plot for PC 1 and 3
plot(project.pca$x[,1], project.pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"))
points(project.pca$x[,1], project.pca$x[,3], col="black", pch=16, cex=1)
#Plots scatter plot for PC 2 and 3
plot(project.pca$x[,2], project.pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"), ylab=paste("PC3, ", round(project.pca.proportionvariances[3], 2), "%"))
points(project.pca$x[,2], project.pca$x[,3], col="black", pch=16, cex=1)
dev.off()


###--- End of PCA
##################
#####  MA  #######
##################

png("MA-default", width = 1128, height = 805 )
DESeq2::plotMA(dds)
dev.off()

png("MA-limits2", width = 1128, height = 805)
DESeq2::plotMA(dds, ylim=c(-2,2))
dev.off()

png("MA-limits.wide", width = 1128, height = 805)
DESeq2::plotMA(dds, ylim=c(-3,4.5))
dev.off()

plotMA(dds, colNonSig = "blue", main="MA plot")
abline(h=c(-1:1), col="red")



# Colors for plots below
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Set1")[1:length(unique(timepoint))])
condition <- factor(timepoint)

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
 col=colorpanel(100, "black", "white"),
 ColSideColors=mycols[condition], RowSideColors=mycols[condition],
 margin=c(10, 10), main="Sample Distance Matrix")
dev.off()




##HeatMap

library(gplots) 
library("RColorBrewer")
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
paste(timepoint, SampleID, sep=" : "))
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))




#CREATE VOLCANO PLOT
#INCLUDES A THRESHOLD FOR P-VALUES THAT TEND TOWARDS ZERO WHICH WOULD
#OTHERWISE BE PLOTTED OUTSIDE THE AXIS LIMITS

#HIGHLIGHTS GENES SIGNIFICANTLY DYSREGULATED AND A LOG2 FOLD CHANGE GREATER THAN 4

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)

volcano <- data.frame(subset(res, select=c(log2FoldChange,pvalue,padj)))
volcano <- rownames_to_column(volcano)

#results = mutate(volcano, sig=ifelse(rv_volcano$padj<0.05 & abs(rv_volcano$log2FoldChange)>1, "Sig","Not Sig"))
results = mutate(volcano, sig=ifelse(volcano$padj<0.05 & abs(volcano$log2FoldChange)>2, "Sig","Not Sig"))
indx <- results$pvalue < abs(1e-80) 
results$pvalue[indx & !is.na(indx)] <- 1e-80



#Subset Data
subset(countdata, select=-c(rnauvb1_37_ba119q_0_S21, rnauvb1_36_ba179f_26_S20))
subset(coldata, !(rownames(coldata) %in% c("rnauvb1_37_ba119q_0_S21", "rnauvb1_36_ba179f_26_S20")))

subset(results, (results$sig %in% c("Sig", "Not Sig")))

ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("black", "red")) +geom_text_repel(data=filter(results, padj<0.05 & abs(log2FoldChange)>2), aes(label=rowname)) + ggtitle("Week 0 vs Week 26") + theme(plot.title = element_text(hjust = 0.5,size=12)) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5))


png("volcano-plot-basic.png", w=1000, h=1000, pointsize=20)
ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("black", "red"))
dev.off()


pdf("volcano-plot-basic.pdf", w=15, h=15, pointsize=20)
ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("black", "red"))
dev.off()


pdf("volcano-plot-annotation.pdf", w=15, h=15, pointsize=20)
ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("black", "red")) +geom_text_repel(data=filter(results, padj<0.05 & abs(log2FoldChange)>2), aes(label=rowname)) + ggtitle("Week 0 vs Week 26") + theme(plot.title = element_text(hjust = 0.5,size=12)) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5))
dev.off()

data=filter(results, padj<0.05 & abs(log2FoldChange)>2)

write.csv(data, file="diffexp-filtered.log2.csv")


#results = mutate(volcano, sig=ifelse(rv_volcano$padj<0.05 & abs(rv_volcano$log2FoldChange)>1, "Sig","Not Sig"))
results = mutate(volcano, sig=ifelse(volcano$padj<0.05 & abs(volcano$log2FoldChange)>1, "Sig","Not Sig"))
indx <- results$pvalue < abs(1e-80) 
results$pvalue[indx & !is.na(indx)] <- 1e-80

ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("black", "red")) +geom_text_repel(data=filter(results, padj<0.05 & abs(log2FoldChange)>2), aes(label=rowname)) + ggtitle("Week 0 vs Week 26") + theme(plot.title = element_text(hjust = 0.5,size=12)) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5))


pdf("volcano-1.5-annotation.pdf", w=15, h=15, pointsize=20)
ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("black", "red")) +geom_text_repel(data=filter(results, padj<0.05 & abs(log2FoldChange)>1.5), aes(label=rowname)) + ggtitle("Week 0 vs Week 26") + theme(plot.title = element_text(hjust = 0.5,size=12)) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5))
dev.off()






## Volcano Plot

## qvalue < 0.01 & logfold <(abs)1

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)

volcano <- data.frame(subset(res, select=c(log2FoldChange,pvalue,padj)))
volcano <- rownames_to_column(volcano)

results = mutate(volcano, sig=ifelse(volcano$padj<0.01 & abs(volcano$log2FoldChange)>1, "Sig","Not Sig"))
indx <- results$pvalue < abs(1e-80) 
results$pvalue[indx & !is.na(indx)] <- 1e-80
subset(results, (results$sig %in% c("Sig", "Not Sig")))


ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("black", "red")) +geom_text_repel(data=filter(results, padj<0.01 & abs(log2FoldChange)>1), aes(label=rowname)) + ggtitle("Week 0 vs Week 26") + theme(plot.title = element_text(hjust = 0.5,size=12)) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5))




## qvalue < 0.05 & logfold <(abs)1


library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)

volcano <- data.frame(subset(res, select=c(log2FoldChange,pvalue,padj)))
volcano <- rownames_to_column(volcano)

results = mutate(volcano, sig=ifelse(volcano$padj<0.05 & abs(volcano$log2FoldChange)>1, "Sig","Not Sig"))
indx <- results$pvalue < abs(1e-80) 
results$pvalue[indx & !is.na(indx)] <- 1e-80
subset(results, (results$sig %in% c("Sig", "Not Sig")))


ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("black", "red")) +geom_text_repel(data=filter(results, padj<0.05 & abs(log2FoldChange)>1), aes(label=rowname)) + ggtitle("Week 0 vs Week 26") + theme(plot.title = element_text(hjust = 0.5,size=12)) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5))




p1 <- ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(labels = c("q-value > 0.05 && logfold <1", "q-value < 0.05 && logfold >1"), values=c("black", "red")) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5)) + ggtitle("Week 0 vs Week 26") 



## qvalue < 0.01 & logfold <(abs)1

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)

volcano <- data.frame(subset(res, select=c(log2FoldChange,pvalue,padj)))
volcano <- rownames_to_column(volcano)

results = mutate(volcano, sig=ifelse(volcano$padj<0.01 & abs(volcano$log2FoldChange)>1, "Sig","Not Sig"))
indx <- results$pvalue < abs(1e-80) 
results$pvalue[indx & !is.na(indx)] <- 1e-80
results <- subset(results, (results$sig %in% c("Sig", "Not Sig")))


ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(labels = c("q-value > 0.01 && logfold <1", "q-value < 0.01 && logfold >1"), values=c("black", "red")) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5)) + ggtitle("Week 0 vs Week 26") 


p2 <- ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(labels = c("q-value > 0.01 && logfold <1", "q-value < 0.01 && logfold >1"), values=c("black", "red")) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5)) + ggtitle("Week 0 vs Week 26") 



ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(labels = c("q-value > 0.01 && logfold <1", "q-value < 0.01 && logfold >1"), values=c("black", "red")) +geom_text_repel(data=filter(results, padj<0.01 & abs(log2FoldChange)>1.5), aes(label=rowname)) + ggtitle("Week 0 vs Week 26") + theme(plot.title = element_text(hjust = 0.5,size=12)) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5))



data=filter(results, padj<0.05 & abs(log2FoldChange)>1)
write.csv(data, file="diffexp.qval<0.05_&&_log1.csv")


data=filter(results, padj<0.01 & abs(log2FoldChange)>1)
write.csv(data, file="diffexp.qval<0.01_&&_log1.csv")




data=filter(results, padj<0.01 & abs(log2FoldChange)>1.5)
write.csv(data, file="diffexp.qval<0.01_&&_log1.5.csv")



data=filter(results, padj<0.05 & abs(log2FoldChange)>1.5)

write.csv(data, file="diffexp-filtered.log1.5.csv")





###HERE


# Plot 1 - colour, no annotation
ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(labels = c("q-value > 0.01 && logfold <1", "q-value < 0.01 && logfold >1"), values=c("black", "red")) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5)) + ggtitle("Week 0 vs Week 26") 

# Plot 2- colour, annotations
ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(labels = c("q-value > 0.01 && logfold <1", "q-value < 0.01 && logfold >1"), values=c("black", "red")) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5)) + ggtitle("Week 0 vs Week 26") +geom_text_repel(data=filter(results, padj<0.01 & abs(log2FoldChange)>1.5), aes(label=rowname))


# Plot 3 - shapes, no annotation


ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(shape=sig)) + scale_shape_manual(labels = c("q-value > 0.01 && logfold <1", "q-value < 0.01 && logfold >1"), values=c(22, 16)) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5)) + ggtitle("Week 0 vs Week 26") 


# Plot 4 - shapes, annotation


ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(shape=sig)) + scale_shape_manual(labels = c("q-value > 0.01 && logfold <1", "q-value < 0.01 && logfold >1"), values=c(22, 16)) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5)) + ggtitle("Week 0 vs Week 26")  +geom_text_repel(data=filter(results, padj<0.01 & abs(log2FoldChange)>1.5), aes(label=rowname))



## qvalue < 0.01 & logfold <(abs)1

library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)

volcano <- data.frame(subset(res, select=c(log2FoldChange,pvalue,padj)))
volcano <- rownames_to_column(volcano)

results = mutate(volcano, sig=ifelse(volcano$padj<0.01 & abs(volcano$log2FoldChange)>1.2, "Sig","Not Sig"))
indx <- results$pvalue < abs(1e-80) 
results$pvalue[indx & !is.na(indx)] <- 1e-80
subset(results, (results$sig %in% c("Sig", "Not Sig")))


ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) +  scale_color_manual(labels = c("q-value > 0.01 && logfold <1.2", "q-value < 0.01 && logfold >1.2"), values=c("black", "red")) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5)) + ggtitle("Week 0 vs Week 26") 





ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(labels = c("q-value > 0.01 && logfold <1", "q-value < 0.01 && logfold >1"), values=c("black", "red")) +geom_text_repel(data=filter(results, padj<0.01 & abs(log2FoldChange)>1.5), aes(label=rowname)) + ggtitle("Week 0 vs Week 26") + theme(plot.title = element_text(hjust = 0.5,size=12)) + scale_x_continuous(breaks = seq(-3, 5, 1), limits=c(-3, 5))












# Get differential expression results
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
 with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
 with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
 if (labelsig) {
 require(calibrate)
 with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
 }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
ggplot(results, aes(log2FoldChange, -log10(pvalue))) + geom_point(aes(col=sig)) + scale_color_manual(values=c("black", "red")) +geom_text_repel(data=filter(results, padj<0.05 & abs(log2FoldChange)>2), aes(label=rowname)) + xlim(-4,5) + ggtitle("Week 0 vs Week 26") + theme(plot.title = element_text(hjust = 0.5,size=12))
