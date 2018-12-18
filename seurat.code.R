#! /usr/bin/env RScript
###   ------------------------ Libraries ------------------------   ### 
library(Seurat)
library(dplyr)
library(DESeq2)
library(ggplot2)

###   ------------------------ End of Libraries ------------------------   ### 



###   ------------------------ Functions for tpm -----------------------   ###  
tpm <- function(counts, lengths) {  
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

## read table from featureCounts output
ftr.cnt <- read.table("2H.count", sep="\t", stringsAsFactors=FALSE, header=TRUE)
library(dplyr)  
library(tidyr)

ftr.tpm <- ftr.cnt %>%  
  gather(sample, cnt, 8:ncol(ftr.cnt)) %>%
  mutate(tpm=tpm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, tpm)
write.table(ftr.tpm, "counts_tpm.txt", sep="\t", row.names=FALSE, quote=FALSE)
###   ------------------------ End of tpm -----------------------   ### 



###   ------------------------ Data Load -----------------------   ### 


# Load in filtered count matrix
rna.2H.tpm <- read.csv("counts_tpm.txt",sep = "\t", header = TRUE, row.names = 1)
rna.2H.tpm <- rna.2H.tpm[ ,7:ncol(rna.2H.tpm)] ## removing meta data
write.table(rna.2H.tpm, "counts_tpm.red.txt", sep="\t", row.names=FALSE, quote=FALSE)


rna.metadata <- read.csv("metadata.tab",sep = "\t", header = TRUE)

###
csf <- CreateSeuratObject(raw.data = csf.data, min.cells = 3, min.genes = 200, project = 'CSF_2018')

csf <- AddMetaData(object = csf, metadata = sc_cell_info, col.name = c('Twin', 'Sort', 'Clones', 'Case', 'CellType'))


# Initialize the Seurat object with the raw (non-normalized data).  Keep all genes expressed in >= 3 cells (~0.1% of the data).
# Keep all cells with at least 200 detected genes
tb.rna <- CreateSeuratObject(raw.data = rna.2H.tpm, min.cells = 3, min.genes = 200, project = "2H_TB")

# Add metadata
tb.rna <- AddMetaData(object = tb.rna, metadata = rna.metadata, col.name = c('Well', 'ID', 'antigen', 'TAM_TB', 'Clinical', 'Cells', 'Gate', 'Sort', 'Study')) 




#####  ------------------------ Seurat -----------------------   ### 

#  Setup a Seurat object, and cluster cells based on RNA expression

# Initialize the Seurat object with the tpm-normalized data. 
tb.rna <- CreateSeuratObject(raw.data = rna.2H.tpm, meta.data = rna.metadata)

# Initialize the Seurat object with the tpm-normalized data).  Keep all genes expressed in >= 3 cells (~0.1% of the data). 
# Keep all cells with at least 200 detected genes
tb.rna.flt <- CreateSeuratObject(raw.data = rna.2H.tpm, min.cells = 3, min.genes = 200, project = "2H_TB")


### QC

mito.genes <- grep(pattern = "^MT-", x = rownames(x = tb.rna@data), value = TRUE)
percent.mito <- Matrix::colSums(tb.rna@raw.data[mito.genes, ])/Matrix::colSums(tb.rna@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
tb.rna <- AddMetaData(object = tb.rna, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = tb.rna, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = tb.rna, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = tb.rna, gene1 = "nUMI", gene2 = "nGene")



# Filter out cells

# We filter out cells that have unique gene counts over INF or less than 200 
# Note that low.thresholds and high.thresholds are used to define a 'gate'.  
# -Inf and Inf should be used if you don't want a lower or upper threshold.
tb.rna <- FilterCells(object = tb.rna, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 0.6))



## Normalizing data
## Normalises the gene exp measurements for each cell by the total expression and multiples this by a scale factor of 10000
tb.rna <- NormalizeData(object = tb.rna, normalization.method = "LogNormalize", scale.factor = 10000)

## Detection of variable genes across the single cells
## The parameters here identify ~2,000 variable genes, and represent typical parameter settings for UMI data that is normalized to a total of 1e4 molecules.
tb.rna <- FindVariableGenes(object = tb.rna, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = FALSE)


# Identifies genes that are outliers on a "mean variability plot", enabling identification of variable genes and controlling for the relationship between variability and average expression. 

# tb.rna <- FindVariableGenes(object = tb.rna,	# Seurat Object
#	mean.function = ExpMean, 		# Computes the x-axis value (average exp), the mean of the non-zero values
#	dispersion.function = LogVMR, 		# Computes the y-axis value (dispersion), 
##	x.low.cutoff = 0.0125, 			# Bottom cutoff on the x-axis for identifying variable genes
#	x.high.cutoff = 3, 			# Top cutoff on the x-axis for identifying variable genes
#	y.cutoff = 0.5)				# Bottom cutoff on the x-axis for identifying variable genes 



length(x = tb.rna@var.genes)

# scale data nad remove unwanted sources of variation
# uninteresting sources of variation - technical noise, batch effects, bio variation (eg cell cycle stages) 
# n this simple example here for post-mitotic blood cells, we regress on the number of detected molecules per cell as well as the percentage mitochondrial gene content.

tb.rna <- ScaleData(object = tb.rna, vars.to.regress = c("nUMI", "percent.mito"))


## Perform linear dimensional reduction
# PCA is performed on the scaled data, using the var.genes as default.  Use of the highly variable genes can improvde performance

tb.rna <- RunPCA(object = tb.rna, pc.genes = tb.rna@var.genes, do.print = TRUE, pcs.print = 1:5,  genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = tb.rna, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

VizPCA(object = tb.rna, pcs.use = 1:2)
# top 9 factors
VizPCA(object = tb.rna, pcs.use = 1:9)
PCAPlot(object = tb.rna, dim.1 = 1, dim.2 = 2)


# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
tb.rna <- ProjectPCA(object = tb.rna, do.print = FALSE)

#PC heatmap 
# Allows easy exploration of the sources of heterogeneity in a dataset useful for deciding which PCs to use in futher analyses
# Both cells and genes are ordered according to their PCA scores. 
# A supervised analysis, however useful for exploring correlated gene sets.

PCHeatmap(object = tb.rna, pc.use = 1, cells.use = 96, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = tb.rna, pc.use = 1:12, cells.use = 96, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

# Determine statistically significant principal components

# To overcome the extensive technical noise in any single gene for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a ‘metagene’ that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step.

# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
tb.rna <- JackStraw(object = tb.rna, num.replicate = 100, display.progress = FALSE)
#The JackStrawPlot function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line). In this case it appears that PCs 1-10 are significant.
JackStrawPlot(object = tb.rna, PCs = 1:12)
PCElbowPlot(object = tb.rna)

## These three graphs should help determine which clusters to take on.
## ? 



## Cluster the cells
# contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.
# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)


tb.rna <- FindClusters(object = tb.rna, reduction.type = "pca", dims.use = 1:10, 
                       resolution = 0.6, print.output = 0, save.SNN = TRUE)

# shows the chosen parameters
PrintFindClustersParams(object = tb.rna)

#Run Non-linear dimensional reduction (tSNE)


tb.rna <- RunTSNE(object = tb.rna, dims.use = 1:10, do.fast = TRUE, perplexity = 10)

#Run t-SNE dimensionality reduction on selected features. Has the option of running in a reduced dimensional space (i.e. spectral tSNE, recommended), or running based on a set of genes. For details about stored TSNE calculation parameters, see PrintTSNEParams.

#tb.rna <- RunTSNE(object = tb.rna, 	# Seurat object
#	dims.use = 1:10, 		# which dimesions to use as input
#	do.fast = TRUE,			# 
#	reduction.use = PCA,		# default is PCA
#	cells.use = all,		# Shows which cells to analyse, default all cells
#	tsne.method = Rtsne		# Rtsne(default), tsne, FIT-SNE

# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = tb.rna)

# Save file :

saveRDS(tb.rna, file = "tb.rna.seurat.rds")

pdf("tsne.pdf", width = 10, height = 10)
TSNEPlot(object = tb.rna)
dev.off()


## Finding differentially expressed genes (cluster biomarkers)

#Seurat can help you find markers that define clusters via differential expression. By default, it identifes positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. FindAllMarkers automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

#The min.pct argument requires a gene to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a gene to be differentially expressed (on average) by some amount between the two groups. You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of genes that are unlikely to be highly discriminatory. As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in power, the speed increases can be significiant and the most highly differentially expressed genes will likely still rise to the top.



# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = tb.rna, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = tb.rna, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))

cluster1.markers <- FindMarkers(object = tb.rna, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)


VlnPlot(object = tb.rna, features.plot = c("MS4A1", "CD79A"))

VlnPlot(object = tb.rna, features.plot = c("NKG7", "PF4"), use.raw = TRUE, y.log = TRUE)


FeaturePlot(object = tb.rna, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", 
                                               "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")


FeaturePlot(object = tb.rna, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14","FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), reduction.use = "tsne")


FeaturePlot(object = tb.rna, features.plot = c("PPD", "GNLY", "CD3E", "CD14","FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), reduction.use = "tsne")


FeaturePlot(object = tb.rna, features.plot = c("PPD"),cols.use = c("grey", "blue"), reduction.use = "tsne")



#####  ------------------------ Seurat end ---------------------------------------------------
