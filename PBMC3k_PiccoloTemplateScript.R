#Piccolo script for 10X PBMC 3k data set

setwd("~/Documents/PBMC3k_Data/")

#Create PiccoloList object
pbmc3k <- Piccolo::CreatePiccoloList(MTX = "10X_PBMC3k_matrix.mtx.gz",Genes = "10X_PBMC3k_features.tsv",Barcodes = "10X_PBMC3k_barcodes.tsv")

#Filter cells based on some basic thresholds
pbmc3k <- Piccolo::FilterCells(PiccoloList = pbmc3k,MinFeaturesPerCell = 100,MT.Perc = 50,RP.Perc = 80,TotalCountsMADLow = NULL,TotalCountsMADHigh = NULL) #Will filter out cells with fewer than 100 genes with non-zero counts, cells with percentage of total counts contributed by mitochondrial genes > 50%, cells with percentage of total counts contributed by ribosomal genes > 80%, and cells whose total counts are more than 3.5 median absolute deviation lesser than the median total count  

#Shortlist HVGs and stable genes

pbmc3k <- Piccolo::SelectFeatures(PiccoloList = pbmc3k) #Default

pbmc3k <- Piccolo::SelectFeatures(PiccoloList = pbmc3k,NoOfHVG = 3000) #Number of HVGs to be shortlisted set to 3000

#Perform normalization

pbmc3k <- Piccolo::Normalize(PiccoloList = pbmc3k) #Default

pbmc3k <- Piccolo::Normalize(PiccoloList = pbmc3k,Transform = "bc") #Box-Cox power law transform

#Apply principal components analysis (PCA)

pbmc3k <- Piccolo::ComputePC(PiccoloList = pbmc3k) #Default - top 50 PCs will be shortlisted

#Get UMAP coordinates

pbmc3k <- Piccolo::UMAPcoords(PiccoloList = pbmc3k)

#Identify k Nearest Neighhours

pbmc3k <- Piccolo::KNearestNeighbors(PiccoloList = pbmc3k) #Default k=10

pbmc3k <- Piccolo::KNearestNeighbors(PiccoloList = pbmc3k,k = 15) #number of nearest neighbours (k) set to 15

#Perform Leiden clustering

pbmc3k <- Piccolo::LeidenClustering(PiccoloList = pbmc3k) #Default Resolution = 1

pbmc3k <- Piccolo::LeidenClustering(PiccoloList = pbmc3k,Resolution = 2)

#Plot UMAP showing the clusters

Piccolo::LeidenUMAP(PiccoloList = pbmc3k,Levels = 1:length(unique(pbmc3k$ClusterLabels)),LegendPosition = "none")

#Perform differential expression analysis between cells in each cluster and the rest

pbmc3k <- Piccolo::PerformDiffExpClusterwise(PiccoloList = pbmc3k)

pbmc3k <- Piccolo::PerformDiffExpClusterwise(PiccoloList = pbmc3k,Method = "wilcoxon", Out = T) #Will generate output .csv files in the working directory

#Obtain marker gene sets
pbmc3k <- Piccolo::GetClusterMarkers(PiccoloList = pbmc3k,MeanDiffSD = 2.5,MaxSharedClusters = 2,Out = T)

#Visualize expression profiles of marker gene sets in UMAP

setwd("~/Documents/PBMC3k_Data/ClusterMarkers_MeanDiffSD2.5_FDR0.005_SharedClusters2/")

ClusterOfInterest <- 6

MarkerGeneSet <- read.csv(paste0("Cluster",ClusterOfInterest,".csv"))$x

#Compute composite z-scores for gene set (using Stouffer's method to combine the z-scores)
pbmc3k <- Piccolo::GeneSetZScores(PiccoloList = pbmc3k,Genes = MarkerGeneSet)

Piccolo::UMAPZScores(PiccoloList = pbmc3k,Name = paste0("Cluster ",ClusterOfInterest," Markers"))
