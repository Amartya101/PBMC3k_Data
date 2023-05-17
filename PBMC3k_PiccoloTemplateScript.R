#Install Piccolo if it hasn't been installed yet
install.packages("devtools")

devtools::install_github("Amartya101/Piccolo")

#Piccolo script for 10X PBMC 3k data set

setwd("~/Documents/PBMC3k_Data/") #Specify the path to the directory that contains the 10X PBMC3k data

#Create PiccoloList object
PiccoloList <- Piccolo::CreatePiccoloList(MTX = "10X_PBMC3k_matrix.mtx.gz",Genes = "10X_PBMC3k_features.tsv",Barcodes = "10X_PBMC3k_barcodes.tsv")
### Note: Make sure all 3 files are in the setwd() folder
###
### Creates PiccoloList object with 3 list items: "CountsOriginal", "GenesOriginal", "BarcodesOriginal","Counts","Genes","Barcodes"
### Counts ---> contains the read counts for each cell and each gene
### Genes ---> contains information on the genes, including identifiers
### Barcodes ---> Barcodes for each cell


#Filter cells based on some basic thresholds
PiccoloList <- Piccolo::FilterCells(PiccoloList = PiccoloList,MinFeaturesPerCell = 100,MT.Perc = 50,RP.Perc = 80,TotalCountsMADLow = NULL,TotalCountsMADHigh = NULL) #Will filter out cells with fewer than 100 genes with non-zero counts, cells with percentage of total counts contributed by mitochondrial genes > 50%, cells with percentage of total counts contributed by ribosomal genes > 80%, and cells whose total counts are more than 3.5 median absolute deviation lesser than the median total count  
### Parameters:
### MinFeaturesPerCell = 100 ---> If < 100 genes have non-zero counts for each cell, filtered out
### MT.Perc = 50 ---> if > 50% of the total counts for a given cell are mitochonodrial genes, filtered out
### RB.Perc = 50 ---> if > 50% of the total counts for a given cell are ribosomal genes, filtered out
### TotalCountsMADLow = 3.5 ---> filter out cells that have a total count < -3.5 SD from the overall cell count distribution (compared to other cells)
### TotalCountsMADLow = 3.5 ---> filter out cells that have a total count > 3.5 SD from the overall cell count distribution (compared to other cells); This may eliminate doublets
###
### Adding 3 list entries: "Counts", "Genes", "Barcodes"
### Counts ---> contains the read counts for each cell that passed the filters and each gene
###### Note: dim(PiccoloList$Counts) returns number of genes / cells after filtration
###### Note: dim(PiccoloList$CountsOriginal) returns number of genes / cells before filtration
###### Note: You can remove the original counts to save up memory by typing: PiccoloList$CountsOriginal <- NULL
### Genes ---> contains information on the genes, including identifiers
### Barcodes ---> Barcodes for each cell that passed the filters


#Shortlist HVGs and stable genes
PiccoloList <- Piccolo::SelectFeatures(PiccoloList = PiccoloList) #Default
### Adds 8 list entries: "RelevantGenes", "RelevantGenes.Ser.Nos", "DiffAlpha", "DispCoef", "HVG", "HVG.Ser.Nos", "StableGenes", "Stable.Ser.Nos"
### RelevantGenes ---> genes that have variance above variance expected under Poisson Distribution (1st basic gene filtration criteria)
### HVG ---> list of highly variable genes (highly variable defined by default as top 90% genes in each bin)
### StableGenes ---> list of genes that are in Relevant Genes but in HVG (used for future normalization); Note: if not enough genes, the genes that did not meet the Poisson criteria above are later (during normalization) included into that list 
###
### Also provided, but not necessary to look at:
### RelevantGenes.Ser.Nos HVG.Ser.Nos Stable.Ser.Nos ---> Serial numbers of the genes for each gene list
### DiffAlpha ---> value used to determine whether a gene is variable (dispersion coefficient - reference dispersion coefficient per bin; each bin contains about 10 genes of similar expression levels)
### DispCoef ---> raw dispersion coefficient for each gene

############### Other option ################
#############################################
# PiccoloList <- Piccolo::SelectFeatures(PiccoloList = PiccoloList,NoOfHVG = 3000) #Number of HVGs to be shortlisted set to 3000
### Selects the top 3000 variable genes based on DiffAlpha
### NoOfHVG = any ---> includes all RelevantGenes into HVG
#############################################
#############################################


#Perform normalization
PiccoloList <- Piccolo::Normalize(PiccoloList = PiccoloList) #Default
### By default, log-normalizes
###
### Adding 3 list entries: "NormCounts", "SizeFactors", "RevisedAlpha"
### NormCounts ---> contains the residuals (z-scores) for the HVG
### SizeFactors ---> total counts across StableGenes for a given cell / mean (total counts across StableGenes for all cells)
###
### Also provided, but not necessary to look at:
### RevisedAlpha ---> dispersion coefficient recalculated after the counts have been adjusted by the SizeFactors
###### Note: some genes that previously passed the Poisson criteria may no longer pass it with the adjusted counts; If that happens, both HVG and NormCounts no longer carry those genes moving forward

############### Other option ################
#############################################
# PiccoloList <- Piccolo::Normalize(PiccoloList = PiccoloList,Transform = "bc") #Box-Cox power law transform
### Transform = bc ---> box-cox transform; works better than log, BUT takes computational time
### Transform = sqrt ---> square-root transform
#############################################
#############################################


#Apply principal components analysis (PCA)
PiccoloList <- Piccolo::ComputePC(PiccoloList = PiccoloList) #Default - top 50 PCs will be shortlisted
### Adds 1 list entry: "PCA"
### PCA is a list that notably contains: "d", "sdev", "x" (also contains the following but not useful: "u", "v", "niter", "nops")
###### PCA$x ---> the coordinates of each cell along the PCA components
###### PCA$d ---> contains the weights of each gene for each PCA component
###### PCA$sdev ---> contains the standard deviation per PCA component 


#Plots ScreePlot of the variance explained per PCA component
Piccolo::ScreePlot(PiccoloList = PiccoloList)
### A ScreePlot is a line plot of the eigenvalues (variance explained) of factors or principal components in an analysis
### uses PCA$sdev to plot the histogram


#Get UMAP coordinates
PiccoloList <- Piccolo::UMAPcoords(PiccoloList = PiccoloList)
### Adds 1 list entry: "UMAP"
### UMAP ---> dataframe contains: CellID (barcodes) and UMAP1 UMAP2 (UMAP coordinates)


#Identify k Nearest Neighhours
PiccoloList <- Piccolo::KNearestNeighbors(PiccoloList = PiccoloList,k = 5) #Default k=10
### Default parameter: kNN = 10
###
### Adds 1 list entry: "kNN"
### kNN ---> matrix that contains the Serial numbers of the k-nearest-neighbor cells for each cell
###### Note: only looks at the Euclidian distance based on the PCs; NO clustering yet

############### Other option ################
#############################################
# PiccoloList <- Piccolo::KNearestNeighbors(PiccoloList = PiccoloList,k = 15) #number of nearest neighbours (k) set to 15
### Can set kNN value
#############################################
#############################################


#Perform Leiden clustering
PiccoloList <- Piccolo::LeidenClustering(PiccoloList = PiccoloList,Resolution = 1) #Default Resolution = 1
### Adds 2 list entries: "Leiden", "ClusterLabels"
### ClusterLabels ---> vector that for each cell associates its cluster label
### Leiden ---> list that contains: "membership"  "nb_clusters" "quality"     "algorithm"   "vcount"  (only first one is useful)
##### Leiden is the method (Traag 2019 Sci. Rep.) used to partition the kNN graph; Note: does not take Euclidian distance into account, only the actual connections between cells
##### Leiden$membership ---> contains the clusters associated with each cell; this is extracted and put into ClusterLabels

############### Other option ################
#############################################
# PiccoloList <- Piccolo::LeidenClustering(PiccoloList = PiccoloList,Resolution = 1.5)
### Resolution = 1.5  ---> adjusts the number of clusters; larger values leads to more clusters
#############################################
#############################################


#Plot UMAP showing the clusters
Piccolo::LeidenUMAP(PiccoloList = PiccoloList,Levels = 1:length(unique(PiccoloList$ClusterLabels)),LegendPosition = "none",Title = "PBMC3k")


#Perform differential expression analysis between cells in each cluster and the rest
PiccoloList <- Piccolo::PerformDiffExpClusterwise(PiccoloList = PiccoloList,Method = "t.test", Out = T)
### Parameters:
### Method ---> specifies which statisticial test is used; 2 options: "t.test" and "wilcoxon"
### Out = T ---> generates the output tables for each cluster; Note: creates them in the working directory, so setwd() elsewhere beforehand if needed
###
### Adds 1 list entry: 'DE.Results'
### DE.Results ---> list that contains tables for each cluster
##### DE.Results$`Cluster 1` ---> Table for Cluster 1 that contains the following columns, most important among them: mean.diff, p.adj(BH), Log2FC.Means, PercNonZero.Group1, PercNonZero.Group2
##### V1	V2 ---> gene identifiers; sometimes V2 does not exist because the input file doesn't contain 2 columns of gene identifiers
##### obs.x	---> number of cells in Cluster 1
##### obs.y ---> number of cells in rest of the clusters
##### obs.tot	---> total number of cells
##### mean.x ---> mean of residuals of cells within Cluster 1 for each HVG (highly variable gene)
##### mean.y ---> mean of residuals of cells within all other clusters for each HVG
##### mean.diff ---> mean.x - mean.y 
##### var.x ---> variance of residuals of cells within Cluster 1 for each HVG
##### var.y ---> variance of residuals of cells within all other clusters for each HVG
#####	var.pooled ---> variance of residuals across all cells for each gene
##### stderr ---> standard error of the test's statistic (ex: if using the t-test, associates the t-value with its standard error)
##### df ---> degrees of freedom for the statistical test used
##### statistic ---> the value of the test's statistic (ex: if using the t-test, corresponds to t)
##### pvalue ---> p-value corresponding to the test's statistic (ex: if using the t-test, returns the p-value associated with the t-value)
##### conf.low	conf.high ---> lower and upper confidence interval of the test's statistic
##### p.adj(BH) ---> corrected p-values based on the Benjamini-Horchberg method (for multiple testing)
##### Mean.Group1	---> mean of the adjusted counts (by SizeFactors) for cells in Cluster 1
##### Mean.Group2	---> mean of the adjusted counts in all other Clusters
##### Log2FC.Means	---> log2( Mean.Group1 / Mean.Group2 )  
##### PercNonZero.Group1 ---> percentage of cells with non-zero counts in Cluster 1
##### PercNonZero.Group2 ---> percentage of cells with non-zero counts in all other clusters

############### Other option ################
#############################################
# PiccoloList <- Piccolo::PerformDiffExpClusterwise(PiccoloList = PiccoloList)
### Edits PiccoloList, but does not output the table
#############################################
#############################################


#Get marker gene sets for each cluster
PiccoloList <- Piccolo::GetClusterMarkers(PiccoloList = PiccoloList,MeanDiffSD = 2.5,FDRcutoff = 0.005,MaxSharedClusters = 2,Out = T)
### Paremeters:
### MeanDiffSD = 2.5 ---> selecting only genes with a MeanDiff > 2.5 SD away from the mean of the MeanDiff distribution across all genes (for each cluster)
### FDRcutoff = 0.005 ---> selecting genes with p.adj(BH) < 0.005
### MaxSharedClusters = 2 ---> filtering out genes that are shared in > 2 clusters 
### Out = T ---> generates .csv files for each clustering containing the list of genes; Note: creates them in the working directory, so setwd() elsewhere beforehand if needed
###
### Adds 1 list entry: 'ClusterMarkers'
### ClusterMarkers ---> list that contains tables for each cluster, tables corresponding to the DE.Results tables filtered for the genes that passed all the cut-offs 

#Save RDS
saveRDS(object = PiccoloList,file = "PBMC3k_DefaultParameters.rds")

#Generates for each cell Z-scores associated with a specific marker gene set 
PiccoloList <- Piccolo::GeneSetZScores(PiccoloList = PiccoloList,Genes = read.csv("ClusterMarkers_MeanDiffSD2.5_FDR0.005_SharedClusters2/Cluster5.csv")[[1]])
### Parameters:
### Genes ---> vector of marker gene set
##### Note: can be extracted from the output of Piccolo::GetClusterMarkers
###
### Adds 1 list entry: 'GeneSetZScore' 
### GeneSetZScore ---> Vector that contains for each cell the computed Stouffer's method z-score based on the marker gene set selected in "Genes"


#Plot z-scores UMAP for marker genes
Piccolo::UMAPZScores(PiccoloList = PiccoloList,Name = "Cluster 5 Marker Genes")
### Parameter: 
### Name ---> Title of the UMAP plot
###
### Outputs a UMAP plot colored for the Z-scores from PiccoloList$GeneSetZScore


