
setwd("/path/to/PBMC900_TruthKnown/")
#######Part1: Default Seurat pipeline############

# read in the counts data in the MTX format
pbmc <- Seurat::Read10X(data.dir = ".")

# create SeuratObject from the counts data
pbmc <- Seurat::CreateSeuratObject(counts = pbmc)

# perform default Seurat logCP10k normalization
pbmc <- Seurat::NormalizeData(object = pbmc)

# shortlist HVGs - 2000 HVGs is the default in Seurat
pbmc <- Seurat::FindVariableFeatures(object = pbmc)

# import the known-labels for the cells
cell_label<-read.table("cell_labels.txt",header = T)
cell_label<-factor(cell_label$cell_type_labels,levels = c("regulatory_t","naive_t","cd34","memory_t","cd56_nk","cytotoxic_t","cd14_monocytes","b_cells","naive_cytotoxic"))

pbmc$cell_label <- cell_label

# standardize the counts matrix for the HVGs - mean = 0 and std. dev. = 1
pbmc<-Seurat::ScaleData(pbmc,verbose=FALSE)

# specify the number of PCs you wish to shortlist for kNN and UMAP
npcs <- 30

# perform PCA
suppressWarnings(pbmc<-Seurat::RunPCA(pbmc,npcs=npcs,verbose=FALSE))

# using PC coordinates of cells as input apply UMAP
pbmc<-Seurat::RunUMAP(pbmc,dims = 1:npcs,verbose=FALSE)

# plot the UMAP
Seurat::DimPlot(pbmc,reduction = "umap",group.by = "cell_label",shuffle = T)

pdf(paste0("umap_default-hvg_Seurat.pdf"),height = 5.5,width = 7)
Seurat::DimPlot(pbmc,reduction = "umap",group.by = "cell_label",shuffle = T)
dev.off() 

#######Part2: Seurat HVGs replaced by Piccolo HVGs#########

# read in the counts data in the MTX format
pbmc <- Seurat::Read10X(data.dir = ".")

# create SeuratObject from the counts data
pbmc <- Seurat::CreateSeuratObject(counts = pbmc)

# perform default Seurat logCP10k normalization
pbmc <- Seurat::NormalizeData(object = pbmc)

# shortlist HVGs - 2000 HVGs is the default in Seurat
pbmc <- Seurat::FindVariableFeatures(object = pbmc)

pbmc$cell_label <- cell_label

# specify Piccolo HVGs as the variable features for scaling
PiccoloHVG <- read.table("Piccolo1500HVG.txt",header = T)$GeneID

Seurat::VariableFeatures(pbmc) <- PiccoloHVG

# standardize the counts matrix for the HVGs - mean = 0 and std. dev. = 1
pbmc <- Seurat::ScaleData(pbmc,verbose=FALSE)

# specify the number of PCs you wish to shortlist for kNN and UMAP
npcs <- 30

# perform PCA
suppressWarnings(pbmc<-Seurat::RunPCA(pbmc,npcs=npcs,verbose=FALSE))

# using PC coordinates of cells as input apply UMAP
pbmc<-Seurat::RunUMAP(pbmc,dims = 1:npcs,verbose=FALSE)

# plot the UMAP
Seurat::DimPlot(pbmc,reduction = "umap",group.by = "cell_label",shuffle = T)

pdf(paste0("umap_Piccolo-hvg_Seurat.pdf"),height = 5.5,width = 7)
Seurat::DimPlot(pbmc,reduction = "umap",group.by = "cell_label",shuffle = T)
dev.off() 