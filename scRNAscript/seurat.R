#create Seurat object
library(dplyr)
library(Seurat)
library(patchwork)
P5_T2.s<-CreateSeuratObject(P5_T2, project = "P5_T2", assay = "RNA")

#Seurat QC
P5_T2.s[["percent.mt"]]<-PercentageFeatureSet(P5_T2.s,pattern = "^MT-")
VlnPlot(P5_T2.s, features = c("nFeature_RNA","nCount_RNA","percent.mt"))
plot1<-FeatureScatter(P5_T2.s, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2<-FeatureScatter(P5_T2.s, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2
P5_T2.s<-subset(P5_T2.s, subset = nCount_RNA<40000 & nFeature_RNA>200 & percent.mt<40)

#Seurat NormalizeData and dim reduction
P5_T2.s<-NormalizeData(P5_T2.s)
P5_T2.s<-FindVariableFeatures(P5_T2.s,selection.method = "vst", nfeatures = 2000)
top10<-head(VariableFeatures(P5_T2.s),10)
plot1<-VariableFeaturePlot(P5_T2.s)
plot2<-LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2
all.genes.P5_T2<-rownames(P5_T2.s)
P5_T2.s<-ScaleData(P5_T2.s, features = all.genes.P5_T2)
P5_T2.s<-RunPCA(P5_T2.s, features = VariableFeatures(object = P5_T2.s))
print(P5_T2.s[["pca"]],dims = 1:5, nfeatures = 5)
VizDimLoadings(P5_T2.s, dims = 1:2, reduction = "pca")
DimPlot(P5_T2.s, reduction = "pca")
DimHeatmap(P5_T2.s, dims = 1, cells = 500, balanced = TRUE)

#dataset dimensionality
P5_T2.s<-JackStraw(P5_T2.s, num.replicate = 100)
P5_T2.s<-ScoreJackStraw(P5_T2.s, dims = 1:20)
JackStrawPlot(P5_T2.s, dims = 1:15)
ElbowPlot(P5_T2.s)

#cluster the cells
P5_T2.s<-FindNeighbors(P5_T2.s, dims = 1:20)
P5_T2.s<-FindClusters(P5_T2.s, resolution = 0.5)
head(Idents(P5_T2.s),5)

#UMAP reduction
P5_T2.s<-RunUMAP(P5_T2.s, dims = 1:20)
DimPlot(P5_T2.s, reduction = "umap", label = TRUE)
