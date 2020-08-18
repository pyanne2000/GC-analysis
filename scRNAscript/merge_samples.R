#integrate samples together
gastric_N.big.normalized<-merge(P3_N1.s, y=c(P3_N2.s, P5_N1.s, P5_N2.s), add.cell.ids=c("P3_N1", "P3_N2", "P5_N1", "P5_N2"), project = "Normal", merge.data = TRUE)
gastric_P.big.normalized<-merge(P3_P1.s, y=c(P3_P2.s, P5_P1.s, P5_P2.s), add.cell.ids=c("P3_P1", "P3_N2", "P5_P1", "P5_P2"), project = "PARI", merge.data = TRUE)
gastric_T.big.normalized<-merge(P3_T1.s, y=c(P3_T2.s, P3_T3.s, P5_T2.s), add.cell.ids=c("P3_T1", "P3_T2", "P5_T3", "P5_T1"), project = "Tumor", merge.data = TRUE)

gastric_N.big.normalized$stim <- "Normal"
gastric_P.big.normalized$stim <- "P"
gastric_T.big.normalized$stim <- "Tumor"
#merge together
gastric_P3P5<-merge(gastric_N.big.normalized, y=c(gastric_P.big.normalized, gastric_T.big.normalized), project = "GASTRIC12", merge.data = TRUE)

#a little bit Quality Control
gastric_P3P5[["percent.mt"]]<-PercentageFeatureSet(gastric_P3P5, pattern = "^MT-")
VlnPlot(gastric_P3P5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

plot1<-FeatureScatter(gastric_P3P5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2<-FeatureScatter(gastric_P3P5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#normalize data
gastric_P3P5<-NormalizeData(gastric_P3P5, normalization.method = "LogNormalize", scale.factor = 50000)

#feature selection
gastric_P3P5<-FindVariableFeatures(gastric_P3P5)
top10<-head(VariableFeatures(gastric_P3P5),10)

#scaling the data
all.genes<-rownames(gastric_P3P5)
gastric_P3P5<-ScaleData(gastric_P3P5, features = all.genes)

#dim reduction--PCA
gastric_P3P5<-RunPCA(gastric_P3P5, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
#gastric_P3P5<-JackStraw(gastric_P3P5, num.replicate = 1000)
gastric_P3P5<-ScoreJackStraw(gastric_P3P5, dims =1:100)
JackStrawPlot(gastric_P3P5, dims = 1:100)
ElowPlot(gastric_P3P5, ndims = 100)
DimHeatmap(gastric_P3P5, dims = c(1:3, 50:60), cells = 500, balanced = TRUE)

#clustering
gastric_P3P5<-FindNeighbors(gastric_P3P5, reduction = "pca", dims = 1:75)
gastric_P3P5<-FindClusters(gastric_P3P5, resolution = 0.4)
head(Idents(gastric_P3P5),5)

#Visualization(UMAP)
gastric_P3P5<-RunUMAP(gastric_P3P5,dims = 1:75)
DimPlot(gastric_P3P5, reduction = "umap", label = TRUE, group.by = "stim")

#Visualization(t-SNE)
gastric_P3P5<-RunTSNE(gastric_P3P5, dims = 1:75, nthreads = 4, max_iter =2000, check_duplicates=FALSE)

#Visualization(t-SNE vs. UMAP)
library(ggplot2)
p1<-DimPlot(gastric_P3P5, reduction = "tsne", label = TRUE) + ggtitle(label="t-SNE")
p2<-DimPlot(gastric_P3P5, reduction = "umap", label = TRUE) + ggtitle(label="UMAP")
p1<-AugmentPlot(plot = p1)
p2<-AugmentPlot(plot = p2)
(p1 + p2) & NoLegend()

#Visualization(batch effect)
library(cowplot)
DimPlot(gastric_P3P5, reduction = "umap", label = TRUE)
p1<-DimPlot(gastric_P3P5, reduction = "umap", group.by = gastric_P3P5)
p2<-DimPlot(gastric_P3P5, reduction = "umap", label = TRUE)
plot_grid(p1,p2)