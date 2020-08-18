SelectIntegrationFeatures(P3_T1.s, P3_T2.s, P3_T3.s, P3_P1.s, P3_P2.s, P3_N1.s, P3_N2.s)
gastric_P3.list.anchors<- FindIntegrationAnchors(object.list=list(P3_T1.s, P3_T2.s,P3_T3.s, P3_P1.s, P3_P2.s, P3_N1.s, P3_N2.s),anchor.features = 6000,dims = 1:20, k.anchor = 5, k.filter = 30)
gastric_P3.integrated<-IntegrateData(anchorset = gastric_P3.list.anchors, dims = 1:20)
SelectIntegrationFeatures(P5_T2.s, P5_P1.s, P5_P2.s, P5_N1.s, P5_N2.s)
gastric_P5.list.anchors<- FindIntegrationAnchors(object.list=list(P5_T2.s, P5_P1.s, P5_P2.s, P5_N1.s, P5_N2.s), dims = 1:20, k.anchor = 5, k.filter = 30)
gastric_P5.integrated<-IntegrateData(anchorset = gastric_P5.list.anchors, dims = 1:20)

#SelectIntegrationFeatures(gastric_P3.integrated, gastric_P5.integrated)

gastric_P3.integrated$stim <- "P3"
gastric_P5.integrated$stim <- "P5"

gastric_P3P5.list.anchors<-FindIntegrationAnchors(object.list = list(gastric_P3.integrated, gastric_P5.integrated),  anchor.features = 4000, dims = 1:20, k.anchor = 5, k.filter = 30)
gastric_P3P5.integrated<-IntegrateData(anchorset = gastric_P3P5.list.anchors, dims = 1:20)

#dataset dimensionality
gastric_P3.integrated<-ScaleData(gastric_P3.integrated, features = all.genes)
gastric_P5.integrated<-ScaleData(gastric_P5.integrated, features = all.genes)
gastric_P3P5.integrated<-ScaleData(gastric_P3P5.integrated, features = all.genes)

#dim reduction--PCA
gastric_P3.integrated<-RunPCA(gastric_P3.integrated, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
gastric_P5.integrated<-RunPCA(gastric_P5.integrated, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
gastric_P3P5.integrated<-RunPCA(gastric_P3P5.integrated, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)

#Visualization(UMAP)
gastric_P3.integrated<-RunUMAP(gastric_P3.integrated,dims = 1:75)
DimPlot(gastric_P3.integrated, reduction = "umap", label = TRUE)
gastric_P5.integrated<-RunUMAP(gastric_P5.integrated,dims = 1:75)
DimPlot(gastric_P5.integrated, reduction = "umap", label = TRUE)
gastric_P3P5.integrated<-RunUMAP(gastric_P3P5.integrated, dims = 1:75)
DimPlot(gastric_P3P5.integrated, reduction = "umap", label = TRUE, group.by = gastric_P3P5.integrated$stim)

#Visualization(t-SNE)
gastric_P3.integrated<-RunTSNE(gastric_P3.integrated, dims = 1:75, nthreads = 4, max_iter =2000, check_duplicates=FALSE)



#cluster the cells
gastric_P3.integrated<-FindNeighbors(gastric_P3.integrated, dims = 1:75)
gastric_P3.integrated<-FindClusters(gastric_P3.integrated, resolution = 0.4)
head(Idents(gastric_P3.integrated),5)

#UMAP reduction
gastric_P3.integrated<-RunUMAP(gastric_P3.integrated, dims = 1:75)
DimPlot(gastric_P3.integrated, reduction = "umap", label = TRUE)
