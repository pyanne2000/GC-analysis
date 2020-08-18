#doublet finder
library(DoubletFinder)
sweep.res.list_P3_N2<-paramSweep_v3(P3_N2.s, PCs = 1:10, sct = FALSE)
sweep.stats_P3_N2<-summarizeSweep(sweep.res.list_P3_N2, GT = FALSE)
bcmvn_P3_N2<-find.pK(sweep.stats_P3_N2)
View(bcmvn_P3_N2)
homotypic.prop<-modelHomotypic(P3_N2.s@meta.data$seurat_clusters)
nExp_poi<-round(0.075*length(rownames(P3_N2.s@meta.data)))
nExp_poi.adj<-round(nExp_poi*(1-homotypic.prop))
P3_N2.s<-doubletFinder_v3(P3_N2.s, PCs = 1:10, pN = 0.25, pK=0.16, nExp = 92, reuse.pANN = FALSE)
P3_N1.s<-doubletFinder_v3(P3_N2.s, PCs = 1:10, pN = 0.25, pK=0.16, nExp = 56, reuse.pANN = "pANN_0.25_0.16_92", sct = FALSE)
table(P3_N2.s@meta.data$DF.classifications_0.25_0.16_92)
#clean the doublet
P3_N2.s_clean<-subset(P3_N2.s, cells = rownames(P3_N2.s@meta.data)[which(P3_N2.s@meta.data$DF.classifications_0.25_0.16_92 == "Singlet")])
#re-cluster the library
DimPlot(P3_N2.s_clean, reduction = "umap",label = TRUE)
