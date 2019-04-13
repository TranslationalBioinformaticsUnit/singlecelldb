

##############################################################################

library(Matrix)
library(Seurat)
library('RcppCNPy')
setwd("/home/qus/SallyProjects/scNetworks/Linnarsson/Seurat/Data/")

# load(file = "linnarsson_matrix_t_unique.rda")  # stop here!!
# linnarsson_matrix <- t(cluster.mtx)
# save(linnarsson_matrix, file =  'linnarsson_matrix_t.rda')
# df <- data.frame(linnarsson_matrix)
# rownames(df) = make.names(rownames(df), unique=TRUE)
# linnarsson_matrix_unique <- data.matrix(df, rownames.force = TRUE)
# save(linnarsson_matrix_unique, file = 'linnarsson_matrix_t_unique.rda')


# linnarsson_metadata <- metadata

# linnarsson <- CreateSeuratObject(counts = linnarsson_matrix_unique,  min.cells = 3, project = "Mouse_Brain")

#################### add meta data###############################
# linnarsson$ClusterNames <-  metadata[,2]
# linnarsson$LeafOrder <-  metadata[,3]
# linnarsson$TaxonomyRank2 <-  metadata[,4]
# linnarsson$TaxonomyRank3 <-  metadata[,5]
# linnarsson$TaxonomyRank4 <-  metadata[,6]


################# read this rdata that is the seurat obeject######
save(linnarsson, file = 'linnarsson_Seurat3Object_with_meta.rda')

################### QC Standard pre-processing #######
# Visualize QC metrics as a violin plot
pdf("VlnFeature_RNA&nCount_RNA.pdf") 
VlnPlot(object = linnarsson, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("ScatterFeature_RNA&nCount_RNA.pdf") 
ScatterPlot <- FeatureScatter(object = linnarsson, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()


##################### add some filters#################
# linnarsson <- subset(x = linnarsson, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# min.features = 200, 
# min.cells = 3, 


################### Normalize Data r##################
# employ a global-scaling normalization method “LogNormalize” that 
# normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# Normalized values are stored in  linnarsson[["RNA"]]@data
linnarsson <- NormalizeData(object = linnarsson, normalization.method = "LogNormalize", scale.factor = 10000)
# Visualize QC metrics as a violin plot

VlnPlot(object = linnarsson, features = , cols= 2)



################# Find variable genes################
linnarsson <- FindVariableFeatures(linnarsson, mean.function = ExpMean, dispersion.function =  LogVMR
                            , do.plot = FALSE)
# linnarsson <- FindVariableFeatures(linnarsson, selection.method = "vst", loess.span = 0.3, clip.max = "auto", mean.function = ExpMean
#                            , dispersion.function = FastLogVMR, num.bin = 20
#                            , do.plot = FALSE)
# hv.genes <- head(rownames(HVFInfo(object = linnarsson)),5000)

# Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = linnarsson), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = linnarsson)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))


#mito.genes <- grep(pattern = "^mt-", x = rownames(x = GetAssayData(object = linnarsson)), value = TRUE)
#percent.mito <- Matrix::colSums(GetAssayData(object = linnarsson, slot = "counts")[mito.genes,])/Matrix::colSums(GetAssayData(object = linnarsson, slot = "counts"))
#linnarsson$name <- percent.mito
# linnarsson <- AddMetaData(object = linnarsson, metadata = percent.mito, col.name = "percent.mito")
# linnarsson <- ScaleData(object = linnarsson, genes.use = hv.genes, display.progress = FALSE,
#                  vars.to.regress = "percent.mito", do.par = TRUE, num.cores = 20)

all.genes <- rownames(x = linnarsson)
hv.genes <- VariableFeatures(object = linnarsson)
top1000 <- head(x = VariableFeatures(object = linnarsson), 1000)


linnarsson <- ScaleData(object = linnarsson, features = all.genes, display.progress = FALSE,
                 do.par = TRUE, num.cores = 20)
#linnarsson <- ScaleData(object = linnarsson, features = hv.genes, display.progress = FALSE,
#                        do.par = TRUE, num.cores = 20)

linnarsson <- RunPCA(object = linnarsson, pc.genes = hv.genes, npcs = 100, do.print = TRUE,
              ndims.print = 1:5, nfeatures.print = 5)

VizDimLoadings(object = linnarsson, dims = 1:2, reduction = "pca")
DimPlot(object = linnarsson, reduction = "pca")
ElbowPlot(object = linnarsson, ndims = 100, reduction = "pca")
# PCHeatmap(linnarsson, pc.use = c(1:3, 70:75), cells.use = 500, do.balanced = TRUE)
# DoHeatmap(object = pbmc, features = heatmap_markers)
DimHeatmap(object = linnarsson, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = linnarsson, dims = c(1:3, 70:75), reduction = "pca", cells = 500)

linnarsson <- JackStraw(object = linnarsson,reduction = "pca", dims = 100, num.replicate = 100)
linnarsson <- ScoreJackStraw(object = linnarsson, dims = 1:100)
JackStrawPlot(object = linnarsson, dims = 1:100)



linnarsson <- FindNeighbors(object = linnarsson, resolution = 1, reduction = "pca", dims = 1:75, nn.eps = 0)
# linnarsson <- FindClusters(object = linnarsson)
linnarsson <- FindClusters(object = linnarsson, algorithm=1, reduction.type = "pca", dims.use = 1:75, graph.name = NULL,resolution = 3, 
                    save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)



##################### 1 FIt-SNE and UMAP #########################
# source('<path to file>', chdir=T)
 

linnarsson <- RunTSNE(object = linnarsson, reduction.use = "pca", dims.use = 1:75, tsne.method = "FIt-SNE", 
               nthreads = 4, reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "PATH/TO//fast_tsne", 
               max_iter = 2000)
linnarsson <- RunUMAP(object = linnarsson, reduction.use = "pca", dims = 1:75, min_dist = 0.75)


library(cowplot)
p1 <- DimPlot(object = linnarsson, reduction = "FItSNE", no.legend = TRUE, do.return = TRUE, 
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("FIt-SNE") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(object = linnarsson, reduction = "umap", no.legend = TRUE, do.return = TRUE, 
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
pdf("Seurat_75pcs_FItSNE.pdf") 
plot_grid(p1)
dev.off()
pdf("Seurat_75pcs_UMAP.pdf") 
plot_grid(p2)
dev.off()

FeaturePlot(linnarsson, c("S100a9", "Sftpc"), reduction.use = "FItSNE", Dark.Theme = TRUE, 
            pt.size = 0.1, vector.friendly = TRUE)


####################### 2. FIt-SNE and UMAP #############################
### set nn.eps to 0 
linnarsson <- FindClusters(object = linnarsson, reduction.type = "pca", dims.use = 1:80, graph.name = NULL,resolution = 3, 
                           save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)


# source('<path to file>', chdir=T)
linnarsson <- RunTSNE(object = linnarsson, reduction.use = "pca", dims.use = 1:80, tsne.method = "FIt-SNE", 
                      nthreads = 4, reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "PATH/To/fast_tsne", 
                      max_iter = 2000)
linnarsson <- RunUMAP(object = linnarsson, reduction.use = "pca", dims = 1:80, min_dist = 0.75)

library(cowplot)
p3 <- DimPlot(object = linnarsson, reduction = "FItSNE", no.legend = TRUE, do.return = TRUE, 
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("FIt-SNE_accurateNeighbors") + theme(plot.title = element_text(hjust = 0.5))
p4 <- DimPlot(object = linnarsson, reduction = "umap", no.legend = TRUE, do.return = TRUE, 
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("UMAP_accurateNeighbors") + theme(plot.title = element_text(hjust = 0.5))
pdf("SeuratClustering_75pcs_acc.pdf") 
plot_grid(p3,p4)
dev.off()


######################### 3. tSNE ##################################
### set nn.eps to 0 
linnarsson <- FindClusters(object = linnarsson, reduction.type = "pca", dims.use = 1:80, graph.name = NULL,resolution = 3, 
                           save.SNN = TRUE, n.start = 100, nn.eps = 0, print.output = FALSE)

linnarsson <- RunTSNE(object = linnarsson, dims = 1:50)
linnarsson <- RunUMAP(object = linnarsson, reduction.use = "pca", dims = 1:50, min_dist = 0.75)

# pdf("SeuratClustering_50pcs.pdf") 
p5 <- DimPlot(object = linnarsson, no.legend = TRUE, do.return = TRUE, 
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("t-SNE") + theme(plot.title = element_text(hjust = 0.5))
p6 <- DimPlot(object = linnarsson, reduction = "umap", no.legend = TRUE, do.return = TRUE, 
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))

pdf("Seurat_50pcs_Acc_FItSNE.pdf") 
plot_grid(p5)
dev.off()
pdf("Seurat_50pcs_Acc_UMAP.pdf") 
plot_grid(p6)
dev.off()





