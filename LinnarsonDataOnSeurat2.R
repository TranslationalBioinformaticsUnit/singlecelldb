### Reading Linnarson file
#install.packages("devtools",lib="/ibex/scratch/thimmamp/R_local_library" )
#devtools::install_github(repo = "hhoeflin/hdf5r",lib="/ibex/scratch/thimmamp/R_local_library" )
#devtools::install_github(repo = "mojaveazure/loomR", ref = "develop", lib="/ibex/scratch/thimmamp/R_local_library/")
#install.packages("loomR", lib="/ibex/scratch/thimmamp/R_local_library")
library(loomR, lib.loc="/ibex/scratch/thimmamp/R_local_library/")
library(Seurat)
library(Matrix)
#cluster.scdata <- connect("/ibex/projects/c2012/DataForSCDBProject/Linnarson/l5_all.loom", mode="rb")
#load("/ibex/projects/c2012/NetworksProject/Data/loom_to_seurat.rda")
#cluster.scdata <- connect("/ibex/scratch/thimmamp/DataForSCDBProject/l5_all.loom", "r+")
cluster.scdata <- connect("~/l5_all.loom", mode="r")

#cluster.scdata

cluster.mtx <- cluster.scdata$matrix
cluster.mtx <-as.matrix(cluster.mtx[,])
save(cluster.mtx, file = "/ibex/projects/c2012/NetworksProject/Data/Data_untrans.rda")
cellID <- cluster.scdata[["col_attrs/CellID"]][]
gene_names <-cluster.scdata[["row_attrs/Gene"]][]

##########finding duplication ID################################
# ENID <- cluster.scdata[["row_attrs/Accession"]][]
# ENID_dup <- as.vector(ENID[group %in% group3])
# gene_dup <- as.vector(group[group %in% group3])
# ENID_gene_dup <- as.data.frame(t(rbind(ENID_dup, gene_dup)))
# colnames(ENID_gene_dup) <- c("Acession", "Gene")
#save(ENID_gene_dup,m3,file ="/ibex/projects/c2012/NetworksProject/Data/Duplicated_ID_genes.rda")

colnames(cluster.mtx) <- gene_names
rownames(cluster.mtx) <- cellID

mouse.mtx <- t(cluster.mtx)
#mtx_forseurat <- Matrix(cluster.mtx, sparse = TRUE) 

Clusters <- cluster.scdata[["col_attrs/Clusters"]][]
ClusterNames <-cluster.scdata[["col_attrs/ClusterName"]][]
LeafOrder <- cluster.scdata[["col_attrs/LeafOrder"]][]
TaxonomyRank2 <-cluster.scdata[["col_attrs/TaxonomyRank2"]][]
TaxonomyRank3 <- cluster.scdata[["col_attrs/TaxonomyRank3"]][]
TaxonomyRank4 <- cluster.scdata[["col_attrs/TaxonomyRank4"]][]
attr <- cbind(Clusters, ClusterNames,LeafOrder,TaxonomyRank2,TaxonomyRank3,TaxonomyRank4)
#cellID <- cluster.scdata[["col_attrs/CellID"]][]
row.names(attr) <- cellID
attr <- as.data.frame(attr)
#write.table(attr, file = "/ibex/projects/c2012/NetworksProject/Linnarson_metadata.txt", sep="\t")
#save(mouse.mtx,attr, file = "/ibex/projects/c2012/NetworksProject/Data/Data_trans.rda")
mouse.mtx[1:5, 1:5]
sampled_mouse_mtx <- mouse.mtx[, 1:20000]
#mca <- CreateSeuratObject(raw.data =mouse.mtx,meta.data = attr[1:20000,], project = "Mousebrain")
mca <- CreateSeuratObject(raw.data =sampled_mouse_mtx,meta.data = attr[1:20000,], project = "Mousebrain")

mca <- MakeSparse(mca)
#mca <- SubsetData(mca, cells.use = rownames(mca@meta.data[!is.na(mca@meta.data$ClusterID), 
#                                                          ]), do.clean = TRUE)
#save(mca, cluster.mtx, file ="/ibex/scratch/thimmamp/loom_to_seurat.rda")
#save(attr, mouse.mtx, mca, file="/ibex/scratch/thimmamp/loom_to_seurat.rda")
#load("/ibex/scratch/thimmamp/loom_to_seurat.rda")

load("/ibex/projects/c2012/NetworksProject/Data/mca_seurat_sparse.rda")

# aggragates vlues from the same gene (with same rowname)
group <- row.names(mca_sparse@data)
group3 <-group[duplicated(group)]
length(group3)
matrix <- mca_sparse@raw.data
matrix2 <- matrix[-which(row.names(matrix)%in% group3),]
matrix3 <-matrix[which(row.names(matrix)%in% group3),]
m3 <- as.matrix(matrix3)
m4 <- rowsum(m3,row.names(m3))
dim(m4)
matrix4 <- as(m4, "dgCMatrix")
matrix_merged <- rbind(matrix2,matrix4)
dim(matrix_merged)
length(unique(row.names(mca_sparse@raw.data)))

mca_sparse_merged <- mca_sparse
mca_sparse_merged@raw.data <- matrix_merged
#save(mca_sparse_merged, file ="/ibex/projects/c2012/NetworksProject/Data/mca_sparse_merged.rda") 


load("/ibex/projects/c2012/NetworksProject/Data/mca_sparse_merged.rda") 
mca_norm <- NormalizeData(object = mca_sparse_merged, normalization.method = "LogNormalize", scale.factor = 10000)
mca_norm <- FindVariableGenes(object = mca_norm)
mca_norm <- ScaleData(object = mca_norm, genes.use = mca_norm@var.genes)
mca_norm <- RunPCA(object = mca_norm, pcs.compute =50, ndims.print = 1:5, nfeatures.print = 5)
PCElbowPlot(object = mca_norm,num.pc = 50)
save(mca_norm, file ="/ibex/projects/c2012/NetworksProject/Data/mca_normalized_PCA.rda") 

mca <- FindClusters(object = mca_norm, reduction.type = "pca", dims.use = 1:50, resolution = 3, 
                    save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE)

mca <-  RunTSNE(mca, reduction.use = "pca", dims.use = 1:50, do.fast = T,perplexity=100,check_duplicates = FALSE)
msca <- RunUMAP(object = mca, reduction.use = "pca", dims.use = 1:25, min_dist = 0.75)
save(msca,file = "/ibex/projects/c2012/NetworksProject/Data/mca__PCA50_UMAP25.rds")
load("/ibex/projects/c2012/NetworksProject/Data/mca.rds")
install.packages("reticulate", lib="/ibex/scratch/thimmamp/R_local_library")

load("/ibex/projects/c2012/NetworksProject/Data/mca__PCA50_UMAP25.rds")
#library(reticulate, lib ="~/R/x86_64-pc-linux-gnu-library/3.4")

############################################remove duplicated cols#################
df <-msca@cell.names
dup_cells <- colnames(msca@data)[duplicated(colnames(msca@data))]

dup_index <- data.frame()
for (i in dup_cells) {
  index <- grep(pattern = i, df)
  cat(index,"\n")
  dup_index <- rbind(dup_index, index)
}
colnames(dup_index) =c("cell1", "cell2")

# Check if all duplicated cols are the same 
# > sum(msca@raw.data[,as.vector(dup_index[,1])] == msca@raw.data[,as.vector(dup_index[,2])])
# [1] 3296094
# > length(msca@raw.data[,as.vector(dup_index[,1])] == msca@raw.data[,as.vector(dup_index[,2])])
# [1] 3296094
# > 3296094/nrow(msca@raw.data)
# [1] 118
# > length(dup_index[,1])
# [1] 118

mca_dedup_col <- SubsetData(object = msca,cells.use = msca@cell.names[-dup_index[,1]])
save(mca_dedup_col,file="/ibex/projects/c2012/NetworksProject/Data/mca__PCA50_UMAP25_dedup.rds")

load("/ibex/projects/c2012/NetworksProject/Data/mca__PCA50_UMAP25_dedup.rds")

DeDup_T4 <- as.data.frame(TaxonomyRank4[-dup_index[,1]])
rownames(DeDup_T4) <- colnames(mca_dedup_col@data)

DeDup_T3 <- as.data.frame(TaxonomyRank3[-dup_index[,1]])
rownames(DeDup_T3) <- colnames(mca_dedup_col@data)

mca_T4 <- AddMetaData(mca_dedup_col,DeDup_T4)
colnames(mca_T4@meta.data)[length(colnames(mca_T4@meta.data))] <- "Rank4"
mca_T4 <- AddMetaData(mca_T4,DeDup_T3, col.name = Rank3)
colnames(mca_T4@meta.data)[length(colnames(mca_T4@meta.data))] <- "Rank3"

mca_T4<- StashIdent(object = mca_T4, save.name = "CellType")
mca_T4<- SetAllIdent(object = mca_T4, id = "Rank3")
###############################################plot###############################################
library(cowplot)
colors <- colorRampPalette(c("blue", "red"))(39)
p1 <- DimPlot(object = mca_T4, reduction.use = "pca", no.legend = F, do.return = TRUE, 
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5))
p2 <- DimPlot(object = mca_T4, reduction.use = "umap", no.legend = F, do.return = TRUE, 
              vector.friendly = TRUE, pt.size = 0.1) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
plot_grid(p1, p2)
plot(p1)
plot(p1)
save(mca_T4, file="/ibex/projects/c2012/NetworksProject/Data/mca__PCA50_UMAP25_T34.rds")
mca_final <- RunUMAP(object = mca_T4, reduction.use = "pca", dims.use = 1:50, min_dist = 0.75)
