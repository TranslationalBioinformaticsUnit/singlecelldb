##load McCarroll data
library(Seurat)
library(DropSeq.util)
df <- loadSparseDge("/ibex/projects/c2012/DataForSCDBProject/McCarroll/temp/GSE116470_F_GRCm38.81.P60SubstantiaNigra.raw.dge.txt")
subnigra <- CreateSeuratObject(raw.data = df, min.cells = 3, min.genes = 200, project = "P60SubstantiaNigra")
cerebellum <- loadSparseDge("/ibex/projects/c2012/DataForSCDBProject/McCarroll/temp/GSE116470_F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt")
cerebellum <- CreateSeuratObject(raw.data = cerebellum, min.cells = 3, min.genes = 200, project = "P60Cerebellum")
subnigra_cereb <- MergeSeurat(object1  = subnigra, object2 = cerebellum, add.cell.id1 ="SN", add.cell.id2 =  "CB", project="McCarroll")
save(subnigra, cerebellum, subnigra_cereb, file="/ibex/projects/c2012/DataForSCDBProject/McCarroll/rdataforMcCarroll/mccarroll_data.rda")
