BiocManager::install("AUCell")

library(dplyr)
library(Seurat)
library(patchwork)
library("org.Hs.eg.db")

library(AUCell)
library(GSEABase)

# Load dataset
melanoma.data <- read.table("D:/McGill/projects/project peroxisome/GSE116237/GSE116237_scRNAseq_expressionMatrix.txt", sep = ",", head = T)
rownames(melanoma.data) <- melanoma.data$X
melanoma.data <- melanoma.data[, -1]

melanoma <- CreateSeuratObject(counts = melanoma.data, project = "melanoma", min.cells = 3, min.features = 200)

# Normalizing the data  
melanoma <- NormalizeData(melanoma, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection)
melanoma <- FindVariableFeatures(melanoma, selection.method = "vst", nfeatures = 200)

# Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(melanoma), 10)

#Scaling the data standardization
#This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
all.genes <- rownames(melanoma)
melanoma <- ScaleData(melanoma, features = all.genes)

melanoma <- RunPCA(melanoma, npcs = 10, features = VariableFeatures(object = melanoma))
# Examine and visualize PCA results a few different ways
print(melanoma[["pca"]], dims = 1:5, nfeatures = 5)

melanoma <- FindNeighbors(melanoma, dims = 1 : 10)
melanoma <- FindClusters(melanoma, resolution = 1)

# choose one
# TSNE
melanoma = RunTSNE(object = melanoma, perplexity = 30, dims.use = 1 : 10, do.fast = TRUE)
Idents(melanoma) <- c(rep("T0", 172), rep("T4", 155), rep("T28", 199), rep("Tres", 148))
DimPlot(melanoma, reduction = "tsne", pt.size = 1.5, shape.by = NULL)
# UMAP
# melanoma = RunUMAP(object = melanoma, reduction = "pca", dims = 1 : 50)
# DimPlot(melanoma, reduction = "umap")

T0 <- subset(melanoma, idents = "T0")
T4 <- subset(melanoma, idents = "T4")
T28 <- subset(melanoma, idents = "T28")
Tres <- subset(melanoma, idents = "Tres")

saveRDS(melanoma, file = "D:/McGill/projects/project peroxisome/melanoma.rds")
saveRDS(T0, file = "D:/McGill/projects/project peroxisome/T0.rds")
saveRDS(T4, file = "D:/McGill/projects/project peroxisome/T4.rds")
saveRDS(T28, file = "D:/McGill/projects/project peroxisome/T28.rds")
saveRDS(Tres, file = "D:/McGill/projects/project peroxisome/Tres.rds")

melanoma <- readRDS("D:/McGill/projects/project peroxisome/GSE116237/melanoma.rds")

FeaturePlot(object = melanoma, features = "ENSG00000135218", reduction = "tsne", split.by = "ident") # CD36

# MITF SOX10 CD36 NGFR AXL
VlnPlot(object = melanoma, features = c("ENSG00000187098", "ENSG00000100146", "ENSG00000135218", "ENSG00000064300", "ENSG00000167601"), pt.size = 0)
FeaturePlot(object = melanoma, features = c("ENSG00000187098", "ENSG00000100146", "ENSG00000135218", "ENSG00000064300", "ENSG00000167601"), reduction = "tsne", pt.size = 2)

FeaturePlot(object = melanoma, features = "ENSG00000187098", reduction = "tsne", pt.size = 2, min.cutoff = 2)
FeaturePlot(object = melanoma, features = "ENSG00000115648", reduction = "tsne", pt.size = 2, min.cutoff = 1)


# UGCG GBA PEX3 AGPS
VlnPlot(object = melanoma, features = c("ENSG00000148154", "ENSG00000177628", "ENSG00000034693", "ENSG00000018510"), pt.size = 0)
FeaturePlot(object = melanoma, features = c("ENSG00000148154", "ENSG00000177628", "ENSG00000034693", "ENSG00000018510"), reduction = "tsne", pt.size = 2)

markers = FindAllMarkers(object = melanoma, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

id2symbol <- toTable(org.Hs.egSYMBOL)
id2ensembl <- toTable(org.Hs.egENSEMBL)
id <- id2ensembl$gene_id[match(markers$gene, id2ensembl$ensembl_id)]
markers$symbol <- id2symbol$symbol[match(id, id2symbol$gene_id)]
write.csv(markers, "D:/McGill/projects/project peroxisome/GSE116237/overall/melanoma.csv")

# dimnames(melanoma[top10$gene])[[1]] <- symbol
top10 <- markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(object = melanoma, features = symbol)

# AUCell
signature <- read.csv("D:/McGill/projects/project peroxisome/GSE116237/signature.csv")
colnames(signature) <- c("mitosis", "pigmentation", "invasion", "neuro", "immune", "MITF-targets", "hypometabolic", "SMC")
cells_rankings <- AUCell_buildRankings(melanoma@assays[["RNA"]]@data)
pigmentation <- GeneSet(signature$pigmentation[1 : 15], setName = "pigmentation")
invasion <- GeneSet(signature$invasion[1 : 49], setName = "invasion")
NCSC <- GeneSet(signature$neuro[1 : 37], setName = "NCSC")
SMC <- GeneSet(signature$SMC[1 : 14], setName = "SMC")
pigmentation_AUC <- AUCell_calcAUC(pigmentation, cells_rankings)
invasion_AUC <- AUCell_calcAUC(invasion, cells_rankings)
NCSC_AUC <- AUCell_calcAUC(NCSC, cells_rankings)
SMC_AUC <- AUCell_calcAUC(SMC, cells_rankings)
melanoma@assays[["RNA"]]@data <- rbind(melanoma@assays[["RNA"]]@data,
                                       pigmentation_AUC@assays@data@listData[["AUC"]],
                                       invasion_AUC@assays@data@listData[["AUC"]],
                                       NCSC_AUC@assays@data@listData[["AUC"]],
                                       SMC_AUC@assays@data@listData[["AUC"]])

FeaturePlot(object = melanoma, features = "pigmentation", reduction = "tsne", pt.size = 2, min.cutoff = 0.5)
FeaturePlot(object = melanoma, features = "ENSG00000167601", reduction = "tsne", pt.size = 2, min.cutoff = 0.5) # invasion AXL
FeaturePlot(object = melanoma, features = "ENSG00000135218", reduction = "tsne", pt.size = 2, min.cutoff = 0.5) # SMC CD36

FeaturePlot(object = melanoma, features = "pigmentation", reduction = "tsne", pt.size = 2, min.cutoff = 0.4)
FeaturePlot(object = melanoma, features = "invasion", reduction = "tsne", pt.size = 2, min.cutoff = 0.1)
FeaturePlot(object = melanoma, features = "NCSC", reduction = "tsne", pt.size = 2, min.cutoff = 0.1)
FeaturePlot(object = melanoma, features = "SMC", reduction = "tsne", pt.size = 2, min.cutoff = 0.05)

pigmentation_id <- names(which(melanoma@assays[["RNA"]]@data["pigmentation",] > 0.5))
invasion_id <- names(which(melanoma@assays[["RNA"]]@data["invasion",] > 0.2))
NCSC_id <- names(which(melanoma@assays[["RNA"]]@data["NCSC",] > 0.3))
SMC_id <- names(which(melanoma@assays[["RNA"]]@data["SMC",] > 0.1))

T4_pigm <- na.omit(match(pigmentation_id, colnames(T4))) + 172
T4_inv <- na.omit(match(invasion_id, colnames(T4))) + 172
T4_NCSC <- na.omit(match(NCSC_id, colnames(T4))) + 172
T4_SMC <- na.omit(match(SMC_id, colnames(T4))) + 172

T28_pigm <- na.omit(match(pigmentation_id, colnames(T28))) + 327
T28_inv <- na.omit(match(invasion_id, colnames(T28))) + 327
T28_NCSC <- na.omit(match(NCSC_id, colnames(T28))) + 327
T28_SMC <- na.omit(match(SMC_id, colnames(T28))) + 327

identity <- c(rep("T0", 172), rep("T4_melanocytic", 155), rep("T28_melanocytic", 199), rep("Tres", 148))
identity[T4_pigm] <- rep("T4_pigmentation", length(T4_pigm))
identity[T4_inv] <- rep("T4_invasion", length(T4_inv))
identity[T4_SMC] <- rep("T4_SMC", length(T4_SMC))
identity[T4_NCSC] <- rep("T4_NCSC", length(T4_NCSC))
identity[T28_pigm] <- rep("T28_pigmentation", length(T28_pigm))
identity[T28_inv] <- rep("T28_invasion", length(T28_inv))
identity[T28_SMC] <- rep("T28_SMC", length(T28_SMC))
identity[T28_NCSC] <- rep("T28_NCSC", length(T28_NCSC))

melanoma <- readRDS("D:/McGill/projects/project peroxisome/GSE116237/melanoma.rds")
Idents(melanoma) <- identity

VlnPlot(object = subset(melanoma, idents = c("T4_pigmentation", "T4_melanocytic", "T4_invasion", "T4_SMC", "T4_NCSC")),
        features = c("ENSG00000148154", "ENSG00000177628", "ENSG00000034693", "ENSG00000018510"), pt.size = 0)
VlnPlot(object = subset(melanoma, idents = c("T28_pigmentation", "T28_melanocytic", "T28_invasion", "T28_SMC", "T28_NCSC")),
        features = c("ENSG00000148154", "ENSG00000177628", "ENSG00000034693", "ENSG00000018510"), pt.size = 0)

VlnPlot(object = melanoma,
        features = c("pigmentation", "invasion", "SMC", "NCSC"), pt.size = 0)
DimPlot(melanoma)

saveRDS(melanoma, file = "D:/McGill/projects/project peroxisome/melanoma_pop.rds")
