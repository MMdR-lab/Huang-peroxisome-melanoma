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

T0 <- subset(melanoma, idents = "T0")
T4 <- subset(melanoma, idents = "T4")
T28 <- subset(melanoma, idents = "T28")
Tres <- subset(melanoma, idents = "Tres")

# Figure S5A
saveRDS(melanoma, file = "D:/McGill/projects/project peroxisome/melanoma.rds")
saveRDS(T0, file = "D:/McGill/projects/project peroxisome/T0.rds")
saveRDS(T4, file = "D:/McGill/projects/project peroxisome/T4.rds")
saveRDS(T28, file = "D:/McGill/projects/project peroxisome/T28.rds")
saveRDS(Tres, file = "D:/McGill/projects/project peroxisome/Tres.rds")

# data processing done!
# ------I am a separation line-------------------------------------------------------------------------------------------------


# AUCell
melanoma <- readRDS("D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/melanoma.rds")
T0 <- readRDS("D:/McGill/projects/peroxisome/GSE116237/T0.rds")
T4 <- readRDS("D:/McGill/projects/peroxisome/GSE116237/T4.rds")
T28 <- readRDS("D:/McGill/projects/peroxisome/GSE116237/T28.rds")
Tres <- readRDS("D:/McGill/projects/peroxisome/GSE116237/Tres.rds")

signature <- read.csv("D:/McGill/projects/peroxisome/Fan-peroxisome/RAWDATA/signature.csv")
colnames(signature) <- c("mitosis", "pigmentation", "invasion", "neuro", "immune", "MITF-targets", "hypometabolic", "SMC", "peroxisome")
cells_rankings <- AUCell_buildRankings(melanoma@assays[["RNA"]]@data)

pigmentation <- GeneSet(signature$pigmentation[1 : 15], setName = "pigmentation")
invasion <- GeneSet(signature$invasion[1 : 49], setName = "invasion")
NCSC <- GeneSet(signature$neuro[1 : 37], setName = "NCSC")
SMC <- GeneSet(signature$SMC[1 : 14], setName = "SMC")
metabolism <- GeneSet(signature$hypometabolic[1 : 27], setName = "metabolism")
peroxisome <- GeneSet(signature$peroxisome[1 : 11], setName = "peroxisome")

pigmentation_AUC <- AUCell_calcAUC(pigmentation, cells_rankings)
invasion_AUC <- AUCell_calcAUC(invasion, cells_rankings)
NCSC_AUC <- AUCell_calcAUC(NCSC, cells_rankings)
SMC_AUC <- AUCell_calcAUC(SMC, cells_rankings)
metabolism_AUC <- AUCell_calcAUC(metabolism, cells_rankings)
peroxisome_AUC <- AUCell_calcAUC(peroxisome, cells_rankings)

melanoma@assays[["RNA"]]@data <- rbind(melanoma@assays[["RNA"]]@data,
                                       pigmentation_AUC@assays@data@listData[["AUC"]],
                                       invasion_AUC@assays@data@listData[["AUC"]],
                                       NCSC_AUC@assays@data@listData[["AUC"]],
                                       SMC_AUC@assays@data@listData[["AUC"]],
                                       metabolism_AUC@assays@data@listData[["AUC"]],
                                       peroxisome_AUC@assays@data@listData[["AUC"]])

# cell signature plot
FeaturePlot(object = melanoma, features = "pigmentation", reduction = "tsne", pt.size = 2, min.cutoff = 0.4)
FeaturePlot(object = melanoma, features = "invasion", reduction = "tsne", pt.size = 2, min.cutoff = 0.1)
FeaturePlot(object = melanoma, features = "NCSC", reduction = "tsne", pt.size = 2, min.cutoff = 0.1)
FeaturePlot(object = melanoma, features = "SMC", reduction = "tsne", pt.size = 2, min.cutoff = 0.05)
FeaturePlot(object = melanoma, features = "peroxisome", reduction = "tsne", pt.size = 2, min.cutoff = 0.05)

# cell state definition
pigmentation_id <- names(which(melanoma@assays[["RNA"]]@data["pigmentation",] > 0.5))
invasion_id <- names(which(melanoma@assays[["RNA"]]@data["invasion",] > 0.1999))
NCSC_id <- names(which(melanoma@assays[["RNA"]]@data["NCSC",] > 0.3))
SMC_id <- names(which(melanoma@assays[["RNA"]]@data["SMC",] > 0.10006))

T0_pigm <- na.omit(match(pigmentation_id, colnames(T0)))
T0_inv <- na.omit(match(invasion_id, colnames(T0)))
T0_NCSC <- na.omit(match(NCSC_id, colnames(T0)))
T0_SMC <- na.omit(match(SMC_id, colnames(T0)))
differ <- union(union(T0_pigm, T0_inv), union(T0_NCSC, T0_SMC))
T0_other <- setdiff(seq(1, 172), differ)

T4_pigm <- na.omit(match(pigmentation_id, colnames(T4))) + 172
T4_inv <- na.omit(match(invasion_id, colnames(T4))) + 172
T4_NCSC <- na.omit(match(NCSC_id, colnames(T4))) + 172
T4_SMC <- na.omit(match(SMC_id, colnames(T4))) + 172
differ <- union(union(T4_pigm, T4_inv), union(T4_NCSC, T4_SMC))
T4_other <- setdiff(seq(173, 327), differ)

T28_pigm <- na.omit(match(pigmentation_id, colnames(T28))) + 327
T28_inv <- na.omit(match(invasion_id, colnames(T28))) + 327
T28_NCSC <- na.omit(match(NCSC_id, colnames(T28))) + 327
T28_SMC <- na.omit(match(SMC_id, colnames(T28))) + 327
differ <- union(union(T28_pigm, T28_inv), union(T28_NCSC, T28_SMC))
T28_other <- setdiff(seq(328, 526), differ)

Tres_pigm <- na.omit(match(pigmentation_id, colnames(Tres))) + 526
Tres_inv <- na.omit(match(invasion_id, colnames(Tres))) + 526
Tres_NCSC <- na.omit(match(NCSC_id, colnames(Tres))) + 526
Tres_SMC <- na.omit(match(SMC_id, colnames(Tres))) + 526
differ <- union(union(Tres_pigm, Tres_inv), union(Tres_NCSC, Tres_SMC))
Tres_other <- setdiff(seq(527, 674), differ)

identity <- c(rep("T0_other", 172), rep("T4_other", 155), rep("T28_other", 199), rep("Tres_other", 148))

identity[T0_pigm] <- rep("T0_pigmentation", length(T0_pigm))
identity[T0_inv] <- rep("T0_invasion", length(T0_inv))
identity[T0_SMC] <- rep("T0_SMC", length(T0_SMC))
identity[T0_NCSC] <- rep("T0_NCSC", length(T0_NCSC))

identity[T4_pigm] <- rep("T4_pigmentation", length(T4_pigm))
identity[T4_inv] <- rep("T4_invasion", length(T4_inv))
identity[T4_SMC] <- rep("T4_SMC", length(T4_SMC))
identity[T4_NCSC] <- rep("T4_NCSC", length(T4_NCSC))

identity[T28_pigm] <- rep("T28_pigmentation", length(T28_pigm))
identity[T28_inv] <- rep("T28_invasion", length(T28_inv))
identity[T28_SMC] <- rep("T28_SMC", length(T28_SMC))
identity[T28_NCSC] <- rep("T28_NCSC", length(T28_NCSC))

identity[Tres_pigm] <- rep("Tres_pigmentation", length(Tres_pigm))
identity[Tres_inv] <- rep("Tres_invasion", length(Tres_inv))
identity[Tres_SMC] <- rep("Tres_SMC", length(Tres_SMC))
identity[Tres_NCSC] <- rep("Tres_NCSC", length(Tres_NCSC))

melanoma <- readRDS("D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/melanoma.rds")
Idents(melanoma) <- identity

# Figure S5C
# melanoma_pop: with cell state definition but without AUCell data
saveRDS(melanoma, file = "D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/melanoma_pop.rds")
write.csv(Idents(melanoma), "D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/AUCell_idents.csv")

# melanoma_pop_AUCell: with cell state definition and AUCell data
saveRDS(melanoma, file = "D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/melanoma_pop_AUCell.rds")


# ------I am a separation line--------------------------------------------------------------------------------------------------------------------------------


# exprs matrix
melanoma <- readRDS("D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/melanoma_pop.rds")
identity <- Idents(melanoma)


T4 <- subset(melanoma, idents = c("T4_other", "T4_invasion", "T4_NCSC", "T4_pigmentation", "T4_SMC"))

pick <- c("ENSG00000187098", "ENSG00000135218", "ENSG00000167601", "ENSG00000064300", "ENSG00000117528", # MITF CD36 AXL NGFR ABCD3
          "ENSG00000116171", "ENSG00000127980", "ENSG00000164751", "ENSG00000034693", "ENSG00000162928", # SCP2 PEX1 PEX2 PEX3 PEX13
          "ENSG00000018510", "ENSG00000197601", "ENSG00000145782", "ENSG00000087470", "ENSG00000109819" # AGPS FAR1 ATG12 DNM1L PPARGC1A
          )

expr_matrix <- melanoma@assays[["RNA"]]@data[pick,]

expr_matrix@Dimnames[[1]] <- c("MITF", "CD36", "AXL", "NGFR", "ABCD3",
                               "SCP2", "PEX1", "PEX2", "PEX3", "PEX13",
                               "AGPS", "FAR1", "ATG12", "DNM1L", "PPARGC1A"
)

expr_matrix@Dimnames[[2]] <- as.character(Idents(melanoma))
order <- c(which(colnames(expr_matrix) == "T4_other"),
           which(colnames(expr_matrix) == "T4_invasion"),
           which(colnames(expr_matrix) == "T4_NCSC"),
           which(colnames(expr_matrix) == "T4_pigmentation"),
           which(colnames(expr_matrix) == "T4_SMC")
)
expr_matrix <- expr_matrix[, order]

# Figure 4A
write.csv(expr_matrix, "D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/expr_matrix_v3.csv")



# ------I am a separation line--------------------------------------------------------------------------------------------------------------------------------




# CD36
melanoma <- readRDS("D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/melanoma.rds")
# Figure S5B
FeaturePlot(object = melanoma, features = "ENSG00000135218", reduction = "tsne", pt.size = 2, min.cutoff = 2,
            cols = c("lightgrey", "red"))
FeaturePlot(object = melanoma, features = "SMC", reduction = "tsne", pt.size = 2, min.cutoff = 0.05,
            cols = c("lightgrey", "blue"))

# separate cell to be CD36+ or CD36-
CD36_id <- names(which(melanoma@assays[["RNA"]]@data["ENSG00000135218",] > 2.2))

T0_CD36 <- na.omit(match(CD36_id, colnames(T0)))
T4_CD36 <- na.omit(match(CD36_id, colnames(T4))) + 172
T28_CD36 <- na.omit(match(CD36_id, colnames(T28))) + 327
Tres_CD36 <- na.omit(match(CD36_id, colnames(Tres))) + 526

identity <- c(rep("T0_CD36neg", 172), rep("T4_CD36neg", 155), rep("T28_CD36neg", 199), rep("Tres_CD36neg", 148))

identity[T0_CD36] <- rep("T0_CD36pos", length(T0_CD36))
identity[T4_CD36] <- rep("T4_CD36pos", length(T4_CD36))
identity[T28_CD36] <- rep("T28_CD36pos", length(T28_CD36))
identity[Tres_CD36] <- rep("Tres_CD36pos", length(Tres_CD36))

Idents(melanoma) <- identity

# melanoma_pop_CD36: with cell state CD36+/-
saveRDS(melanoma, file = "D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/melanoma_pop_CD36.rds")

# exprs matrix
melanoma <- readRDS("D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/melanoma_pop_CD36.rds")

pick <- c(# "metabolism",
          "ENSG00000148154", "ENSG00000177628", "ENSG00000034693", "ENSG00000018510", # UGCG GBA PEX3 AGPS
          "ENSG00000116171", "ENSG00000164751", "ENSG00000117528", "ENSG00000127980", # SCP2 PEX2 ABCD3 PEX1
          "ENSG00000135218",
          "ENSG00000109819", "ENSG00000117139", "ENSG00000187140" # PPARGC1A KDM5B FOXD3
)

expr_matrix <- melanoma@assays[["RNA"]]@data[pick,]
expr_matrix@Dimnames[[1]] <- c(# "metabolism",
                               "UGCG", "GBA", "PEX3", "AGPS",
                               "SCP2", "PEX2", "ABCD3", "PEX1",
                               "CD36",
                               "PPARGC1A", "KDM5B", "FOXD3"
)

expr_matrix@Dimnames[[2]] <- as.character(Idents(melanoma))
order <- c(which(colnames(expr_matrix) == "T0_CD36neg"),
           which(colnames(expr_matrix) == "T0_CD36pos"),
           which(colnames(expr_matrix) == "T28_CD36neg"),
           which(colnames(expr_matrix) == "T28_CD36pos"))
expr_matrix <- expr_matrix[, order]

# Figure 4D
write.csv(expr_matrix, "D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/expr_matrix_CD36_v2.csv")

# Figure 4C
# obtain metabolic data from melanoma_pop_AUCell.rds and repeat the code above

# -----I am a separation line------------------------------------------------------------------------------------------

# diapause
diapause <- read.csv("D:/McGill/projects/peroxisome/Fan-peroxisome/RAWDATA/diapause_ENSG.csv")

# diapause score after subset
# CD36
melanoma <- readRDS("D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/melanoma_pop_CD36.rds")
melanoma <- subset(melanoma, idents = c("T0_CD36neg", "T0_CD36pos", "T28_CD36neg", "T28_CD36pos"))

# Z-score
z_score <- t(apply(melanoma@assays[["RNA"]]@data, 1, scale))

filter <- match(diapause$Gene, rownames(z_score))
matrix_pick <- z_score[diapause$Gene[-which(is.na(filter))],]
diapause_score <- diapause$Weight[-which(is.na(filter))] %*% matrix_pick
diapause_score <- diapause_score / length(diapause$Weight[-which(is.na(filter))])

rownames(diapause_score) <- "diapause"
melanoma@assays[["RNA"]]@data <- rbind(melanoma@assays[["RNA"]]@data,
                                       diapause_score)
VlnPlot(object = melanoma,
        features = c("diapause"), pt.size = 0)



rownames(diapause_score) <- "diapause"

colnames(diapause_score) <- as.character(Idents(melanoma))
order <- c(which(colnames(diapause_score) == "T0_CD36neg"),
           which(colnames(diapause_score) == "T0_CD36pos"),
           
           which(colnames(diapause_score) == "T28_CD36neg"),
           which(colnames(diapause_score) == "T28_CD36pos"))
expr_matrix <- as.matrix(diapause_score[, order])

# Figure S8H
write.csv(t(expr_matrix), "D:/McGill/projects/peroxisome/Fan-peroxisome/outputs/CD36_diapause_028.csv")
