library(reticulate)
library(nbHelpers)

remotes::install_github("barkasn/nbHelpers")

setwd("/Users/mscavino/Thèse/ITMO/")
set.seed(123)

seu <- readRDS("data/MERGED_SEU_REANOTATE_V3.RDS")
pd <- import("pandas")
test <- pd$read_pickle("data/merge_velocity.pickle")



spliced <- transpose_dgRMatrix(test$layers["spliced"])
unspliced <- transpose_dgRMatrix(test$layers["unspliced"])
ambiguous <- transpose_dgRMatrix(test$layers["ambiguous"])


colnames(spliced) <- colnames(seu)
colnames(unspliced) <- colnames(seu)
colnames(ambiguous) <- colnames(seu)

rownames(spliced) <- rownames(test$var)
rownames(unspliced) <- rownames(test$var)
rownames(ambiguous) <- rownames(test$var)

seu@assays["spliced"] <- CreateAssayObject(spliced)
seu@assays["unspliced"] <- CreateAssayObject(unspliced)
seu@assays["ambiguous"] <- CreateAssayObject(ambiguous)


library(SCP)
library(Seurat)
library(dplyr)
library(gprofiler2)
library(harmony)
library(ggplot2)
library(SeuratWrappers)
library(velocyto.R)

seu$Stage <- factor(seu$Stage, levels = c("Healthy", "Early_Stages", "Tumor"))


Stage_colors = list(Healthy = "dodgerblue", Early_Stages = "darkorange", Tumor = "tomato2")


DC <- subset(seu, Population %in% c("cDC1", "cDC2", "Mature_DC"))

DC <- NormalizeData(DC)
DC <- FindVariableFeatures(DC)
DC <- ScaleData(DC)
DC <- RunPCA(DC)
ElbowPlot(DC, ndims = 50)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

DC <- CellCycleScoring(DC, s.features = mmus_s, g2m.features = mmus_g2m)

CellDimPlot(srt = DC, group.by = "Population",
            reduction = "PCA", theme_use = "theme_blank")

CellDimPlot(srt = DC, group.by = "Stage",
            reduction = "PCA", theme_use = "theme_blank", palcolor = Stage_colors)

CellDimPlot(srt = DC, group.by = "Phase",
            reduction = "PCA", theme_use = "theme_blank")

DC <- FindNeighbors(DC, dims = 1:15)
DC <- RunUMAP(DC, dims = 1:15, n.neighbors = 5)

CellDimPlot(srt = DC, group.by = "Population",
            reduction = "UMAP", theme_use = "theme_blank")
CellDimPlot(srt = DC, group.by = "Stage",
            reduction = "UMAP", theme_use = "theme_blank", palcolor = Stage_colors)
CellDimPlot(srt = DC, group.by = "Phase",
            reduction = "UMAP", theme_use = "theme_blank")
CellDimPlot(srt = DC, group.by = "orig.ident",
            reduction = "UMAP", theme_use = "theme_blank")



## On perd encore des gènes, à voir comment régler ce probleme
## Peut être virer les gènes non présents

RNA_velo<- DC@assays$RNA[which(rownames(DC@assays$RNA) %in% rownames(DC@assays$spliced)),]
DC@assays["RNA_velo"] <- CreateAssayObject(RNA_velo)

DC <- RunSCVELO(DC, group_by = "Population", linear_reduction = "pca", nonlinear_reduction = "umap", assay_X = "RNA_velo")

VelocityPlot(DC, plot_type = "stream", reduction = "umap", group_by = "Population")
CellDimPlot(DC, group.by = "Population", reduction = "umap", velocity = "stochastic")


# cDC1s, which are under the transcriptional control of IRF8, ID2 and BATF3
#CD103, which is in the migratory cell subset

#DC2 differentiation is driven by the transcription factors IRF4, ID2, ZEB2 and NOTCH2/KLF4 

# pDC differentiation is guided by the transcription factors IRF8, RUNX1, and TCF4 [38]. pDCs 

# CD5+CD163–CD14– cDC2s and CD5–CD163+CD14+ cells, which have been renamed cDC3
# The CD5–CD163+CD14+ cDC2 population was shown to expand in inflammatory diseases, and these cells were renamed as ‘DC3s’.

# Check velocity aussi dans la tumeur uniquement pour regarder si on a des DC qui vont vers un état activé


