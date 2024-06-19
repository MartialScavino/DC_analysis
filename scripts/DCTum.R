library(SCP)
library(Seurat)
library(dplyr)
library(gprofiler2)
library(harmony)
library(ggplot2)

setwd("/Users/mscavino/Thèse/ITMO/")
set.seed(123)

seu <- readRDS("data/MERGED_SEU_REANOTATE_V3.RDS")

seu$Stage <- factor(seu$Stage, levels = c("Healthy", "Early_Stages", "Tumor"))


Stage_colors = list(Healthy = "dodgerblue", Early_Stages = "darkorange", Tumor = "tomato2")


DCtum <- subset(seu, Population %in% c("cDC1", "cDC2", "Mature_DC") & Stage == 'Tumor')

DCtum <- NormalizeData(DCtum)
DCtum <- FindVariableFeatures(DCtum)
DCtum <- ScaleData(DCtum)
DCtum <- RunPCA(DCtum)
ElbowPlot(DCtum, ndims = 50)


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

DCtum <- CellCycleScoring(DCtum, s.features = mmus_s, g2m.features = mmus_g2m)

CellDimPlot(srt = DCtum, group.by = "Population",
            reduction = "PCA", theme_use = "theme_blank")

CellDimPlot(srt = DCtum, group.by = "Stage",
            reduction = "PCA", theme_use = "theme_blank")

CellDimPlot(srt = DCtum, group.by = "Phase",
            reduction = "PCA", theme_use = "theme_blank")

DCtum <- RunHarmony(DCtum, "orig.ident")
DCtum <- FindNeighbors(DCtum, dims = 1:15, reduction = "harmony")
DCtum <- RunUMAP(DCtum, dims = 1:15, reduction = "harmony")
DCtum <- FindClusters(DCtum, resolution = 0.35)

CellDimPlot(srt = DCtum, group.by = "Population",
            reduction = "UMAP", theme_use = "theme_blank")
CellDimPlot(srt = DCtum, group.by = "Phase",
            reduction = "UMAP", theme_use = "theme_blank")
CellDimPlot(srt = DCtum, group.by = "orig.ident",
            reduction = "UMAP", theme_use = "theme_blank")
CellDimPlot(srt = DCtum, group.by = "seurat_clusters",
            reduction = "UMAP", theme_use = "theme_blank")
