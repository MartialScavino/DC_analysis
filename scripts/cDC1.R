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


DC <- subset(seu, Population %in% c("cDC1", "cDC2", "pDC", "Mature_DC", "Langerhans_cells"))

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
            reduction = "PCA", theme_use = "theme_blank")

CellDimPlot(srt = DC, group.by = "Phase",
            reduction = "PCA", theme_use = "theme_blank")

DC <- FindNeighbors(DC, dims = 1:15)
DC <- RunUMAP(DC, dims = 1:15, n.neighbors = 10)

CellDimPlot(srt = DC, group.by = "Population",
            reduction = "UMAP", theme_use = "theme_blank")
CellDimPlot(srt = DC, group.by = "Stage",
            reduction = "UMAP", theme_use = "theme_blank", palcolor = Stage_colors)
CellDimPlot(srt = DC, group.by = "Phase",
            reduction = "UMAP", theme_use = "theme_blank")
CellDimPlot(srt = DC, group.by = "orig.ident",
            reduction = "UMAP", theme_use = "theme_blank")

# Langherhans cells
FeatureDimPlot(DC, feature = "Cd207")
FeatureDimPlot(DC, feature = "Siglech")

# cDC2 Tumor deviennent peut être des Langherans
FeatureDimPlot(DC, feature = "Epcam")


CellStatPlot(DC, stat.by = "Stage", group.by = "Population", label = T, palcolor = Stage_colors)
CellStatPlot(DC, stat.by = "Phase", group.by = "Population", label = T)

# Mature DC moins en phase S -> prolifèrent moins ?
# LC quasi uniquement tumorales -> Comme au FACS


#################### cDC1

cDC1 <- subset(DC, Population == "cDC1")

cDC1 <- NormalizeData(cDC1) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
ElbowPlot(cDC1, ndims = 50)

CellDimPlot(srt = cDC1, group.by = "Stage",
            reduction = "PCA", theme_use = "theme_blank", palcolor = Stage_colors)

CellDimPlot(srt = cDC1, group.by = "Phase",
            reduction = "PCA", theme_use = "theme_blank")


cDC1 <- FindNeighbors(cDC1, dims = 1:20)
cDC1 <- RunUMAP(cDC1, dims = 1:20)

CellDimPlot(srt = cDC1, group.by = "Stage",
            reduction = "UMAP", theme_use = "theme_blank", palcolor = Stage_colors)
CellDimPlot(srt = cDC1, group.by = "Phase",
            reduction = "UMAP", theme_use = "theme_blank")

cDC1 <- FindClusters(cDC1, resolution = 0.5)

CellDimPlot(srt = cDC1, group.by = "seurat_clusters",
            reduction = "UMAP", theme_use = "theme_blank")
CellStatPlot(cDC1, stat.by = "Phase", group.by = "seurat_clusters", label = T)
CellStatPlot(cDC1, stat.by = "Stage", group.by = "seurat_clusters", label = T, palcolor = Stage_colors)

### Cluster spé Phase, on va intégrer avec harmony

cDC1 <- RunHarmony(cDC1, "Phase")
cDC1 <- FindNeighbors(cDC1, dims = 1:20, reduction = "harmony")
cDC1 <- RunUMAP(cDC1, dims = 1:20, reduction = "harmony")
cDC1 <- FindClusters(cDC1, resolution = 0.55)

CellDimPlot(srt = cDC1, group.by = "seurat_clusters",
            reduction = "UMAP", theme_use = "theme_blank")

CellDimPlot(srt = cDC1, group.by = "Stage",
            reduction = "UMAP", theme_use = "theme_blank", palcolor = Stage_colors)

CellStatPlot(cDC1, stat.by = "Phase", group.by = "seurat_clusters", label = T)
CellStatPlot(cDC1, stat.by = "Stage", group.by = "seurat_clusters", label = T, palcolor = Stage_colors)


cDC1 <- RunDEtest(cDC1, group_by = 'seurat_clusters', fc.threshold = 1, only.pos = FALSE)
cDC1 <- RunGSEA(cDC1, group_by = "seurat_clusters", species = "Mus_musculus", db = 'GO_BP')

GSEAPlot(cDC1, group_by = "seurat_clusters", group_use = "6", db = "GO_BP")
GSEAPlot(cDC1, group_by = "seurat_clusters", group_use = "6", db = "GO_BP", plot_type = "bar",
direction = "both", topTerm = 20)
GSEAPlot(cDC1, group_by = "seurat_clusters", db = "GO_BP",plot_type = "comparison")

cDC1@tools$DEtest_seurat_clusters$AllMarkers_wilcox

UMAP <- data.frame(cDC1[["umap"]]@cell.embeddings)
cDC1@meta.data["UMAP_1"] <- UMAP$UMAP_1
cDC1@meta.data["UMAP_2"] <- UMAP$UMAP_2

to_export <- cDC1@meta.data[,c("seurat_clusters", "UMAP_1", "UMAP_2")]