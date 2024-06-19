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

VlnPlot(DC, "Havcr2", group.by = "Stage", split.by = "Population")
VlnPlot(DC, "Havcr2", group.by = "Population")

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
CellStatPlot(cDC1, stat.by = "Phase", group.by = "seurat_clusters", label = T)
CellStatPlot(cDC1, stat.by = "Stage", group.by = "seurat_clusters", label = T, palcolor = Stage_colors)



t <- FindMarkers(cDC1, "Early_Stages", "Healthy", group.by = "Stage", subset.ident = 6)

test_top20 <- rownames(t%>% 
                    arrange(-avg_log2FC) %>% 
                    head(20))

cDC1_6_1 <-  FindMarkers(cDC1, 6, 1, only.pos = T)
cDC1_6_1 <- cDC1_6_1 %>% filter(p_val_adj < 0.01)

cDC1_6_1

top50 <- rownames(cDC1_6_1 %>% 
  filter(p_val_adj < 0.01) %>% 
  arrange(-avg_log2FC) %>% 
  head(50))


FeatureDimPlot(cDC1, top50[1:5], palette = "RdYlBu")

ht <- GroupHeatmap(
  srt = cDC1,
  top50,
  group.by = c("seurat_clusters", "Stage"),
  heatmap_palette = "RdYlBu",
  cell_annotation = c("Phase", "Stage"),
  cell_annotation_palette = c("Dark2", "rickandmorty"),
  group_palcolor = list(seurat_clusters = NULL, Stage = Stage_colors),
  show_column_names = TRUE, add_dot = TRUE, add_reticle = TRUE, flip = T, column_names_rot = 45,
)
ht$plot


FeatureStatPlot(DC, c("Ccr7"), group.by = "Population", plot_type = "box", add_trend = T)

FeaturePlot(DC, "Ccr7")



FeatureStatPlot(cDC1, c("S100a10", "Abracl", "Tmsb10"), group.by = "Stage", plot_type = "box", add_trend = T)
FeatureStatPlot(cDC1, top50[1:10], group.by = "Stage", plot_type = "box")
FeatureStatPlot(cDC1, c("S100a10", "Abracl", "Tmsb10"), group.by = "seurat_clusters", plot_type = "box")


m <- FindMarkers(cDC1, ident.1 = "Early_Stages", group.by = "Stage", only.pos = T)

top50es <- rownames(m %>% 
                    filter(p_val_adj < 0.01) %>% 
                    arrange(-avg_log2FC) %>% 
                    head(50))


ht <- GroupHeatmap(
  srt = cDC1,
  top50es,
  group.by = c("seurat_clusters", "Stage"),
  heatmap_palette = "RdYlBu",
  cell_annotation = c("Phase", "Stage"),
  cell_annotation_palette = c("Dark2", "rickandmorty"),
  group_palcolor = list(seurat_clusters = NULL, Stage = Stage_colors),
  show_column_names = TRUE, add_dot = TRUE, add_reticle = TRUE, flip = T, column_names_rot = 45,
)
ht$plot

m %>% 
  filter(p_val_adj < 0.01) %>% 
  arrange(-avg_log2FC) %>% 
  head(50)

FeatureStatPlot(cDC1, c("S100a10", "Srsf5", "Scand1", "S100a6", "Ubac2", "Fcer1g"), group.by = "Stage", plot_type = "box", add_trend = T)



m <- FindMarkers(cDC1, ident.1 = "Early_Stages", group.by = "Stage", ident.2 = "Healthy", only.pos = T)

top50es_t <- rownames(m %>% 
                      filter(p_val_adj < 0.01) %>% 
                      arrange(-avg_log2FC) %>% 
                      head(50))

ht <- GroupHeatmap(
  srt = cDC1,
  top50es_t,
  group.by = c("seurat_clusters", "Stage"),
  heatmap_palette = "RdYlBu",
  cell_annotation = c("Phase", "Stage"),
  cell_annotation_palette = c("Dark2", "rickandmorty"),
  group_palcolor = list(seurat_clusters = NULL, Stage = Stage_colors),
  show_column_names = TRUE, add_dot = TRUE, add_reticle = TRUE, flip = T, column_names_rot = 45,
)
ht$plot


library(enrichR)
library(enrichplot)
library(plotly)

db <- enrichR::listEnrichrDbs()
db$libraryName

test <- enrichr(rownames(cDC1_6_1), databases = "Reactome_2022")
temp <- unlist(strsplit(test$Reactome_2022$Overlap, "/"))
gene_ratio <- as.integer(temp[seq(1,length(temp),2)]) / as.integer(temp[seq(2,length(temp),2)])
test$Reactome_2022["gene_ratio"] <- gene_ratio
new_df <- test$Reactome_2022 %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  arrange(-gene_ratio)

new_df$Term <- factor(new_df$Term, levels = rev(new_df$Term))

ggplotly(ggplot(new_df[1:30,], aes(x = gene_ratio, y = Term, color = Adjusted.P.value)) +
  geom_point(aes(text = Genes)) + 
  theme_light())
# Ya des trucs interessants -> Go faire un script spé DC1
## MsigDB aussi

database = "GO_Biological_Process_2023"
test <- enrichr(rownames(cDC1_6_1), databases = database)
temp <- unlist(strsplit(test[[database]]$Overlap, "/"))
gene_ratio <- as.integer(temp[seq(1,length(temp),2)]) / as.integer(temp[seq(2,length(temp),2)])
test[[database]]["gene_ratio"] <- gene_ratio
new_df <- test[[database]] %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  arrange(-gene_ratio)

new_df$Term <- factor(new_df$Term, levels = rev(new_df$Term))

ggplotly(ggplot(new_df[1:30,], aes(x = gene_ratio, y = Term, color = Adjusted.P.value)) +
           geom_point(aes(text = Genes)) + 
           theme_light())

db$libraryName
database = "HDSigDB_Mouse_2021"
test <- enrichr(rownames(cDC1_6_1), databases = database)
temp <- unlist(strsplit(test[[database]]$Overlap, "/"))
gene_ratio <- as.integer(temp[seq(1,length(temp),2)]) / as.integer(temp[seq(2,length(temp),2)])
test[[database]]["gene_ratio"] <- gene_ratio
new_df <- test[[database]] %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  arrange(-gene_ratio)

new_df$Term <- factor(new_df$Term, levels = rev(new_df$Term))

ggplotly(ggplot(new_df[1:30,], aes(x = gene_ratio, y = Term, color = Adjusted.P.value)) +
           geom_point(aes(text = Genes)) + 
           theme_light())












FeaturePlot(DC, "Xcr1")

markers_cDC1_clusters <- FindAllMarkers(cDC1, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = 0.1)

top5markers_cDC1_clusters <- markers_cDC1_clusters %>% 
  filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 5)

top5markers_cDC1_clusters




ht <- GroupHeatmap(
  srt = cDC1,
  top5markers_cDC1_clusters$gene,
  group.by = c("seurat_clusters", "Stage"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase"),
  cell_annotation_palette = c("Dark2"),
  group_palcolor = list(seurat_clusters = NULL, Stage = Stage_colors),
  show_column_names = TRUE, add_dot = TRUE, add_reticle = TRUE, flip = T, column_names_rot = 45,
)
ht$plot


## Up cDC1 Tumor vs others

up_cDC1_Stage <- FindAllMarkers(cDC1, only.pos = T, idents = "Stage")


cDC1 <- RunDEtest(cDC1, group_by = "Stage", fc.threshold = 1.5, only.pos = T)
cDC1
cDC1 <- RunEnrichment(cDC1, group_by = "Stage", 
              DE_threshold = "avg_log2FC > log2(0.58) & p_val_adj < 0.05",
              species = "Mus_musculus", db = c("GO_BP", "Reactome", "KEGG"))

cDC1@tools$DEtest_custom$AllMarkers_wilcox

EnrichmentPlot(
  srt = cDC1, group_by = "Stage", plot_type = "comparison"
)

EnrichmentPlot(
  srt = cDC1, group_by = "Stage", group_use = "Tumor",
  plot_type = "network"
)
EnrichmentPlot(
  srt = cDC1, group_by = "Stage", group_use = "Tumor",
  plot_type = "enrichmap"
)


EnrichmentPlot(
  srt = cDC1, group_by = "Stage", group_use = "Early_Stages",
  plot_type = "network"
)
EnrichmentPlot(
  srt = cDC1, group_by = "Stage", group_use = "Early_Stages",
  plot_type = "enrichmap"
)


EnrichmentPlot(
  srt = cDC1, group_by = "Stage", group_use = "Healthy",
  plot_type = "network"
)
EnrichmentPlot(
  srt = cDC1, group_by = "Stage", group_use = "Healthy",
  plot_type = "enrichmap"
)



cDC1 <- RunDEtest(cDC1, group_by = "seurat_clusters", fc.threshold = 1.5, only.pos = T)
cDC1 <- RunEnrichment(cDC1, group_by = "seurat_clusters",
                      DE_threshold = "avg_log2FC > log2(0.58) & p_val_adj < 0.05",
                      species = "Mus_musculus", db = c("GO_BP", "Reactome", "KEGG"))

EnrichmentPlot(
  srt = cDC1, group_by = "seurat_clusters",
  plot_type = "comparison"
)


EnrichmentPlot(
  srt = cDC1, group_by = "seurat_clusters",group_use = 6,
  plot_type = "enrichmap"
)

EnrichmentPlot(
  srt = cDC1, group_by = "seurat_clusters",group_use = 6,
  plot_type = "network"
)



################# cDC2


cDC2 <- subset(DC, Population == "cDC2")

cDC2 <- NormalizeData(cDC2) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
ElbowPlot(cDC2, ndims = 50)

CellDimPlot(srt = cDC2, group.by = "Stage",
            reduction = "PCA", theme_use = "theme_blank", palcolor = Stage_colors)

CellDimPlot(srt = cDC2, group.by = "Phase",
            reduction = "PCA", theme_use = "theme_blank")


cDC2 <- FindNeighbors(cDC2, dims = 1:20)
cDC2 <- RunUMAP(cDC2, dims = 1:20)

CellDimPlot(srt = cDC2, group.by = "Stage",
            reduction = "UMAP", theme_use = "theme_blank", palcolor = Stage_colors)
CellDimPlot(srt = cDC2, group.by = "Phase",
            reduction = "UMAP", theme_use = "theme_blank")

cDC2 <- FindClusters(cDC2, resolution = 0.6)

CellDimPlot(srt = cDC2, group.by = "seurat_clusters",
            reduction = "UMAP", theme_use = "theme_blank", label = TRUE, label_insitu = TRUE, label.size = 5)
CellStatPlot(cDC2, stat.by = "Phase", group.by = "seurat_clusters", label = T)
CellStatPlot(cDC2, stat.by = "Stage", group.by = "seurat_clusters", label = T, palcolor = Stage_colors)

markers_cDC2_clusters <- FindAllMarkers(cDC2, only.pos = T, min.pct = 0.1, logfc.threshold = 0.25, min.diff.pct = 0.1)

top5markers_cDC2_clusters <- markers_cDC2_clusters %>% 
  filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 5)


ht <- GroupHeatmap(
  srt = cDC2,
  top5markers_cDC2_clusters$gene,
  group.by = c("seurat_clusters", "Stage"),
  heatmap_palette = "YlOrRd",
  cell_annotation = c("Phase"),
  cell_annotation_palette = c("Dark2"),
  group_palcolor = list(seurat_clusters = NULL, Stage = Stage_colors),
  show_column_names = TRUE, add_dot = TRUE, add_reticle = TRUE, flip = T, column_names_rot = 45,
)
ht$plot


cDC2_2_6 <- FindMarkers(cDC2, 2, 6, only.pos = T)

top20 <- rownames(cDC2_2_6 %>% 
  filter(p_val_adj < 0.01) %>% 
  arrange(-avg_log2FC) %>% 
  head(20))


ht <- GroupHeatmap(
  srt = cDC2,
  top20,
  group.by = c("seurat_clusters", "Stage"),
  heatmap_palette = "RdYlBu",
  cell_annotation = c("Phase", "Stage"),
  cell_annotation_palette = c("Dark2", "rickandmorty"),
  group_palcolor = list(seurat_clusters = NULL, Stage = Stage_colors),
  show_column_names = TRUE, add_dot = TRUE, add_reticle = TRUE, flip = T, column_names_rot = 45,
)
ht$plot


cDC2_5_6 <- FindMarkers(cDC2, 5, 6, only.pos = T)

top20 <- rownames(cDC2_5_6 %>% 
                    filter(p_val_adj < 0.01) %>% 
                    arrange(-avg_log2FC) %>% 
                    head(20))


ht <- GroupHeatmap(
  srt = cDC2,
  top20,
  group.by = c("seurat_clusters", "Stage"),
  heatmap_palette = "RdYlBu",
  cell_annotation = c("Phase", "Stage"),
  cell_annotation_palette = c("Dark2", "rickandmorty"),
  group_palcolor = list(seurat_clusters = NULL, Stage = Stage_colors),
  show_column_names = TRUE, add_dot = TRUE, add_reticle = TRUE, flip = T, column_names_rot = 45,
)
ht$plot


################ pDC

pDC <- subset(DC, Population == "pDC")

pDC <- NormalizeData(pDC) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
ElbowPlot(pDC, ndims = 50)

CellDimPlot(srt = pDC, group.by = "Stage",
            reduction = "PCA", theme_use = "theme_blank", palcolor = Stage_colors)

CellDimPlot(srt = pDC, group.by = "Phase",
            reduction = "PCA", theme_use = "theme_blank")


pDC <- FindNeighbors(pDC, dims = 1:20)
pDC <- RunUMAP(pDC, dims = 1:20)

CellDimPlot(srt = pDC, group.by = "Stage",
            reduction = "UMAP", theme_use = "theme_blank", palcolor = Stage_colors)
CellDimPlot(srt = pDC, group.by = "Phase",
            reduction = "UMAP", theme_use = "theme_blank")

pDC <- FindClusters(pDC, resolution = 0.6)

CellDimPlot(srt = pDC, group.by = "seurat_clusters",
            reduction = "UMAP", theme_use = "theme_blank", label = TRUE, label_insitu = TRUE, label.size = 5)
CellStatPlot(pDC, stat.by = "Phase", group.by = "seurat_clusters", label = T)
CellStatPlot(pDC, stat.by = "Stage", group.by = "seurat_clusters", label = T, palcolor = Stage_colors)


pDC_7_5 <- FindMarkers(pDC, 7, 5, only.pos = T)

top20 <- rownames(pDC_7_5 %>% 
                    filter(p_val_adj < 0.01) %>% 
                    arrange(-avg_log2FC) %>% 
                    head(20))


ht <- GroupHeatmap(
  srt = pDC,
  top20,
  group.by = c("seurat_clusters", "Stage"),
  heatmap_palette = "RdYlBu",
  cell_annotation = c("Phase", "Stage"),
  cell_annotation_palette = c("Dark2", "rickandmorty"),
  group_palcolor = list(seurat_clusters = NULL, Stage = Stage_colors),
  show_column_names = TRUE, add_dot = TRUE, add_reticle = TRUE, flip = T, column_names_rot = 45,
)
ht$plot


pDC_0_5 <- FindMarkers(pDC, 0, 5, only.pos = T)

top20 <- rownames(pDC_0_5 %>% 
                    filter(p_val_adj < 0.01) %>% 
                    arrange(-avg_log2FC) %>% 
                    head(20))


ht <- GroupHeatmap(
  srt = pDC,
  top20,
  group.by = c("seurat_clusters", "Stage"),
  heatmap_palette = "RdYlBu",
  cell_annotation = c("Phase", "Stage"),
  cell_annotation_palette = c("Dark2", "rickandmorty"),
  group_palcolor = list(seurat_clusters = NULL, Stage = Stage_colors),
  show_column_names = TRUE, add_dot = TRUE, add_reticle = TRUE, flip = T, column_names_rot = 45,
)
ht$plot



########### Mature DC

Mature <- subset(DC, Population == "Mature_DC")

Mature <- NormalizeData(Mature) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()
ElbowPlot(Mature, ndims = 50)

CellDimPlot(srt = Mature, group.by = "Stage",
            reduction = "PCA", theme_use = "theme_blank", palcolor = Stage_colors)

CellDimPlot(srt = Mature, group.by = "Phase",
            reduction = "PCA", theme_use = "theme_blank")


Mature <- FindNeighbors(Mature, dims = 1:20)
Mature <- RunUMAP(Mature, dims = 1:20)

CellDimPlot(srt = Mature, group.by = "Stage",
            reduction = "UMAP", theme_use = "theme_blank", palcolor = Stage_colors, pt.size = 5)
CellDimPlot(srt = Mature, group.by = "Phase",
            reduction = "UMAP", theme_use = "theme_blank")

Mature <- FindClusters(Mature, resolution = 0.6)

CellDimPlot(srt = Mature, group.by = "seurat_clusters",
            reduction = "UMAP", theme_use = "theme_blank", label = TRUE, label_insitu = TRUE, label.size = 5)
CellStatPlot(Mature, stat.by = "Phase", group.by = "seurat_clusters", label = T)
CellStatPlot(Mature, stat.by = "Stage", group.by = "seurat_clusters", label = T, palcolor = Stage_colors)


Mature_4_5 <- FindMarkers(Mature, 4, 5, only.pos = T)

top20 <- rownames(Mature_4_5 %>% 
                    filter(p_val_adj < 0.01) %>% 
                    arrange(-avg_log2FC) %>% 
                    head(20))


ht <- GroupHeatmap(
  srt = Mature,
  top20,
  group.by = c("seurat_clusters", "Stage"),
  heatmap_palette = "RdYlBu",
  cell_annotation = c("Phase", "Stage"),
  cell_annotation_palette = c("Dark2", "rickandmorty"),
  group_palcolor = list(seurat_clusters = NULL, Stage = Stage_colors),
  show_column_names = TRUE, add_dot = TRUE, add_reticle = TRUE, flip = T, column_names_rot = 45,
)
ht$plot

