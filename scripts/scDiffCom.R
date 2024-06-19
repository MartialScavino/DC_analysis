library(Seurat)
library(anndata)
library(scDiffCom)
library(data.table)
library(dplyr)
library(ggplot2)
library(gdata) # Pour la fonction startswith


setwd("/Users/mscavino/Thèse/ITMO/")

data <- read_h5ad("data/adata_liana.h5ad")
seu <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)

tab <- read.csv('data/adata_epi_metadata.csv', header = T, row.names = 1)
cells = rownames(tab)


seu$Stage <- as.vector(seu$Stage)
seu$Population <- as.vector(seu$Population)


x = rownames(seu@meta.data)
liste_pop <- ifelse(seu@meta.data[x, "Population"] %in% c("Normal_Epith_Cells", "Pre-neoplastic_cells", "Tumor_cells"), "EpithelialCell", seu@meta.data[x, "Population"])
liste_stage <- ifelse(seu@meta.data[x, "Population"] %in% c("Normal_Epith_Cells", "Pre-neoplastic_cells", "Tumor_cells"), tab[x, "Stage"], seu@meta.data[x, "Stage"])


seu[["PopDiff"]] <- liste_pop
seu[["StageDiff"]] <- liste_stage
seu$StageDiff[which(seu$StageDiff == "Early_Stages")] <- "EarlyStages"


### Reduire le nombre de type cellulaire

## Strategie Enlever prolif + Groupes grossiers pour NK et Virer luminal et basal dans epi

names(table(seu$PopDiff))

# NK, Nk1.1, NKT, ILC1, ILC2, Tgd -> Innate lymphocytes
# T CD4, Th17
# Treg
# T CD8
# Epithelial cells
# Endothelial cells 
# B_cells et plasma
# Cdc2 et langherans
# Cdc1
# Mature DC
# Neutro
# Basos


# Retirer prolif et prendre que HES
liste_prolif <- ifelse(startsWith(seu@meta.data[x, "PopDiff"], "Prolif"), TRUE, FALSE)

seu$is_prolif <- liste_prolif


seu.HES <- subset(seu, is_prolif == FALSE  & Stage %in% c("Healthy", "Early_Stages"))

unique(seu.HES$PopDiff)

dico = list(
  CD8 = "T CD8",
  CD4 = "T CD4",
  Treg = "Treg",
  Mo_Mac = "Mo Mac",
  cDC2 = "cDC2",
  Fibroblast = "Fibroblast",
  Tgd = "Innate lymphocytes",
  NK_cells = "Innate lymphocytes",
  cDC1 = "cDC1",
  Basophils = "Basophils",
  B_cells = "B cells",
  Plasma_cells = "B cells",
  Neutrophils = "Neutrophils",
  Endothelial_cells = "Endothelial cells",
  ILC2 = "Innate lymphocytes",
  pDC = "pDC",
  ILC1 = "Innate lymphocytes",
  NK_NK1.1_gd = "Innate lymphocytes",
  Th17 = "T CD4",
  EpithelialCell = "Epithelial cells",
  Luminal_cells = "Epithelial cells",
  Mature_DC = "Mature DC",
  NKT = "Innate lymphocytes",
  Basal_cells = "Epithelial cells",
  Langerhans_cells = "cDC2"
)


annotation <- sapply(seu.HES$PopDiff, function(i){
  return(dico[[i]])
})

unique(annotation)

seu.HES$annotation <- annotation


# Test healthy vs early

#seu.HES <- subset(seu, Stage != "Tumor")
#seu.HES

table(seu.HES[["Stage"]])


seu.HES <- NormalizeData(seu.HES)
seu.HES <- FindVariableFeatures(seu.HES)

seu.HES@assays$RNA@data@x[is.na(seu.HES@assays$RNA@data@x)] <- 0

seu.HES <- ScaleData(seu.HES)
 
# seu.HES <- FindNeighbors(seu.HES, dims = 1:10)
# seu.HES <- RunUMAP(seu.HES, dims = 1:10)
# seu.HES <- FindClusters(seu.HES)
# 
# DimPlot(seu.HES, group.by = "Stage")


# Peut être une nouvelle stratégie : Noter les epi en epithelial cells et 
# changer le Stage des epi pour avoir le vrai stade (en faisant une nouvelle colonne)


# 
# seu.HES$Stage <- as.vector(seu.HES$Stage)
# seu.HES$Population <- as.vector(seu.HES$Population)


# liste_stage = c()
# liste_pop = c()
# for (cell in rownames(seu.HES@meta.data)){ # Si c'est pas une tumorale
#   if (cell %in% cells){ # Si c'est une cellule epitheliale
#     liste_stage = c(liste_stage, tab[cell, "Stage"])
#     liste_pop = c(liste_pop, "EpithelialCell")
#   }
#   else{ # Si c'est pas une cellule epithéliale
#     liste_stage = c(liste_stage, seu.HES@meta.data[cell, "Stage"])
#     liste_pop = c(liste_pop, seu.HES@meta.data[cell, "Population"])
#   }
#   
# }
# 
# seu.HES[["PopDiff"]] <- liste_pop
# seu.HES[["StageDiff"]] <- liste_stage
# seu.HES$StageDiff[which(seu.HES$StageDiff == "Early_Stages")] <- "EarlyStages"


# DimPlot(seu.HES, group.by = "StageDiff")
# DimPlot(seu.HES, group.by = "PopDiff")


###

library(future)
plan(multisession, workers = 4)

scdiffcom_object <- run_interaction_analysis(iterations = 10000,
  seurat_object = seu.HES,
  LRI_species = "mouse",
  seurat_celltype_id = "annotation",
  seurat_condition_id = list(
    column_name = "StageDiff",
    cond1_name = "Healthy",
    cond2_name = "EarlyStages"
  )
)

saveRDS(scdiffcom_object, "data/scdiffcom_H_ES.rds")


CCI_detected <- GetTableCCI(scdiffcom_object, type = "detected", simplified = TRUE)
table(CCI_detected$REGULATION)


ORA_results <- GetTableORA(scdiffcom_object, categories = "all", simplified = TRUE)
names(ORA_results)


ggplot(
  CCI_detected,
  aes(
    x = LOGFC,
    y = -log10(BH_P_VALUE_DE + 1E-2),
    colour = REGULATION
  )
) + geom_point(
) + scale_colour_manual(
  values = c("UP" = "red", "DOWN" = "blue", "FLAT" = "green", "NSC" = "grey")
) + xlab(
  "log(FC)"
) + ylab(
  "-log10(Adj. p-value)"
)



PlotORA(
  object = scdiffcom_object,
category = "LRI",
  regulation = "UP"
) + theme(
  legend.position = c(0.85, 0.4),
  legend.key.size = unit(0.4, "cm")
)


if (!require("visNetwork")) install.packages("visNetwork")
if (!require("igraph")) install.packages("igraph")
if (!require("kableExtra")) install.packages("kableExtra")
if (!require("RColorBrewer")) install.packages("RColorBrewer")


BuildNetwork(
  object = scdiffcom_object,layout_type = "conventional" 
)

reduced_go_terms <- ReduceGO(scdiffcom_object) # optional

BuildShiny(
  scdiffcom_object
)

VlnPlot(seu.HES, "Fpr2", group.by = "Stage")


scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION, CCI_SCORE_Healthy, CCI_SCORE_EarlyStages) %>% 
  filter(P_VALUE_DE < 0.01, RECEIVER_CELLTYPE == "Epithelial cells", REGULATION != "NSC") %>% 
  arrange(-LOGFC) %>% 
  head(20)

names(scdiffcom_object@cci_table_detected)

LRI_HES <- scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(REGULATION != "NSC", P_VALUE_DE < 0.01) %>% 
  arrange(-LOGFC)

LRI_HES %>% head(20)


# Regarder interactions up entre T et ES et entre T et H et faire venn diagramm

seu.HT <- subset(seu, is_prolif == FALSE  & Stage %in% c("Healthy", "Tumor"))


dico = list(
  CD8 = "T CD8",
  CD4 = "T CD4",
  Treg = "Treg",
  Mo_Mac = "Mo Mac",
  cDC2 = "cDC2",
  Fibroblast = "Fibroblast",
  Tgd = "Innate lymphocytes",
  NK_cells = "Innate lymphocytes",
  cDC1 = "cDC1",
  Basophils = "Basophils",
  B_cells = "B cells",
  Plasma_cells = "B cells",
  Neutrophils = "Neutrophils",
  Endothelial_cells = "Endothelial cells",
  ILC2 = "Innate lymphocytes",
  pDC = "pDC",
  ILC1 = "Innate lymphocytes",
  NK_NK1.1_gd = "Innate lymphocytes",
  Th17 = "T CD4",
  EpithelialCell = "Epithelial cells",
  Luminal_cells = "Epithelial cells",
  Mature_DC = "Mature DC",
  NKT = "Innate lymphocytes",
  Basal_cells = "Epithelial cells",
  Langerhans_cells = "cDC2"
)


annotation <- sapply(seu.HT$PopDiff, function(i){
  return(dico[[i]])
})


seu.HT$annotation <- annotation


seu.HT <- NormalizeData(seu.HT)
seu.HT <- FindVariableFeatures(seu.HT)

seu.HT@assays$RNA@data@x[is.na(seu.HT@assays$RNA@data@x)] <- 0

seu.HT <- ScaleData(seu.HT)


seu.HT@meta.data[which(seu.HT$StageDiff == "EarlyStages"), "StageDiff"] <- "Healthy"




scdiffcom_object <- run_interaction_analysis(iterations = 10000,
                                             seurat_object = seu.HT,
                                             LRI_species = "mouse",
                                             seurat_celltype_id = "annotation",
                                             seurat_condition_id = list(
                                               column_name = "StageDiff",
                                               cond1_name = "Healthy",
                                               cond2_name = "Tumor"
                                             )
)

saveRDS(scdiffcom_object, "data/scdiffcom_H_T.rds")


CCI_detected <- GetTableCCI(scdiffcom_object, type = "detected", simplified = TRUE)
table(CCI_detected$REGULATION)


ORA_results <- GetTableORA(scdiffcom_object, categories = "all", simplified = TRUE)
names(ORA_results)


ggplot(
  CCI_detected,
  aes(
    x = LOGFC,
    y = -log10(BH_P_VALUE_DE + 1E-2),
    colour = REGULATION
  )
) + geom_point(
) + scale_colour_manual(
  values = c("UP" = "red", "DOWN" = "blue", "FLAT" = "green", "NSC" = "grey")
) + xlab(
  "log(FC)"
) + ylab(
  "-log10(Adj. p-value)"
)



PlotORA(
  object = scdiffcom_object,
  category = "LRI",
  regulation = "UP"
) + theme(
  legend.position = c(0.85, 0.4),
  legend.key.size = unit(0.4, "cm")
)


if (!require("visNetwork")) install.packages("visNetwork")
if (!require("igraph")) install.packages("igraph")
if (!require("kableExtra")) install.packages("kableExtra")
if (!require("RColorBrewer")) install.packages("RColorBrewer")


BuildNetwork(
  object = scdiffcom_object,layout_type = "conventional" 
)




scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(P_VALUE_DE < 0.05, RECEIVER_CELLTYPE == "Treg", EMITTER_CELLTYPE == "cDC2", REGULATION != "NSC") %>% 
  arrange(-LOGFC)

names(scdiffcom_object@cci_table_detected)

LRI_HES <- scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(REGULATION != "NSC", P_VALUE_DE < 0.01) %>% 
  arrange(-LOGFC)

LRI_HES %>% head(20)





### Early vs Tumor

seu.EST <- subset(seu, is_prolif == FALSE  & Stage %in% c("Early_Stages", "Tumor"))

dico = list(
  CD8 = "T CD8",
  CD4 = "T CD4",
  Treg = "Treg",
  Mo_Mac = "Mo Mac",
  cDC2 = "cDC2",
  Fibroblast = "Fibroblast",
  Tgd = "Innate lymphocytes",
  NK_cells = "Innate lymphocytes",
  cDC1 = "cDC1",
  Basophils = "Basophils",
  B_cells = "B cells",
  Plasma_cells = "B cells",
  Neutrophils = "Neutrophils",
  Endothelial_cells = "Endothelial cells",
  ILC2 = "Innate lymphocytes",
  pDC = "pDC",
  ILC1 = "Innate lymphocytes",
  NK_NK1.1_gd = "Innate lymphocytes",
  Th17 = "T CD4",
  EpithelialCell = "Epithelial cells",
  Luminal_cells = "Epithelial cells",
  Mature_DC = "Mature DC",
  NKT = "Innate lymphocytes",
  Basal_cells = "Epithelial cells",
  Langerhans_cells = "cDC2"
)


annotation <- sapply(seu.EST$PopDiff, function(i){
  return(dico[[i]])
})

seu.EST$annotation <- annotation



seu.EST <- NormalizeData(seu.EST)
seu.EST <- FindVariableFeatures(seu.EST)

seu.EST@assays$RNA@data@x[is.na(seu.EST@assays$RNA@data@x)] <- 0

seu.EST <- ScaleData(seu.EST)

seu.EST@meta.data[which(seu.EST$StageDiff == "Healthy"), "StageDiff"] <- "EarlyStages"

scdiffcom_object <- run_interaction_analysis(iterations = 10000,
                                             seurat_object = seu.EST,
                                             LRI_species = "mouse",
                                             seurat_celltype_id = "annotation",
                                             seurat_condition_id = list(
                                               column_name = "StageDiff",
                                               cond1_name = "EarlyStages",
                                               cond2_name = "Tumor"
                                             )
)

saveRDS(scdiffcom_object, "data/scdiffcom_ES_T.rds")

CCI_detected <- GetTableCCI(scdiffcom_object, type = "detected", simplified = TRUE)
table(CCI_detected$REGULATION)


ORA_results <- GetTableORA(scdiffcom_object, categories = "all", simplified = TRUE)
names(ORA_results)


ggplot(
  CCI_detected,
  aes(
    x = LOGFC,
    y = -log10(BH_P_VALUE_DE + 1E-2),
    colour = REGULATION
  )
) + geom_point(
) + scale_colour_manual(
  values = c("UP" = "red", "DOWN" = "blue", "FLAT" = "green", "NSC" = "grey")
) + xlab(
  "log(FC)"
) + ylab(
  "-log10(Adj. p-value)"
)



PlotORA(
  object = scdiffcom_object,
  category = "LRI",
  regulation = "UP"
) + theme(
  legend.position = c(0.85, 0.4),
  legend.key.size = unit(0.4, "cm")
)


if (!require("visNetwork")) install.packages("visNetwork")
if (!require("igraph")) install.packages("igraph")
if (!require("kableExtra")) install.packages("kableExtra")
if (!require("RColorBrewer")) install.packages("RColorBrewer")


BuildNetwork(
  object = scdiffcom_object,layout_type = "conventional" 
)




scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(P_VALUE_DE < 0.05, RECEIVER_CELLTYPE == "Treg", EMITTER_CELLTYPE == "cDC2", REGULATION != "NSC") %>% 
  arrange(-LOGFC)

names(scdiffcom_object@cci_table_detected)

LRI_HES <- scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(REGULATION != "NSC", P_VALUE_DE < 0.01) %>% 
  arrange(-LOGFC)

LRI_HES %>% head(20)



# Refaire les analyses mais en ne gardant que des types cellulaires candidats


# Celltypes à garder : 
## cDC1
## cDC2
## pDC
## TCD8
## TCD4
## Treg
## Epi
## Endo


keep_cells <- c("cDC1", "cDC2", "pDC", "T CD8", "T CD4", "Treg", "Epithelial cells", "Endothelial cells")


seu.HES <- subset(seu, is_prolif == FALSE  & Stage %in% c("Healthy", "Early_Stages"))
unique(seu.HES$PopDiff)

dico = list(
  CD8 = "T CD8",
  CD4 = "T CD4",
  Treg = "Treg",
  Mo_Mac = "Mo Mac",
  cDC2 = "cDC2",
  Fibroblast = "Fibroblast",
  Tgd = "Innate lymphocytes",
  NK_cells = "Innate lymphocytes",
  cDC1 = "cDC1",
  Basophils = "Basophils",
  B_cells = "B cells",
  Plasma_cells = "B cells",
  Neutrophils = "Neutrophils",
  Endothelial_cells = "Endothelial cells",
  ILC2 = "Innate lymphocytes",
  pDC = "pDC",
  ILC1 = "Innate lymphocytes",
  NK_NK1.1_gd = "Innate lymphocytes",
  Th17 = "T CD4",
  EpithelialCell = "Epithelial cells",
  Luminal_cells = "Epithelial cells",
  Mature_DC = "Mature DC",
  NKT = "Innate lymphocytes",
  Basal_cells = "Epithelial cells",
  Langerhans_cells = "cDC2"
)


annotation <- sapply(seu.HES$PopDiff, function(i){
  return(dico[[i]])
})

unique(annotation)

seu.HES$annotation <- annotation


seu.HES.trim <- subset(seu.HES, annotation %in% keep_cells)

table(seu.HES.trim[["Stage"]])
table(seu.HES.trim[["annotation"]])

seu.HES.trim <- NormalizeData(seu.HES.trim)
seu.HES.trim <- FindVariableFeatures(seu.HES.trim)

seu.HES.trim@assays$RNA@data@x[is.na(seu.HES.trim@assays$RNA@data@x)] <- 0

seu.HES.trim <- ScaleData(seu.HES.trim)

library(future)
plan(multisession, workers = 4)

scdiffcom_object <- run_interaction_analysis(iterations = 10000,
                                             seurat_object = seu.HES.trim,
                                             LRI_species = "mouse",
                                             seurat_celltype_id = "annotation",
                                             seurat_condition_id = list(
                                               column_name = "StageDiff",
                                               cond1_name = "Healthy",
                                               cond2_name = "EarlyStages"
                                             )
)

saveRDS(scdiffcom_object, "data/scdiffcom_H_ES_trim.rds")




seu.HT <- subset(seu, is_prolif == FALSE  & Stage %in% c("Healthy", "Tumor"))

dico = list(
  CD8 = "T CD8",
  CD4 = "T CD4",
  Treg = "Treg",
  Mo_Mac = "Mo Mac",
  cDC2 = "cDC2",
  Fibroblast = "Fibroblast",
  Tgd = "Innate lymphocytes",
  NK_cells = "Innate lymphocytes",
  cDC1 = "cDC1",
  Basophils = "Basophils",
  B_cells = "B cells",
  Plasma_cells = "B cells",
  Neutrophils = "Neutrophils",
  Endothelial_cells = "Endothelial cells",
  ILC2 = "Innate lymphocytes",
  pDC = "pDC",
  ILC1 = "Innate lymphocytes",
  NK_NK1.1_gd = "Innate lymphocytes",
  Th17 = "T CD4",
  EpithelialCell = "Epithelial cells",
  Luminal_cells = "Epithelial cells",
  Mature_DC = "Mature DC",
  NKT = "Innate lymphocytes",
  Basal_cells = "Epithelial cells",
  Langerhans_cells = "cDC2"
)


annotation <- sapply(seu.HT$PopDiff, function(i){
  return(dico[[i]])
})

unique(annotation)

seu.HT$annotation <- annotation


seu.HT.trim <- subset(seu.HT, annotation %in% keep_cells)

table(seu.HT.trim[["Stage"]])
table(seu.HT.trim[["StageDiff"]])

table(seu.HT.trim[["annotation"]])

plan("sequential")

seu.HT.trim <- NormalizeData(seu.HT.trim)
seu.HT.trim <- FindVariableFeatures(seu.HT.trim)

seu.HT.trim@assays$RNA@data@x[is.na(seu.HT.trim@assays$RNA@data@x)] <- 0

seu.HT.trim <- ScaleData(seu.HT.trim)

plan(multisession, workers = 4)

seu.HT.trim$StageDiff[which(seu.HT.trim$StageDiff == "EarlyStages")] <- "Healthy"


scdiffcom_object <- run_interaction_analysis(iterations = 10000,
                                             seurat_object = seu.HT.trim,
                                             LRI_species = "mouse",
                                             seurat_celltype_id = "annotation",
                                             seurat_condition_id = list(
                                               column_name = "StageDiff",
                                               cond1_name = "Healthy",
                                               cond2_name = "Tumor"
                                             )
)

saveRDS(scdiffcom_object, "data/scdiffcom_H_T_trim.rds")




####

seu.EST <- subset(seu, is_prolif == FALSE  & Stage %in% c("Early_Stages", "Tumor"))

dico = list(
  CD8 = "T CD8",
  CD4 = "T CD4",
  Treg = "Treg",
  Mo_Mac = "Mo Mac",
  cDC2 = "cDC2",
  Fibroblast = "Fibroblast",
  Tgd = "Innate lymphocytes",
  NK_cells = "Innate lymphocytes",
  cDC1 = "cDC1",
  Basophils = "Basophils",
  B_cells = "B cells",
  Plasma_cells = "B cells",
  Neutrophils = "Neutrophils",
  Endothelial_cells = "Endothelial cells",
  ILC2 = "Innate lymphocytes",
  pDC = "pDC",
  ILC1 = "Innate lymphocytes",
  NK_NK1.1_gd = "Innate lymphocytes",
  Th17 = "T CD4",
  EpithelialCell = "Epithelial cells",
  Luminal_cells = "Epithelial cells",
  Mature_DC = "Mature DC",
  NKT = "Innate lymphocytes",
  Basal_cells = "Epithelial cells",
  Langerhans_cells = "cDC2"
)


annotation <- sapply(seu.EST$PopDiff, function(i){
  return(dico[[i]])
})

unique(annotation)

seu.EST$annotation <- annotation


seu.EST.trim <- subset(seu.EST, annotation %in% keep_cells)

table(seu.EST.trim[["Stage"]])
table(seu.EST.trim[["StageDiff"]])

table(seu.EST.trim[["annotation"]])

plan("sequential")

seu.EST.trim <- NormalizeData(seu.EST.trim)
seu.EST.trim <- FindVariableFeatures(seu.EST.trim)

seu.EST.trim@assays$RNA@data@x[is.na(seu.EST.trim@assays$RNA@data@x)] <- 0

seu.EST.trim <- ScaleData(seu.EST.trim)

plan(multisession, workers = 4)


seu.EST.trim$StageDiff[which(seu.EST.trim$StageDiff == "Healthy")] <- "EarlyStages"

scdiffcom_object <- run_interaction_analysis(iterations = 10000,
                                             seurat_object = seu.EST.trim,
                                             LRI_species = "mouse",
                                             seurat_celltype_id = "annotation",
                                             seurat_condition_id = list(
                                               column_name = "StageDiff",
                                               cond1_name = "EarlyStages",
                                               cond2_name = "Tumor"
                                             )
)

saveRDS(scdiffcom_object, "data/scdiffcom_ES_T_trim.rds")

