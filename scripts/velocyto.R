library(Seurat)
library(anndata)
library(velocyto.R)
library(SeuratWrappers)



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


# Load loom file
# ldat <- ReadVelocity(file = "outs/velocyto/merged_combine_no_demuxtiplex.loom")


ldat <- readRDS("data/ldata_R.RDS")


rownames(ldat$spliced) <- make.unique(rownames(ldat$spliced))
rownames(ldat$unspliced) <- make.unique(rownames(ldat$unspliced))
rownames(ldat$ambiguous) <- make.unique(rownames(ldat$ambiguous))

# Convert loom to Seurat object
ldat_ser <- as.Seurat(ldat)

# Rename loom object to same convention as integrated object

colnames(ldat_ser)



#new_names <- gsub('Part To Remove', '', colnames(ldat_ser))


# To confirm that we have the same cell names in our intergrated object and our loom object
intersect(colnames(seu), colnames(ldat_ser[, colnames(seu)])) # This works

# Merge objects, creates duplicates of all cells
merge_ser <- merge(seu, ldat_ser)


colnames(seu)[which(endsWith(colnames(seu), "4"))]

test <- sapply(colnames(seu), function(x){
  
  return(substr(x, 1,6))
})


new_names <- gsub('CD45pos', '', colnames(seu))
new_names <- gsub('CD45neg', '', new_names)
new_names <- gsub('45p', '', new_names)
new_names <- gsub('45n', '', new_names)
new_names <- gsub('_', '', new_names)
new_names <- gsub('-', '', new_names)
new_names <- gsub('1', '', new_names)
new_names <- gsub('2', '', new_names)
new_names <- gsub('3', '', new_names)
new_names <- gsub('4', '', new_names)
new_names <- gsub('5', '', new_names)
new_names <- gsub('6', '', new_names)
new_names <- gsub('7', '', new_names)
new_names <- gsub('8', '', new_names)
new_names <- gsub('9', '', new_names)


new_names <- make.unique(new_names)

seu <- RenameCells(seu, new.names = new_names, for.merge = T)


new_names_ldat <- unlist(strsplit(colnames(ldat_ser), ":"))

new_names_ldat <- new_names_ldat[seq(2, length(new_names_ldat), 2)]


new_names_ldat <- make.unique(new_names_ldat)

ldat_ser <- RenameCells(ldat_ser, new.names = new_names_ldat, for.merge = T)



colnames(ldat_ser)
colnames(seu)

merge_ser <- merge(seu, ldat_ser)

# To confirm that we have the same cell names in our intergrated object and our loom object
length(intersect(colnames(seu), colnames(ldat_ser))) # This works


length(intersect(new_names_ldat, new_names))


ldat_ser
seu


length(colnames(seu)) + length(colnames(ldat_ser)) == length(colnames(merge_ser))

merge_ser@meta.data[which(is.na(merge_ser$nFeature_ambiguous) == FALSE & is.na(merge_ser$nFeature_RNA) == FALSE) , ]

