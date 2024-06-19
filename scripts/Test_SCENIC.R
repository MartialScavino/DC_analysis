setwd("/Users/mscavino/Thèse/ITMO/")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version()
# If your bioconductor version is previous to 4.0, see the section bellow

## Required
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if(!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)


if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")

dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather")

dir.create("data/SCENIC/cisTarget_databases"); setwd("data/SCENIC/cisTarget_databases") # if needed
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}

setwd("/Users/mscavino/Thèse/ITMO/")
