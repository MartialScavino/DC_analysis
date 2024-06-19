library(Seurat)
library(anndata)
library(scDiffCom)
library(data.table)
library(dplyr)
library(ggplot2)
library(gdata) # Pour la fonction startswith


setwd("/Users/mscavino/Thèse/ITMO/")


scdiffcom_object <- readRDS("data/scdiffcom_H_ES.rds")



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


BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "bipartite"
)

BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "conventional"
)


BuildShiny(
  scdiffcom_object
)



scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION, CCI_SCORE_Healthy, CCI_SCORE_EarlyStages) %>% 
  filter(P_VALUE_DE < 0.01, RECEIVER_CELLTYPE == "Epithelial cells", REGULATION != "NSC") %>% 
  arrange(-LOGFC) %>% 
  head(20)

names(scdiffcom_object@cci_table_detected)

LRI_HES <- scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(REGULATION %in% c("UP", "DOWN"), P_VALUE_DE < 0.01) %>% 
  arrange(-LOGFC)

LRI_HES %>% head(20)



scdiffcom_object <- readRDS("data/scdiffcom_H_T.rds")



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



BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "bipartite"
)

BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "conventional"
)

BuildShiny(
  scdiffcom_object
)



scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION, CCI_SCORE_Healthy, CCI_SCORE_Tumor) %>% 
  filter(P_VALUE_DE < 0.01, RECEIVER_CELLTYPE == "Epithelial cells", REGULATION != "NSC") %>% 
  arrange(-LOGFC) %>% 
  head(20)

names(scdiffcom_object@cci_table_detected)

LRI_HES <- scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(REGULATION %in% c("UP", "DOWN"), P_VALUE_DE < 0.01) %>% 
  arrange(-LOGFC)

LRI_HES %>% head(20)







scdiffcom_object <- readRDS("data/scdiffcom_ES_T.rds")



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



BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "bipartite"
)

BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "conventional"
)

BuildShiny(
  scdiffcom_object
)



scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION, CCI_SCORE_EarlyStages, CCI_SCORE_Tumor) %>% 
  filter(P_VALUE_DE < 0.01, RECEIVER_CELLTYPE == "Epithelial cells", REGULATION != "NSC") %>% 
  arrange(-LOGFC) %>% 
  head(20)

names(scdiffcom_object@cci_table_detected)

LRI_HES <- scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(REGULATION %in% c("UP", "DOWN"), P_VALUE_DE < 0.01) %>% 
  arrange(-LOGFC)

LRI_HES %>% head(20)



################## TRIMMED

scdiffcom_object <- readRDS("data/scdiffcom_H_ES_trim.rds")


scdiffcom_object



CCI_detected <- GetTableCCI(scdiffcom_object, type = "detected", simplified = TRUE)
table(CCI_detected$REGULATION)
CCI_detected

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



BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "bipartite"
)

BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "conventional"
)

BuildShiny(
  scdiffcom_object
)



scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION, CCI_SCORE_Healthy, CCI_SCORE_EarlyStages) %>% 
  filter(P_VALUE_DE < 0.01, RECEIVER_CELLTYPE == "Epithelial cells", REGULATION != "NSC") %>% 
  arrange(-LOGFC) %>% 
  head(20)

names(scdiffcom_object@cci_table_detected)

LRI_HES <- scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(REGULATION %in% c("UP", "DOWN"), P_VALUE_DE < 0.01) %>% 
  arrange(LOGFC)

LRI_HES %>% head(20)




#### H AND T
scdiffcom_object <- readRDS("data/scdiffcom_H_T_trim.rds")


scdiffcom_object



CCI_detected <- GetTableCCI(scdiffcom_object, type = "detected", simplified = TRUE)
table(CCI_detected$REGULATION)
CCI_detected

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



BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "bipartite"
)

BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "conventional"
)

BuildShiny(
  scdiffcom_object
)



scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION, CCI_SCORE_Healthy, CCI_SCORE_Tumor) %>% 
  filter(P_VALUE_DE < 0.01, RECEIVER_CELLTYPE == "Epithelial cells", REGULATION != "NSC") %>% 
  arrange(-LOGFC) %>% 
  head(20)

names(scdiffcom_object@cci_table_detected)

LRI_HES <- scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(REGULATION %in% c("UP", "DOWN"), P_VALUE_DE < 0.01) %>% 
  arrange(LOGFC)

LRI_HES %>% head(20)




scdiffcom_object <- readRDS("data/scdiffcom_ES_T_trim.rds")

CCI_detected <- GetTableCCI(scdiffcom_object, type = "detected", simplified = TRUE)
table(CCI_detected$REGULATION)
CCI_detected

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



BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "bipartite"
)

BuildNetwork(
  object = scdiffcom_object,
  layout_type =  "conventional"
)

BuildShiny(
  scdiffcom_object
)



scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION, CCI_SCORE_Tumor, CCI_SCORE_EarlyStages) %>% 
  filter(P_VALUE_DE < 0.01, RECEIVER_CELLTYPE == "Epithelial cells", REGULATION != "NSC") %>% 
  arrange(-LOGFC) %>% 
  head(20)

names(scdiffcom_object@cci_table_detected)

LRI_HES <- scdiffcom_object@cci_table_detected %>% 
  select(LRI, EMITTER_CELLTYPE, RECEIVER_CELLTYPE, LOGFC, P_VALUE_DE, REGULATION) %>% 
  filter(REGULATION %in% c("UP", "DOWN"), P_VALUE_DE < 0.01) %>% 
  arrange(-LOGFC)

LRI_HES %>% head(20)







