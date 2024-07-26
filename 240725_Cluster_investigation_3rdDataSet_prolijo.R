library(Seurat)
library(tidyverse)
library(genefilter)
library(ggcorrplot)
library(readxl)
library(GEOquery)
library(Biobase)

wd = "/data_husihuilke/shared/SMB/carcu/Parkinsons_project/"
setwd(wd)

## Define date and filter parameters of seurat objects to import
# date = 221025
# UMIs.p.cell_lowthreshold = 1000 #2000
# UMIs.p.cell_highthreshold = 60000
# percent.rb_threshold = 5
## Import Seurat objects
# SeuObj <- readRDS(paste0(paste0("Outputs/", date, "_SeuInt_fil_mayt", UMIs.p.cell_lowthreshold, "lowt", UMIs.p.cell_highthreshold,"_lowt", percent.rb_threshold, ".rds")))

SeuObj <- readRDS("Outputs/231129_SeuInt_3rd4th_fil_int_mayt1000lowt60000_lowtInfrb_24clusters.rds")

## Reconfigure date for outputs
#date = 221028
date = 240725

# ## Get the average expression of SeuObj cells
# SeuObj@meta.data <- select(SeuObj@meta.data, seurat_clusters, Treat, EGFP, hSNCA,ActiveGenes.p.cell, UMIs.p.cell)
# DefaultAssay(SeuObj) <- "RNA"
# genenames <- rownames(SeuObj)
# clusters_all <- AverageExpression(SeuObj, slot = "data")[["RNA"]] %>% as.data.frame() # average expression of SeuObj cells
# clusters_all$gene_name <- rownames(clusters_all)
# saveRDS(clusters_all, paste0("Outputs/", date, "_Clusters_avg_exp_int_all.rds"))

# ## Get the average xpression of SeuObj cells of Ctrl treatment (we will asses the identity of clusters with the expression of the CTrl treatment)
# clusters_ctrl <- AverageExpression(subset(SeuObj, subset = Treat == "Ctrl"), slot = "data")[["RNA"]] %>% as.data.frame()
# clusters_ctrl$gene_name <- rownames(clusters_ctrl)
# saveRDS(clusters_ctrl, paste0("Outputs/", date, "_Clusters_avg_exp_int_ctrl.rds"))


DefaultAssay(SeuObj) <- "RNA"

## Get the average xpression of SeuObj cells of Ctrl treatment (we will asses the identity of clusters with the expression of the CTrl treatment)
clusters_ctrl <- AverageExpression(subset(SeuObj, subset = Treat == "Ctrl" & orig.ident == "`3rd`"), slot = "data", group.by = "seurat_clusters")[["RNA"]] %>% as.data.frame()
clusters_ctrl$gene_name <- rownames(clusters_ctrl)
saveRDS(clusters_ctrl, paste0("Outputs/", date, "_Clusters_avg_exp_int_ctrl_3rd.rds"))

## Get the average xpression of SeuObj cells of Ctrl treatment (we will asses the identity of clusters with the expression of the CTrl treatment)
clusters_ctrl <- AverageExpression(subset(SeuObj, subset = Treat == "Ctrl" & orig.ident == "`4th`"), slot = "data", group.by = "seurat_clusters")[["RNA"]] %>% as.data.frame()
clusters_ctrl$gene_name <- rownames(clusters_ctrl)
saveRDS(clusters_ctrl, paste0("Outputs/", date, "_Clusters_avg_exp_int_ctrl_4th.rds"))

## Get the average xpression of SeuObj cells of Ctrl treatment (we will asses the identity of clusters with the expression of the CTrl treatment)
clusters_ctrl <- AverageExpression(subset(SeuObj, subset = Treat == "Ctrl"), slot = "data", group.by = "seurat_clusters")[["RNA"]] %>% as.data.frame()
clusters_ctrl$gene_name <- rownames(clusters_ctrl)
saveRDS(clusters_ctrl, paste0("Outputs/", date, "_Clusters_avg_exp_int_ctrl_all.rds"))
