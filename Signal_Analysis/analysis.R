Sys.setenv(RETICULATE_PYTHON = "/projects/molonc/aparicio_lab/smeans/miniconda3/envs/signal/bin/python")
Sys.setenv(RETICULATE_USE_UV = "FALSE")

library(anndata)
library(reticulate)
library(Matrix) 
library(SIGNAL)
library(irlba)
library(uwot)
library(ggpubr)
library(randomcoloR)
library(cowplot)
library(Seurat)

set.seed(2024)

file_path <- "signal_prepared.h5ad" 

adata <- read_h5ad(file_path)

X <- adata$X
meta <- adata$obs
meta$Donor <- as.factor(meta$sample_id)
meta$Condition <- as.factor(meta$Condition)
meta<- meta[, c("Condition", "Donor")]

X <- t(X)

#Visualization of Raw Data:
Colors = distinctColorPalette(13)
pca_res = irlba(t(X), nv = 50)
raw_emb = as.matrix(pca_res$u %*% diag(pca_res$d))
raw_umap = as.data.frame(umap(raw_emb))
colnames(raw_umap) = c("UMAP1", "UMAP2")
raw_umap = cbind.data.frame(meta, raw_umap)
p1 = ggscatter(raw_umap, x = "UMAP1", y = "UMAP2", size = 0.1, color = "Condition", palette = Colors, legend = "right") + 
  guides(colour = guide_legend(override.aes = list(size = 2)))
p2 = ggscatter(raw_umap, x = "UMAP1", y = "UMAP2", size = 0.1, palette = Colors, legend = "right") + 
  guides(colour = guide_legend(override.aes = list(size = 2))) + facet_wrap(~Donor)

ggsave("raw_umap_condition.png", plot=p1)
ggsave("raw_umap_donor.png", plot=p2)

signal_emb = Run.gcPCA(X, meta, g_factor = "Condition", b_factor = "Donor")


#Results Plots:
signal_umap = as.data.frame(umap(t(signal_emb)))
colnames(signal_umap) = c("UMAP1", "UMAP2")
signal_umap = cbind.data.frame(meta, signal_umap)
q1 = ggscatter(signal_umap, x = "UMAP1", y = "UMAP2", size = 0.1, color = "Condition", palette = Colors, legend = "right") + 
  guides(colour = guide_legend(override.aes = list(size = 2)))
q2 = ggscatter(signal_umap, x = "UMAP1", y = "UMAP2", size = 0.1, palette = Colors, legend = "right") + 
  guides(colour = guide_legend(override.aes = list(size = 2))) + facet_wrap(~Donor)

ggsave("signal_umap_condition.png", plot=q1)
ggsave("signal_umap_donor.png", plot=q2)


b