# Supp Figure 4

# Load packages
library(Seurat)
library(Matrix)
library(R.utils)
library(Matrix.utils)
library(ggplot2)
library(patchwork)
library(Rsamtools)
library(reshape2)
library(dplyr)
library(tidyr)

#### SUPP FIG. 4A ####
# Load data
sr.pbmc <- readRDS("SR.PBMC.S3.rds")
lr.pbmc <- readRDS("LR.PBMC.S3.rds")

# Create a new metadata column to store filtered celltypes
sr.pbmc$celltype.filtered <- dplyr::case_when(
  
  # Monocytes
  sr.pbmc$predicted.celltype.l2 %in% c("CD14 Mono") ~ "CD14 Mono",
  sr.pbmc$predicted.celltype.l2 %in% c("CD16 Mono") ~ "CD16 Mono",
  
  # DC (combine ASDC, cDC1, cDC2)
  sr.pbmc$predicted.celltype.l2 %in% c("ASDC", "cDC1", "cDC2") ~ "DC",
  
  # Plasmablast
  sr.pbmc$predicted.celltype.l2 %in% c("Plasmablast") ~ "Plasmablast",
  
  # B cells (combine naive, memory, intermediate)
  sr.pbmc$predicted.celltype.l2 %in% c("B naive", "B memory", "B intermediate") ~ "B Cells",
  
  # HSPC
  sr.pbmc$predicted.celltype.l2 %in% c("HSPC") ~ "HSPC",
  
  # Erythroid
  sr.pbmc$predicted.celltype.l2 %in% c("Eryth") ~ "Erythroid",
  
  # Platelet
  sr.pbmc$predicted.celltype.l2 %in% c("Platelet") ~ "Platelet",
  
  # pDC
  sr.pbmc$predicted.celltype.l2 %in% c("pDC") ~ "pDC",
  
  # CD4 T cells (combine naive, TCM, TEM, CTL)
  sr.pbmc$predicted.celltype.l2 %in% c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL") ~ "CD4 T Cells",
  
  # Treg
  sr.pbmc$predicted.celltype.l2 %in% c("Treg") ~ "Treg",
  
  # CD8 T cells (combine naive, TCM, TEM)
  sr.pbmc$predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TCM", "CD8 TEM") ~ "CD8 T Cells",
  
  # MAIT
  sr.pbmc$predicted.celltype.l2 %in% c("MAIT") ~ "MAIT",
  
  # CD4 proliferating
  sr.pbmc$predicted.celltype.l2 %in% c("CD4 Proliferating") ~ "CD4 Proliferating",
  
  # NK
  sr.pbmc$predicted.celltype.l2 %in% c("NK", "NK Proliferating", "NK_CD56bright") ~ "NK",
  
  # Everything else → ignored (NA)
  TRUE ~ NA_character_
)

# Add filtered celltypes to LR object
sr.cell.types <- sr.pbmc@meta.data[rownames(lr.pbmc@meta.data), "celltype.filtered", drop = FALSE]
lr.pbmc@meta.data[["celltype.filtered"]] <- sr.cell.types$celltype.filtered
umap.cells <- WhichCells(lr.pbmc, expression = !is.na(celltype.filtered))

# LR transcript UMAP
lr.transcript.umap <- DimPlot(lr.pbmc, 
                              reduction = "transcript.umap", 
                              group.by = "celltype.filtered", 
                              cells = umap.cells,
                              label = F, 
                              repel = T, 
                              label.size = 5.5, 
                              pt.size = 1.0) + 
  NoLegend() + 
  xlab("UMAP_1") +
  ylab("UMAP_2") + 
  theme(plot.title = element_blank())
lr.transcript.umap <- LabelClusters(lr.transcript.umap, id = "celltype.filtered", fontface = "bold", size = 5.5, repel = TRUE, color = "black")
lr.transcript.umap



#### SUPP FIG. 4B - 4G ####
# Load data
sr.pbmc <- readRDS("SR.PBMC.S3.rds")
lr.pbmc <- readRDS("LR.PBMC.S3.rds")

# Create a new metadata column to store filtered celltypes
sr.pbmc$celltype.filtered <- dplyr::case_when(
  
  # Monocytes
  sr.pbmc$predicted.celltype.l2 %in% c("CD14 Mono") ~ "CD14 Mono",
  sr.pbmc$predicted.celltype.l2 %in% c("CD16 Mono") ~ "CD16 Mono",
  
  # DC (combine ASDC, cDC1, cDC2)
  sr.pbmc$predicted.celltype.l2 %in% c("ASDC", "cDC1", "cDC2") ~ "DC",
  
  # Plasmablast
  sr.pbmc$predicted.celltype.l2 %in% c("Plasmablast") ~ "Plasmablast",
  
  # B cells (combine naive, memory, intermediate)
  sr.pbmc$predicted.celltype.l2 %in% c("B naive", "B memory", "B intermediate") ~ "B Cells",
  
  # HSPC
  sr.pbmc$predicted.celltype.l2 %in% c("HSPC") ~ "HSPC",
  
  # Erythroid
  sr.pbmc$predicted.celltype.l2 %in% c("Eryth") ~ "Erythroid",
  
  # Platelet
  sr.pbmc$predicted.celltype.l2 %in% c("Platelet") ~ "Platelet",
  
  # pDC
  sr.pbmc$predicted.celltype.l2 %in% c("pDC") ~ "pDC",
  
  # CD4 T cells (combine naive, TCM, TEM, CTL)
  sr.pbmc$predicted.celltype.l2 %in% c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL") ~ "CD4 T Cells",
  
  # Treg
  sr.pbmc$predicted.celltype.l2 %in% c("Treg") ~ "Treg",
  
  # CD8 T cells (combine naive, TCM, TEM)
  sr.pbmc$predicted.celltype.l2 %in% c("CD8 Naive", "CD8 TCM", "CD8 TEM") ~ "CD8 T Cells",
  
  # MAIT
  sr.pbmc$predicted.celltype.l2 %in% c("MAIT") ~ "MAIT",
  
  # CD4 proliferating
  sr.pbmc$predicted.celltype.l2 %in% c("CD4 Proliferating") ~ "CD4 Proliferating",
  
  # NK
  sr.pbmc$predicted.celltype.l2 %in% c("NK", "NK Proliferating", "NK_CD56bright") ~ "NK",
  
  # Everything else → ignored (NA)
  TRUE ~ NA_character_
)

# Add filtered celltypes to LR object
sr.cell.types <- sr.pbmc@meta.data[rownames(lr.pbmc@meta.data), "celltype.filtered", drop = FALSE]
lr.pbmc@meta.data[["celltype.filtered"]] <- sr.cell.types$celltype.filtered
umap.cells <- WhichCells(lr.pbmc, expression = !is.na(celltype.filtered))
DefaultAssay(lr.pbmc) <- "transcript"
Idents(lr.pbmc) <- "celltype.filtered"

### IL7R
# Define celltype and transcript
celltype <- "CD4 T Cells"
gname <- "IL7R"
tname <- "ENST00000303115"

# Create feature plot
tx.umap <- FeaturePlot(lr.pbmc, 
                       reduction = "transcript.umap",
                       feature = tname, 
                       cells = umap.cells,
                       repel = T, 
                       label = F, 
                       pt.size = 0.75) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(paste(tname, "/", gname)) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
tx.umap <- LabelClusters(tx.umap, 
                         id = "ident", 
                         clusters = celltype,
                         size = 7,
                         repel = TRUE,
                         fontface = "bold", 
                         color = "black")

# Plot 
tx.umap

### CD8A
# Define celltype and transcript
celltype <- "CD8 T Cells"
gname <- "CD8A"
tname <- "ENST00000352580"

# Create feature plot
tx.umap <- FeaturePlot(lr.pbmc, 
                       reduction = "transcript.umap",
                       feature = tname, 
                       repel = T, 
                       label = F, 
                       pt.size = 0.75) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(paste(tname, "/", gname)) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
tx.umap <- LabelClusters(tx.umap, 
                         id = "ident", 
                         clusters = celltype,
                         size = 7,
                         repel = TRUE,
                         fontface = "bold", 
                         color = "black")

# Plot 
tx.umap

### NKG7
# Define celltype and transcript
celltype <- "NK"
gname <- "NKG7"
tname <- "ENST00000221978"

# Create feature plot
tx.umap <- FeaturePlot(lr.pbmc, 
                       reduction = "transcript.umap",
                       feature = tname, 
                       repel = T, 
                       label = F, 
                       pt.size = 0.75) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(paste(tname, "/", gname)) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
tx.umap <- LabelClusters(tx.umap, 
                         id = "ident", 
                         clusters = celltype,
                         size = 7,
                         repel = TRUE,
                         fontface = "bold", 
                         color = "black")

# Plot 
tx.umap

### SERPINA1
# Define celltype and transcript
celltype <- c("CD14 Mono", "CD16 Mono")
gname <- "SERPINA1"
tname <- "ENST00000636712"

# Create feature plot
tx.umap <- FeaturePlot(lr.pbmc, 
                       reduction = "transcript.umap",
                       feature = tname, 
                       repel = T, 
                       label = F, 
                       pt.size = 0.75) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(paste(tname, "/", gname)) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
tx.umap <- LabelClusters(tx.umap, 
                         id = "ident", 
                         clusters = celltype,
                         size = 7,
                         repel = TRUE,
                         fontface = "bold", 
                         color = "black")

# Plot 
tx.umap

### CST3
# Define celltype and transcript
celltype <- c("DC", "pDC")
gname <- "CST3"
tname <- "ENST00000376925"

# Create feature plot
tx.umap <- FeaturePlot(lr.pbmc, 
                       reduction = "transcript.umap",
                       feature = tname, 
                       repel = T, 
                       label = F, 
                       pt.size = 0.75) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(paste(tname, "/", gname)) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
tx.umap <- LabelClusters(tx.umap, 
                         id = "ident", 
                         clusters = celltype,
                         size = 7,
                         repel = TRUE,
                         fontface = "bold", 
                         color = "black")

# Plot 
tx.umap

### CD79A
# Define celltype and transcript
celltype <- "B Cells"
gname <- "CD79A"
tname <- "ENST00000221972"

# Create feature plot
tx.umap <- FeaturePlot(lr.pbmc, 
                       reduction = "transcript.umap",
                       feature = tname, 
                       repel = T, 
                       label = F, 
                       pt.size = 0.75) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle(paste(tname, "/", gname)) +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
tx.umap <- LabelClusters(tx.umap, 
                         id = "ident", 
                         clusters = celltype,
                         size = 7,
                         repel = TRUE,
                         fontface = "bold", 
                         color = "black")

# Plot 
tx.umap