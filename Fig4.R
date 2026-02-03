# Figure 4

# Load packages
rm(list = ls())
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(reshape2)
library(dplyr)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(ensembldb)
library(EnsDb.Hsapiens.v86)

#### FIG. 4A (LR GENE-LEVEL UMAP) ####
# Load data
sr.pbmc <- readRDS("SR.PBMC.S3.rds")
lr.pbmc <- readRDS("LR.PBMC.S3.rds")

# Create a new metadata column to store celltypes
sr.pbmc$celltype.grouped <- sr.pbmc$predicted.celltype.l2
sr.pbmc$celltype.grouped[sr.pbmc$celltype.grouped %in% c("cDC1", "cDC2")] <- "DC"
sr.pbmc$celltype.grouped[sr.pbmc$celltype.grouped %in% c("NK", "NK Proliferating", "NK_CD56bright")] <- "NK"
sr.pbmc$celltype.grouped[sr.pbmc$celltype.grouped %in% c("B memory", "B naive", "B intermediate")] <- "B Cells"
sr.pbmc$celltype.grouped[sr.pbmc$celltype.grouped %in% c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL")] <- "CD4 T Cells"
sr.pbmc$celltype.grouped[sr.pbmc$celltype.grouped %in% c("CD8 Naive", "CD8 TCM", "CD8 TEM")] <- "CD8 T Cells"
sr.pbmc$celltype.grouped[sr.pbmc$celltype.grouped %in% c("Eryth")] <- "Erythroid"
sr.pbmc$celltype.grouped[sr.pbmc$celltype.grouped %in% c("dnT", "ILC")] <- NA

# Add SR gene labels as a metadata column for LR
sr.cell.types <- sr.pbmc@meta.data[rownames(lr.pbmc@meta.data), "celltype.grouped", drop = FALSE]
lr.pbmc@meta.data[["celltype.grouped"]] <- sr.cell.types$celltype.grouped

# Create LR gene UMAP
lr.gene.umap <- DimPlot(lr.pbmc, reduction = "ref.umap", group.by = "celltype.grouped", label = F, repel = T, label.size = 5.5, pt.size = 1.5) + 
  NoLegend() + 
  xlab("UMAP_1") +
  ylab("UMAP_2") + 
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
lr.gene.umap <- LabelClusters(lr.gene.umap, id = "celltype.grouped", fontface = "bold", size = 5.5, repel = TRUE, color = "black")
lr.gene.umap



#### FIG. 4B - 4I (MARKER TRANSCRIPT PLOTS) ####
# Load data
sr.pbmc <- readRDS("SR.PBMC.S3.rds")
lr.pbmc <- readRDS("LR.PBMC.S3.rds")

# Add celltypes to LR object
sr.pbmc <- subset(sr.pbmc, subset = predicted.celltype.l1 %in% c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK"))
lr.pbmc <- subset(lr.pbmc, cells = intersect(colnames(lr.pbmc), colnames(sr.pbmc)))
lr.pbmc$sr.celltype <- sr.pbmc$predicted.celltype.l1[match(colnames(lr.pbmc), colnames(sr.pbmc))]

# Set identities
Idents(lr.pbmc) <- "sr.celltype"

# Find transcript-level markers (positive only)
transcript.markers <- FindAllMarkers(lr.pbmc, assay = "transcript", only.pos = TRUE)

# Filter to highly up-regulated genes
transcript.markers.filtered <- transcript.markers %>%
  group_by(cluster) %>%
  mutate(spec.score = avg_log2FC * (pct.1 - pct.2)) %>%
  dplyr::filter(avg_log2FC > 0.4, pct.1 > 0.25, pct.2 < 0.15)

# Initialize empty DF to store top unique transcripts per cluster
top3.unique <- data.frame()

# Keep track of transcripts already selected
used.transcripts <- c()

# Loop over clusters
for (clust in unique(transcript.markers.filtered$cluster)) {
  
  # Get markers for this cluster, ordered by spec.score
  df <- transcript.markers.filtered %>%
    dplyr::filter(cluster == clust) %>%
    arrange(desc(spec.score))
  
  # Keep only transcripts not already used
  df <- df[!df$gene %in% used.transcripts, ]
  
  # Use top n transcripts
  top_n <- head(df, 3)
  
  # Append to results
  top3.unique <- rbind(top3.unique, top_n)
  
  # Add selected transcripts to used list
  used.transcripts <- c(used.transcripts, top_n$gene)
}

# Use these transcripts for features.vec
features.vec <- top3.unique$gene

# Read in transcript to gene name mapping
t2g <- read.table("t2gnames.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(t2g) <- c("transcript.id", "gene.name")

# Create labeled features
features.labeled <- data.frame(transcript.id.full = features.vec, stringsAsFactors = FALSE) %>%
  left_join(t2g, by = c("transcript.id.full" = "transcript.id")) %>%
  mutate(label = paste0(transcript.id.full, " / ", gene.name),
         cluster = top3.unique$cluster)

View(features.labeled)

# Define gene of interest
gname <- "CD74"

# Create gene annotation plot
obj.cov <- readRDS("sr.pbmc.coverage.rds")
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(obj.cov) <- annotations
annotation.plot <- AnnotationPlot(obj.cov, region = gname, mode = "transcript")

# Define metadata column
split.by <- "predicted.celltype.l1"

# Access gene expression (SR) 
DefaultAssay(sr.pbmc) <- "RNA"
gene.expr <- as.numeric(GetAssayData(sr.pbmc, layer = "data")[gname, ])
meta <- sr.pbmc@meta.data
gene.df <- data.frame(expression = gene.expr)
gene.df[["predicted.celltype.l1"]] = meta[["predicted.celltype.l1"]]

# Plot gene expression
gene.plot <- ggplot(
  gene.df,
  aes_string(x = "predicted.celltype.l1", y = "expression", fill = "predicted.celltype.l1")
) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = "predicted.celltype.l1", y = "Expression", title = gname)

# Access gene expression (LR()
DefaultAssay(lr.pbmc) <- "transcript"

# Load transcript-to-gene mapping
t2g <- read.table("t2gnames.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(t2g) <- c("transcript.id", "gene.name")
tx.mat <- GetAssayData(lr.pbmc, layer = "data")
tx.present <- rownames(tx.mat)
gname.tx <- t2g$transcript.id[t2g$gene.name == gname]
gname.tx <- intersect(gname.tx, tx.present)
tx.data <- tx.mat[gname.tx, , drop = FALSE]
filtered.tx <- gname.tx[Matrix::rowMeans(tx.data) > 0.005]
tx.data <- tx.data[filtered.tx, , drop = FALSE]
meta.lr <- lr.pbmc@meta.data
query.df <- as.data.frame(t(as.matrix(tx.data)))
query.df[["sr.celltype"]] <- meta.lr[["sr.celltype"]]

# Plot transcript expression
plot.df <- reshape2::melt(query.df, id.vars = "sr.celltype", variable.name = "transcript", value.name = "expression")
plot.df[["sr.celltype"]] <- factor(plot.df[["sr.celltype"]], levels = sort(unique(as.character(plot.df[["sr.celltype"]]))))

tx.plot <- ggplot(
  plot.df,
  aes_string(x = "sr.celltype", y = "expression", fill = "sr.celltype")
) +
  geom_boxplot() +
  facet_wrap(~ transcript, ncol = length(filtered.tx)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(x = split.by, y = "Expression")

# Combine plots and visualize
gene.plot / annotation.plot / tx.plot