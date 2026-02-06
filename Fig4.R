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

#### FIG. 4A - 4F ####
# Load data
sr.pbmc <- readRDS("SR.PBMC.S3.rds")
lr.pbmc <- readRDS("LR.PBMC.S3.rds")

# Add celltypes to LR object
sr.pbmc <- subset(sr.pbmc, subset = predicted.celltype.l1 %in% c("B", "CD4 T", "CD8 T", "DC", "Mono", "NK"))
lr.pbmc <- subset(lr.pbmc, cells = intersect(colnames(lr.pbmc), colnames(sr.pbmc)))
lr.pbmc$sr.celltype <- sr.pbmc$predicted.celltype.l1[match(colnames(lr.pbmc), colnames(sr.pbmc))]

# Find transcript-level markers (positive only)
Idents(lr.pbmc) <- "sr.celltype"
transcript.markers <- FindAllMarkers(lr.pbmc, assay = "transcript", only.pos = TRUE)

# Filter to highly up-regulated genes
transcript.markers.filtered <- transcript.markers %>%
  group_by(cluster) %>%
  mutate(spec.score = avg_log2FC * (pct.1 - pct.2)) %>%
  dplyr::filter(avg_log2FC > 0.4, pct.1 > 0.25, pct.2 < 0.15)

# Top 3 unique transcripts per cell type
top3.unique <- data.frame()
used.transcripts <- c()

# Loop
for (clust in unique(transcript.markers.filtered$cluster)) {

  df <- transcript.markers.filtered %>%
    dplyr::filter(cluster == clust) %>%
    arrange(desc(spec.score))
  
  df <- df[!df$gene %in% used.transcripts, ]
  top_n <- head(df, 3)
  top3.unique <- rbind(top3.unique, top_n)
  
  used.transcripts <- c(used.transcripts, top_n$gene)
}

# Select features
features.vec <- top3.unique$gene

# Transcript to gene name mapping
t2g <- read.table("t2gnames.txt", header = F, sep = "\t", stringsAsFactors = F)
colnames(t2g) <- c("transcript.id", "gene.name")

# Create labeled features
features.labeled <- data.frame(transcript.id.full = features.vec, stringsAsFactors = FALSE) %>%
  left_join(t2g, by = c("transcript.id.full" = "transcript.id")) %>%
  mutate(label = paste0(transcript.id.full, " / ", gene.name),
         cluster = top3.unique$cluster)
View(features.labeled)

# Gene
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

# Transcript to gene mapping
t2g <- read.table("t2gnames.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(t2g) <- c("transcript.id", "gene.name")
tx.mat <- GetAssayData(lr.pbmc, layer = "data")
tx.present <- rownames(tx.mat)
gname.tx <- t2g$transcript.id[t2g$gene.name == gname]
gname.tx <- intersect(gname.tx, tx.present)
tx.data <- tx.mat[gname.tx, , drop = FALSE]
filtered.tx <- gname.tx[Matrix::rowMeans(tx.data) > 0.005] # Filter
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

# Plot
gene.plot / annotation.plot / tx.plot