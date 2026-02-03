# Supplementary Figure 2

# Load packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(viridis)
library(dplyr)

#### SUPPLEMENTARY FIG. 2B (SINGLE-CELL METRICS) ####
# Load data
lr.pbmc <- readRDS("LR.PBMC.S3.rds")

# Format metadata
lr.pbmc$sample <- "BenchDrop-seq"
tx.meta <- lr.pbmc@meta.data %>%
  dplyr::select(
    nCount_RNA,
    nFeature_RNA,
    nFeature_transcript,
    sample
  )

# Visualize number of UMIs per cell
umis.per.cell <- ggplot(
  tx.meta,
  aes(x = sample, y = nCount_RNA, fill = sample)
) +
  geom_violin(scale = "width", trim = TRUE) +
  # scale_y_log10() +
  theme_minimal(base_size = 25) +
  scale_fill_manual(values = c("BenchDrop-seq" = "#F8766D")) +
  labs(
    x = NULL,
    y = "UMIs per cell",
    title = "UMIs per Cell"
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 20)
  )

# Visualize number of genes per cell
genes.per.cell <- ggplot(
  tx.meta,
  aes(x = sample, y = nFeature_RNA, fill = sample)
) +
  geom_violin(scale = "width", trim = TRUE) +
  # scale_y_log10() +
  theme_minimal(base_size = 25) +
  scale_fill_manual(values = c("BenchDrop-seq" = "#00BFC4")) +
  labs(
    x = NULL,
    y = "Genes per cell",
    title = "Genes per Cell"
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 20)
  )

# Visualize number of transcripts per cell
txs.per.cell <- ggplot(
  tx.meta,
  aes(x = sample, y = nFeature_transcript, fill = sample)
) +
  geom_violin(scale = "width", trim = TRUE) +
  # scale_y_log10() +
  theme_minimal(base_size = 25) +
  scale_fill_manual(values = c("BenchDrop-seq" = "gray")) +
  labs(
    x = NULL,
    y = "Transcripts per cell",
    title = "Transcripts per Cell"
  ) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 20)
  )

# Visualize together
umis.per.cell + genes.per.cell + txs.per.cell



#### SUPPLEMENTARY FIG. 2C (SHORT VS BULK CORRELATION) ####
# Load data
sr.pbmc <- readRDS("SR.PBMC.S3.rds")

# Load bulk data
bulk.tx.quants <- read.table("quant.sf", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract columns of interest
bulk.tx.quants <- bulk.tx.quants %>% dplyr::select(transcript.id = Name, TPM)
bulk.tx.quants <- setNames(bulk.tx.quants$TPM, bulk.tx.quants$transcript.id)

# Read in transcript to gene mapping
t2g <- read.delim("t2gnames.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(t2g) <- c("transcript.id", "gene.name")

# Align names
t2g.sub <- t2g[t2g$transcript.id %in% names(bulk.tx.quants), ]

# Aggregate bulk transcript TPM counts to bulk gene TPM counts
bulk.gene.quants <- tapply(bulk.tx.quants[t2g.sub$transcript.id], t2g.sub$gene.name, sum)

# Access SR gene matrix
sr.gene.mat <- GetAssayData(sr.pbmc, assay = "RNA", layer = "data")

# Compute pseudobulk gene counts 
sr.gene.rowsums <- rowSums(sr.gene.mat)
sr.gene.bulk <- sr.gene.rowsums*1e6 / sum(sr.gene.rowsums)

# Subset for ALL genes
all.genes <- union(names(sr.gene.bulk), names(bulk.gene.quants))
sr.gene.bulk.full <- setNames(numeric(length(all.genes)), all.genes)
bulk.gene.quants.full <- setNames(numeric(length(all.genes)), all.genes)
sr.gene.bulk.full[names(sr.gene.bulk)] <- sr.gene.bulk
bulk.gene.quants.full[names(bulk.gene.quants)] <- bulk.gene.quants
sr.gene.bulk <- sr.gene.bulk.full
bulk.gene.quants <- bulk.gene.quants.full

# Compute correlations
pearson.cor <- cor(sr.gene.bulk, bulk.gene.quants, method = "pearson")
spearman.cor <- cor(sr.gene.bulk, bulk.gene.quants, method = "spearman")

# Create DF for plotting
df.cor <- data.frame(Short = log1p(sr.gene.bulk), Bulk = log1p(bulk.gene.quants))

# Plot with hex density
p <- ggplot(df.cor, aes(x = Bulk, y = Short)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", trans = "log10") + 
  labs(x = "Bulk", 
       y = "Short", 
       fill = "Density") +
  theme_minimal(base_size = 25) +
  theme(axis.line = element_line(color = "black", linewidth = 0.8),
        panel.grid = element_blank())

p + annotate("text", 
             x = min(df.cor$Short, na.rm = TRUE), 
             y = max((df.cor$Short + 1), na.rm = TRUE), 
             hjust = 0, vjust = 1, size = 6,
             label = paste0("Spearman = ", round(spearman.cor, 3)))



#### SUPPLEMENTARY FIG. 2D (LONG VS BULK CORRELATION) ####
# Load data
lr.pbmc <- readRDS("LR.PBMC.S3.rds")

# Load bulk data
bulk.tx.quants <- read.table("quant.sf", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract columns of interest
bulk.tx.quants <- bulk.tx.quants %>% dplyr::select(transcript.id = Name, TPM)
bulk.tx.quants <- setNames(bulk.tx.quants$TPM, bulk.tx.quants$transcript.id)

# Read in transcript to gene mapping
t2g <- read.delim("t2gnames.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(t2g) <- c("transcript.id", "gene.name")

# Align names
t2g.sub <- t2g[t2g$transcript.id %in% names(bulk.tx.quants), ]

# Aggregate bulk transcript TPM counts to bulk gene TPM counts
bulk.gene.quants <- tapply(bulk.tx.quants[t2g.sub$transcript.id], t2g.sub$gene.name, sum)

# Access SR gene matrix
lr.gene.mat <- GetAssayData(lr.pbmc, assay = "RNA", layer = "data")

# Compute pseudobulk gene counts 
lr.gene.rowsums <- rowSums(lr.gene.mat)
lr.gene.bulk <- lr.gene.rowsums*1e6 / sum(lr.gene.rowsums)

# Subset for ALL genes
all.genes <- union(names(lr.gene.bulk), names(bulk.gene.quants))
lr.gene.bulk.full <- setNames(numeric(length(all.genes)), all.genes)
bulk.gene.quants.full <- setNames(numeric(length(all.genes)), all.genes)
lr.gene.bulk.full[names(lr.gene.bulk)] <- lr.gene.bulk
bulk.gene.quants.full[names(bulk.gene.quants)] <- bulk.gene.quants
lr.gene.bulk <- lr.gene.bulk.full
bulk.gene.quants <- bulk.gene.quants.full

# Compute correlations
pearson.cor <- cor(lr.gene.bulk, bulk.gene.quants, method = "pearson")
spearman.cor <- cor(lr.gene.bulk, bulk.gene.quants, method = "spearman")

# Create DF for plotting
df.cor <- data.frame(`BenchDrop-seq` = log1p(lr.gene.bulk), Bulk = log1p(bulk.gene.quants), check.names = FALSE)

# Plot with hex density
p <- ggplot(df.cor, aes(x = Bulk, y = `BenchDrop-seq`)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", trans = "log10") + 
  labs(x = "Bulk", 
       y = "BenchDrop-seq", 
       fill = "Density") +
  theme_minimal(base_size = 25) +
  theme(axis.line = element_line(color = "black", linewidth = 0.8),
        panel.grid = element_blank())

p + annotate("text", 
             x = min(df.cor$`BenchDrop-seq`, na.rm = TRUE), 
             y = max((df.cor$`BenchDrop-seq` + 1), na.rm = TRUE), 
             hjust = 0, vjust = 1, size = 6,
             label = paste0("Spearman = ", round(spearman.cor, 3)))



#### SUPPLEMENTARY FIG. 2E (SHORT VS LONG CORRELATION) ####
# Load data
sr.pbmc <- readRDS("SR.PBMC.S3.rds")
lr.pbmc <- readRDS("LR.PBMC.S3.rds")

# Access SR gene matrix
sr.gene.mat <- GetAssayData(sr.pbmc, assay = "RNA", layer = "data")

# Access LR gene matrix
lr.gene.mat <- GetAssayData(lr.pbmc, assay = "RNA", layer = "data")

# Pseudobulk each matrix
sr.gene.bulk <- rowSums(sr.gene.mat) 
lr.gene.bulk <- rowSums(lr.gene.mat)

# Subset for ALL genes
all.genes <- union(names(sr.gene.bulk), names(lr.gene.bulk))
sr.gene.bulk.full <- setNames(numeric(length(all.genes)), all.genes)
lr.gene.bulk.full <- setNames(numeric(length(all.genes)), all.genes)
sr.gene.bulk.full[names(sr.gene.bulk)] <- sr.gene.bulk
lr.gene.bulk.full[names(lr.gene.bulk)] <- lr.gene.bulk
sr.gene.bulk <- sr.gene.bulk.full
lr.gene.bulk <- lr.gene.bulk.full

# Compute correlations
pearson.cor <- cor(sr.gene.bulk, lr.gene.bulk, method = "pearson")
spearman.cor <- cor(sr.gene.bulk, lr.gene.bulk, method = "spearman")

# Create DF for plotting
df.cor <- data.frame(Short = log1p(sr.gene.bulk), `BenchDrop-seq` = log1p(lr.gene.bulk), check.names = FALSE)

# Plot with hex density
p <- ggplot(df.cor, aes(x = Short, y = `BenchDrop-seq`)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", trans = "log10") + 
  labs(x = "Short", 
       y = "BenchDrop-seq",
       fill = "Density") +
  theme_minimal(base_size = 18) +
  theme(axis.line = element_line(color = "black", linewidth = 0.8),
        panel.grid = element_blank())

p + annotate("text", 
             x = min(df.cor$`BenchDrop-seq`, na.rm = TRUE), 
             y = max((df.cor$`BenchDrop-seq` + 1), na.rm = TRUE), 
             hjust = 0, vjust = 1, size = 6,
             label = paste0("Spearman = ", round(spearman.cor, 3)))