# Figure 2

# Load packages
library(Seurat)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(dplyr)
library(patchwork)
library(pbapply)
library(reshape2)
library(viridis)

#### FIG. 2B (SINGLE-CELL METRICS) ####
# Load Seurat objects
sr.k562 <- readRDS("SR.K562.S3.rds")
lr.k562 <- readRDS("LR.K562.S3.rds")

# Add metadata label & combine
sr.k562$sample <- "Short"
lr.k562$sample <- "BenchDrop-seq"
metadata <- bind_rows(
  sr.k562@meta.data %>% dplyr::select(nCount_RNA, nFeature_RNA, sample),
  lr.k562@meta.data %>% dplyr::select(nCount_RNA, nFeature_RNA, sample),
)

# Visualize number of UMIs per cell
umis.cell <- ggplot(metadata, aes(x = sample, y = nCount_RNA, fill = sample)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_y_log10() +
  theme_minimal(base_size = 25) +
  labs(x = "Sample", 
       y = "UMIs per cell") +
  theme(legend.position = "none", 
        panel.grid = element_blank())

# Visualize number of genes per cell
genes.cell <- ggplot(metadata, aes(x = sample, y = nFeature_RNA, fill = sample)) +
  geom_violin(scale = "width", trim = TRUE) +
  scale_y_log10() +
  theme_minimal(base_size = 25) +
  labs(x = "Sample", 
       y = "Genes per cell") +
  theme(legend.position = "none", 
        panel.grid = element_blank())

# Plot together
umis.cell + genes.cell



#### FIG. 2C (EQUIVALENCE CLASS LENGTH) ####
# Load length data
lf <- read.table("lr.lens.txt")
sf <- read.table("sr.lens.txt")

# Store data as DataFrame
lf <- data.frame(table(lf$V1))
sf <- data.frame(table(sf$V1))

# Add platform label
lf$kind <- "BenchDrop-seq"
sf$kind <- "Short"

# Compute percent
lf$percent <- lf$Freq * 100.0 / sum(lf$Freq)
sf$percent <- sf$Freq * 100.0 / sum(sf$Freq)

# Combine DF for plotting
df <- rbind(lf, sf)

# Visualize
ggplot(df, aes(x = factor(kind), y = Var1, group = factor(kind), color = kind, size = percent)) +
  geom_point() +
  theme_classic() +
  labs(title = "Equivalence Class Length",
       x = "Platform", 
       y = "Length") + 
  scale_size_continuous(range = c(2, 12)) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))



#### FIG. 2E (SR vs BULK CORRELATION) #####
# Load data
sr.k562 <- readRDS("SR.K562.S3.rds")
lr.k562 <- readRDS("LR.K562.S3.rds")

# Load bulk data
bulk.tx.quants <- read.table("k562.pe.bulk.quant.sf", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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
sr.gene.mat <- GetAssayData(sr.k562, assay = "RNA", layer = "data")

# Compute pseudobulk gene counts 
sr.gene.rowsums <- rowSums(sr.gene.mat)
sr.gene.bulk <- sr.gene.rowsums*1e6 / sum(sr.gene.rowsums)

# Initialize vectors with all transcript names
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
        panel.grid = element_blank(), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p + annotate("text", 
             x = min(df.cor$Short, na.rm = TRUE), 
             y = max(df.cor$Bulk, na.rm = TRUE), 
             hjust = 0, vjust = 1, size = 6,
             label = paste0("Spearman = ", round(spearman.cor, 3)))



#### FIG. 2F (LR vs BULK CORRELATION) #####
# Load data
sr.k562 <- readRDS("SR.K562.S3.rds")
lr.k562 <- readRDS("LR.K562.S3.rds")

# Load bulk data
bulk.tx.quants <- read.table("k562.pe.bulk.quant.sf", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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

# Access LR gene matrix
lr.gene.mat <- GetAssayData(lr.k562, assay = "RNA", layer = "data")

# Compute pseudobulk gene counts 
lr.gene.rowsums <- rowSums(lr.gene.mat)
lr.gene.bulk <- lr.gene.rowsums*1e6 / sum(lr.gene.rowsums)

# Initialize vectors with all transcript names
all.genes <- union(names(lr.gene.bulk), names(bulk.gene.quants))
lr.gene.bulk.full <- setNames(numeric(length(all.genes)), all.genes)
bulk.gene.quants.full <- setNames(numeric(length(all.genes)), all.genes)

# Fill values where present
lr.gene.bulk.full[names(lr.gene.bulk)] <- lr.gene.bulk
bulk.gene.quants.full[names(bulk.gene.quants)] <- bulk.gene.quants

# Replace originals
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
        panel.grid = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p + annotate("text", 
             x = min(df.cor$`BenchDrop-seq`, na.rm = TRUE), 
             y = max(df.cor$Bulk, na.rm = TRUE), 
             hjust = 0, vjust = 1, size = 6,
             label = paste0("Spearman = ", round(spearman.cor, 3)))



##### FIG. 2G (SR vs LR CORRELATION) #####
# Load data
sr.k562 <- readRDS("SR.K562.S3.rds")
lr.k562 <- readRDS("LR.K562.S3.rds")

# Access SR gene matrix
sr.gene.mat <- GetAssayData(sr.k562, assay = "RNA", layer = "data")

# Access LR gene matrix
lr.gene.mat <- GetAssayData(lr.k562, assay = "RNA", layer = "data")

# Pseudobulk each matrix
sr.gene.bulk <- rowSums(sr.gene.mat) 
lr.gene.bulk <- rowSums(lr.gene.mat)

# Initialize vectors with all transcript names
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
       y = "BenchDrop-seq") +
  theme_minimal(base_size = 25) +
  theme(axis.line = element_line(color = "black", linewidth = 0.8),
        panel.grid = element_blank(), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

p + annotate("text", 
             x = min(df.cor$Short, na.rm = TRUE), 
             y = max(df.cor$`BenchDrop-seq`, na.rm = TRUE), 
             hjust = 0, vjust = 1, size = 6,
             label = paste0("Spearman = ", round(spearman.cor, 3)))



#### FIG. 2H (DOWNSAMPLING) ####
# Load SR data
sfmat <- Read10X("sensitivity_3")

# Read in transcript to gene mapping
fread <- read.table("t2gnames.txt", header = FALSE, stringsAsFactors = FALSE)
t2g <- as.vector(fread$V2)
names(t2g) <- as.vector(fread$V1)

# Initialize paths to downsampled data
fpaths <- c("v3/",
            paste0("Downsampling/", seq(10, 90, 10), "/"))

# Function to load matrix
GetMat <- function(bpath, keep.cells) {
  mat <- readMM(paste0(bpath, "matrix.mtx"))
  bcs <- read.table(paste0(bpath, "barcodes.tsv"))$V1
  feats <- read.table(paste0(bpath, "features.tsv"))$V1
  mat <- t(mat)
  
  rownames(mat) <- feats
  colnames(mat) <- bcs
  mat <- mat[, keep.cells]
  gmat <- aggregate.Matrix(mat, groupings = as.character(t2g[rownames(mat)]), 
                           fun = "sum", MARGIN = 2)
  gmat
}

# Function to access matrix from each level of downsampling
sub.mats <- lapply(fpaths, function(x) {
  print(x)
  GetMat(x, colnames(sfmat))
})

# Access all transcripts and format into a list and remove duplicates
all.txps <- unlist(lapply(sub.mats, function(x) rownames(x)))
all.txps <- all.txps[!duplicated(all.txps)]

# Add transcripts from SR K562 object
all.txps <- union(all.txps, rownames(sfmat))

# Loop over each cell
corrs <- pblapply(colnames(sfmat), function(cname){
  
  # Initialize a DF from all transcripts
  df <- data.frame(all.txps)
  
  # Initialize columns 
  # df$Short <- 0
  rownames(df) <- df$all.txps
  df$all.txps <- NULL
  df$`BenchDrop-seq` <- 0
  df$"0.1" <- 0
  df$"0.2" <- 0
  df$"0.3" <- 0
  df$"0.4" <- 0
  df$"0.5" <- 0
  df$"0.6" <- 0
  df$"0.7" <- 0
  df$"0.8" <- 0
  df$"0.9" <- 0
  
  # Fill in each transcript value for that level lf of downsampling
  df[rownames(sub.mats[[1]]),"BenchDrop-seq"] <-  sub.mats[[1]][,cname]
  df[rownames(sub.mats[[2]]),"0.1"] <-  sub.mats[[2]][,cname]
  df[rownames(sub.mats[[3]]),"0.2"] <-  sub.mats[[3]][,cname]
  df[rownames(sub.mats[[4]]),"0.3"] <-  sub.mats[[4]][,cname]
  df[rownames(sub.mats[[5]]),"0.4"] <-  sub.mats[[5]][,cname]
  df[rownames(sub.mats[[6]]),"0.5"] <-  sub.mats[[6]][,cname]
  df[rownames(sub.mats[[7]]),"0.6"] <-  sub.mats[[7]][,cname]
  df[rownames(sub.mats[[8]]),"0.7"] <-  sub.mats[[8]][,cname]
  df[rownames(sub.mats[[9]]),"0.8"] <-  sub.mats[[9]][,cname]
  df[rownames(sub.mats[[10]]),"0.9"] <- sub.mats[[10]][,cname]
  # df[rownames(sfmat),"Short"] <-  sfmat[,cname]
  
  # Compute correlation
  cor(df, method = "pearson")[1,]
})

# Combine results across downsampling levels and reshape for plotting
corrs <- do.call(rbind, corrs)
corrs <- reshape2::melt(corrs)
colnames(corrs) <- c("var1", "Subsampling", "Correlation")

# Order the downsampling levels
corrs$Subsampling <- factor(
  corrs$Subsampling,
  # levels = c("Short", "BenchDrop-seq", "0.9", "0.8", "0.7", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1"))
  levels = c("BenchDrop-seq", "0.9", "0.8", "0.7", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1"))

# Plot
ggplot(corrs, aes(x = Subsampling, y = Correlation, fill = Subsampling)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_minimal(base_size = 25) +
  labs(x = "Subsampling Level",
       y = "Correlation") + 
  theme(
    panel.grid = element_blank(), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20) 
  )