# Supp Figure 2

# Load packages
library(Seurat)
library(Matrix)
library(ggplot2)
library(viridis)

#### SUPP FIG.1D ####
# Load LR data
lr.k562 <- readRDS("LR.K562.S3.rds")
lr.mat <- GetAssayData(lr.k562, assay = "RNA", layer = "data")

# Load Alevin quants
alevin.mat <- readMM("quants_mat.mtx")
alevin.mat <- t(alevin.mat)
alevin.feats <- read.table("/quants_mat_cols.txt")$V1
rownames(alevin.mat) <- alevin.feats
alevin.bcs <- read.table("quants_mat_rows.txt")$V1
colnames(alevin.mat) <- alevin.bcs

# Subset Alevin quants for the same cells as BenchDrop-seq
common.cols <- intersect(colnames(lr.mat), colnames(alevin.mat))
alevin.mat <- alevin.mat[, common.cols, drop = FALSE]

# Compute pseudobulk gene counts
lr.mat.bulk <- rowSums(lr.mat)
alevin.bulk <- rowSums(alevin.mat)

# Initialize vectors with all transcript names
all.genes <- union(names(lr.mat.bulk), names(alevin.bulk))
lr.mat.bulk.full <- setNames(numeric(length(all.genes)), all.genes)
alevin.bulk.full <- setNames(numeric(length(all.genes)), all.genes)
lr.mat.bulk.full[names(lr.mat.bulk)] <- lr.mat.bulk
alevin.bulk.full[names(alevin.bulk)] <- alevin.bulk
lr.mat.bulk <- lr.mat.bulk.full
alevin.bulk <- alevin.bulk.full

# Compute correlations
pearson.cor <- cor(lr.mat.bulk, alevin.bulk, method = "pearson")
spearman.cor <- cor(lr.mat.bulk, alevin.bulk, method = "spearman")

# DF for plotting
df.cor <- data.frame(`BenchDrop-seq` = log1p(lr.mat.bulk), `Salmon-Alevin` = log1p(alevin.bulk), check.names = FALSE)

# Plot
p <- ggplot(df.cor, aes(x = `Salmon-Alevin`, y = `BenchDrop-seq`)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", trans = "log10") + 
  labs(x = "Salmon-Alevin", 
       y = "BenchDrop-seq", 
       fill = "Density") +
  theme_minimal(base_size = 25) +
  theme(axis.line = element_line(color = "black", linewidth = 0.8),
        panel.grid = element_blank())

p + annotate("text", 
             x = min(df.cor$`BenchDrop-seq`, na.rm = TRUE), 
             y = max(df.cor$`Salmon-Alevin`, na.rm = TRUE), 
             hjust = 0, vjust = 1, size = 6,
             label = paste0("Spearman = ", round(spearman.cor, 3)))



#### SUPP FIG.1D ####
# Load LR data
lr.k562 <- readRDS("LR.K562.S3.rds")
lr.mat <- GetAssayData(lr.k562, assay = "RNA", layer = "data")

# Load Kallisto-Bustools quants
kb.mat <- readMM("cells_x_genes.mtx")
kb.mat <- t(kb.mat)
kb.feats <- read.table("cells_x_genes.genes.names.txt")$V1
rownames(kb.mat) <- kb.feats
kb.bcs <- read.table("cells_x_genes.barcodes.txt")$V1
colnames(kb.mat) <- kb.bcs

# Subset Kallisto-Bustools quants for the same cells as BenchDrop-seq
common.cols <- intersect(colnames(lr.mat), colnames(kb.mat))
kb.mat <- kb.mat[, common.cols, drop = FALSE]

# Compute pseudobulk gene counts
lr.mat.bulk <- rowSums(lr.mat)
kb.bulk <- rowSums(kb.mat)

# Initialize vectors with all transcript names
all.genes <- union(names(lr.mat.bulk), names(kb.bulk))
lr.mat.bulk.full <- setNames(numeric(length(all.genes)), all.genes)
kb.bulk.full <- setNames(numeric(length(all.genes)), all.genes)
lr.mat.bulk.full[names(lr.mat.bulk)] <- lr.mat.bulk
kb.bulk.full[names(kb.bulk)] <- kb.bulk
lr.mat.bulk <- lr.mat.bulk.full
kb.bulk <- kb.bulk.full

# Compute correlations
pearson.cor <- cor(lr.mat.bulk, kb.bulk, method = "pearson")
spearman.cor <- cor(lr.mat.bulk, kb.bulk, method = "spearman")

# DF for plotting
df.cor <- data.frame(`BenchDrop-seq` = log1p(lr.mat.bulk), `Kallisto-Bustools` = log1p(kb.bulk), check.names = FALSE)

# Plot
p <- ggplot(df.cor, aes(x = `Kallisto-Bustools`, y = `BenchDrop-seq`)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", trans = "log10") + 
  labs(x = "Kallisto-Bustools", 
       y = "BenchDrop-seq", 
       fill = "Density") +
  theme_minimal(base_size = 25) +
  theme(axis.line = element_line(color = "black", linewidth = 0.8),
        panel.grid = element_blank())

p + annotate("text", 
             x = min(df.cor$`BenchDrop-seq`, na.rm = TRUE), 
             y = max(df.cor$`Kallisto-Bustools`, na.rm = TRUE), 
             hjust = 0, vjust = 1, size = 6,
             label = paste0("Spearman = ", round(spearman.cor, 3)))
