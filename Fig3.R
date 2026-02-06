# Figure 3

# Load packages
rm(list = ls())
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

#### FIG. 3A - 3C ####
source("Fig3.Source")

{# Load data
  sr.pip <- GetObject("sr.pbmc.sens3")
  sr <- Getcpm(sr.pip)
  
  lr <- GetObject("lr.pbmc.95")
  lr <- subset(lr, cells = Cells(sr.pip))
  lr <- NormalizeData(lr, normalization.method = "LogNormalize")
  lr <- Gettpm(lr, kind = "lr")
  
  bulk <- GetObject("bulk.pbmc")
  bulk <- Gettpm(bulk, kind = "bulk")
}

{ 
  keep.genes <- intersect(
    intersect(names(lr$cpm), 
              names(sr)), 
    names(bulk$cpm))
  length(keep.genes)
  
  df <- data.frame(
    "long" = lr$cpm[keep.genes],
    "short" = sr[keep.genes],
    "bulk" = bulk$cpm[keep.genes]
  )
  cor(df, method = "spearman")
  df <- FillGeneMetadata(df)
  
  df$fc <- log((df$short + 1) / (df$long + 1))
  df <- df[order(df$fc),]
  
  sdf <- df#[df$type == "protein_coding",]
  sdf$gname <- rownames(sdf)
  sdf$sfc <- log10((sdf$short + 1) / (sdf$bulk + 1))
  sdf$lfc <- log10((sdf$long + 1) / (sdf$bulk + 1))
  
  # Add metadata
  gtf <- import("genes.gtf")
  
  # GENE LEVEL INFORMATION
  genes <- gtf[gtf$type == "gene"]
  gene.info <- data.frame(
    gname = mcols(genes)$gene_name,
    gene_id = mcols(genes)$gene_id,
    chr = seqnames(genes),
    gene_length = width(ranges(genes))
  ) %>% 
    distinct(gname, .keep_all = TRUE)
  
  # TRANSCRIPT LEVEL INFORMATION
  transcripts <- gtf[gtf$type == "transcript"]
  transcript.info <- data.frame(
    transcript_id = mcols(transcripts)$transcript_id,
    gname = mcols(transcripts)$gene_name,
    transcript_length = width(ranges(transcripts))
  )
  
  # EXON LEVEL INFORMATION
  exons <- gtf[gtf$type == "exon"]
  exon.info <- data.frame(
    transcript_id = mcols(exons)$transcript_id,
    exon_length = width(ranges(exons))
  )
  
  # Sum of exon lengths per transcript
  exon.sum <- exon.info %>%
    group_by(transcript_id) %>%
    summarise(total_exon_length = sum(exon_length, na.rm = TRUE)) %>%
    ungroup()
  
  # Compute intron lengths per transcript (intron = transcript_length - sum(exon_lengths))
  transcript.full <- transcript.info %>%
    left_join(exon.sum, by = "transcript_id") %>%
    mutate(
      intron_length = transcript_length - total_exon_length
    )
  
  # Average transcript and intron metrics per gene
  transcript.summary <- transcript.full %>%
    group_by(gname) %>%
    summarise(
      n_transcripts = n(),
      mean_transcript_length = mean(transcript_length, na.rm = TRUE),
      mean_exon_length = mean(total_exon_length, na.rm = TRUE),
      mean_intron_length = mean(intron_length, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Merge all information
  gene.info <- gene.info %>%
    left_join(transcript.summary, by = "gname")
  
  # Add to sdf DF
  sdf <- sdf %>%
    left_join(
      gene.info %>% 
        dplyr::select(gname, chr, gene_length,
                      n_transcripts, mean_transcript_length,
                      mean_exon_length, mean_intron_length),
      by = "gname"
    )
  
  # Remove duplicate rows
  sdf <- sdf[!duplicated(sdf), ]
  
  # Remove circle around origin
  radius <- 0.25
  sdf.p1 <- sdf[sqrt(sdf$lfc^2 + sdf$fc^2) > radius, ]
  sdf.p2 <- sdf[sqrt(sdf$sfc^2 + sdf$fc^2) > radius, ]
  sdf.p3 <- sdf[sqrt(sdf$lfc^2 + sdf$sfc^2) > radius, ]
  
  # Metrics
  genes_plotted_unique <- bind_rows(
    sdf.p1 %>% dplyr::select(gname),
    sdf.p2 %>% dplyr::select(gname),
    sdf.p3 %>% dplyr::select(gname)
  ) %>%
    distinct(gname)
  
  n_genes_plotted_unique <- nrow(genes_plotted_unique)

  genes_removed_unique <- sdf %>%
    filter(
      sqrt(lfc^2 + fc^2) <= radius &
        sqrt(sfc^2 + fc^2) <= radius &
        sqrt(lfc^2 + sfc^2) <= radius
    ) %>%
    distinct(gname)
  
  n_genes_removed_unique <- nrow(genes_removed_unique)
  
  # Final values
  cat("Unique genes plotted:", n_genes_plotted_unique, "\n")
  cat("Unique genes removed:", n_genes_removed_unique, "\n")
}

# Plot Fig. 3A 
p1 <- ggplot(sdf.p1, aes(y = fc, x = lfc, label = gname)) +
  geom_point() +
  geom_hex(bins = 100) +
  ylim(-8,8)+
  xlim(-8,8)+
  scale_fill_viridis_c(option = "plasma", trans = "log10") +
  geom_label_repel(
    aes(label = ifelse(gname %in% c("HIST1H4F", "HIST2H2AB",
                                    "RPS4Y1", "DDX3Y",
                                    "PTPRCAP", "PTPRD", "RBFOX1", "CORO1B"), as.character(gname), "")),
    size = 7, 
    fontface = "bold",
    box.padding   = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps  = Inf
  ) +
  theme_minimal()

p1 <- ggMarginal(p1, yparams = list(fill = "#5C1A8A", color = "black"),
                 xparams = list(fill = "#D5722A", color = "black")) # Customize x-axis histogram

p2 <- ggplot(sdf.p2, aes(y=fc, x=sfc, label = gname)) +
  geom_point() +
  geom_hex(bins = 100) +
  ylim(-8,8)+
  xlim(-8,8)+
  scale_fill_viridis_c(option = "plasma", trans = "log10") +
  geom_label_repel(
    aes(label = ifelse(gname %in% c("HIST1H4F", "HIST2H2AB", 
                                    "RPS4Y1", "DDX3Y",
                                    "PTPRCAP", "PTPRD", "RBFOX1", "CORO1B"), as.character(gname), "")),
    size = 7, 
    fontface = "bold",
    box.padding   = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps  = Inf
  ) +
  theme_minimal()

p2 <- ggMarginal(p2, yparams = list(fill = "#5C1A8A", color = "black"),
                 xparams = list(fill = "#D5722A", color = "black")) # Customize x-axis histogram
wrap_plots(list(p1, p2))

p3 <- ggplot(sdf.p3, aes(y=lfc, x=sfc, label = gname)) +
  geom_point() +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(option = "plasma", trans = "log10") +
  geom_label_repel(
    aes(label = ifelse(gname %in% c("HIST1H4F", "HIST2H2AB",
                                    "RPS4Y1", "DDX3Y",
                                    "PTPRCAP", "PTPRD", "RBFOX1", "CORO1B"), as.character(gname), "")),
    size = 7, 
    fontface = "bold",
    box.padding   = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps  = Inf
  ) +
  theme_minimal()

ggMarginal(p3, yparams = list(fill = "#5C1A8A", color = "black"),
           xparams = list(fill = "#D5722A", color = "black")) # Customize x-axis histogram

# Plot Fig. 3B (log(mean intron length))
red_genes <- c("PTPRD", "RBFOX1")
col_red  <- "black"  

add_label_colors_red <- function(df) {
  df %>%
    mutate(
      label = ifelse(gname %in% red_genes, gname, ""),
      label_col = ifelse(gname %in% red_genes, col_red, "black")
    )
}

sdf.p1 <- add_label_colors_red(sdf.p1)
sdf.p2 <- add_label_colors_red(sdf.p2)
sdf.p3 <- add_label_colors_red(sdf.p3)

p1 <- ggplot(sdf.p1, aes(y=fc, x=lfc)) +
  geom_point(alpha = 0.5) +
  stat_summary_hex(aes(z = log10(mean_intron_length)), fun = mean, bins = 100) +
  scale_fill_viridis_c(option = "plasma", name = " ") +
  ylim(-8,8) + xlim(-8,8) +
  theme_minimal() +
  geom_label_repel(
    aes(label = label, color = label_col),
    size = 7,
    fontface = "bold",
    box.padding   = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps  = Inf,
    na.rm = TRUE
  ) +
  scale_color_identity() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.key.size = unit(1.2, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

p2 <- ggplot(sdf.p2, aes(y=fc, x=sfc)) +
  geom_point(alpha = 0.5) +
  stat_summary_hex(aes(z = log10(mean_intron_length)), fun = mean, bins = 100) +
  scale_fill_viridis_c(option = "plasma", name = " ") +
  ylim(-8,8) + xlim(-8,8) +
  theme_minimal() +
  geom_label_repel(
    aes(label = label, color = label_col),
    size = 7,
    fontface = "bold",
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  scale_color_identity() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.key.size = unit(1.2, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

p3 <- ggplot(sdf.p3, aes(y=lfc, x=sfc)) +
  geom_point(alpha = 0.5) +
  stat_summary_hex(aes(z = log10(mean_intron_length)), fun = mean, bins = 100) +
  scale_fill_viridis_c(option = "plasma", name = " ") +
  geom_abline(slope = 1, intercept = 0, color = "black", alpha = 0.6, linetype = "dashed") +
  theme_minimal() +
  geom_label_repel(
    aes(label = label, color = label_col),
    size = 7,
    fontface = "bold",
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  scale_color_identity() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.key.size = unit(1.2, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

p3

# Plot Fig. 3B (log(gene length))
red_genes <- c("PTPRD", "RBFOX1")
col_red  <- "black"   # same red as before

add_label_colors_red <- function(df) {
  df %>%
    mutate(
      label = ifelse(gname %in% red_genes, gname, ""),
      label_col = ifelse(gname %in% red_genes, col_red, "black")
    )
}

sdf.p1 <- add_label_colors_red(sdf.p1)
sdf.p2 <- add_label_colors_red(sdf.p2)
sdf.p3 <- add_label_colors_red(sdf.p3)

p1 <- ggplot(sdf.p1, aes(y=fc, x=lfc)) +
  geom_point(alpha = 0.5) +
  stat_summary_hex(aes(z = log10(gene_length)), fun = mean, bins = 100) +
  scale_fill_viridis_c(option = "plasma", name = " ") +
  ylim(-8,8) + xlim(-8,8) +
  theme_minimal() +
  geom_label_repel(
    aes(label = label, color = label_col),
    size = 7,
    fontface = "bold",
    box.padding   = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps  = Inf,
    na.rm = TRUE
  ) +
  scale_color_identity() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.key.size = unit(1.2, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

p2 <- ggplot(sdf.p2, aes(y=fc, x=sfc)) +
  geom_point(alpha = 0.5) +
  stat_summary_hex(aes(z = log10(gene_length)), fun = mean, bins = 100) +
  scale_fill_viridis_c(option = "plasma", name = " ") +
  ylim(-8,8) + xlim(-8,8) +
  theme_minimal() +
  geom_label_repel(
    aes(label = label, color = label_col),
    size = 7,
    fontface = "bold",
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  scale_color_identity() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.key.size = unit(1.2, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

p3 <- ggplot(sdf.p3, aes(y=lfc, x=sfc)) +
  geom_point(alpha = 0.5) +
  stat_summary_hex(aes(z = log10(gene_length)), fun = mean, bins = 100) +
  scale_fill_viridis_c(option = "plasma", name = " ") +
  geom_abline(slope = 1, intercept = 0, color = "black", alpha = 0.6, linetype = "dashed") +
  theme_minimal() +
  geom_label_repel(
    aes(label = label, color = label_col),
    size = 7,
    fontface = "bold",
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps = Inf,
    na.rm = TRUE
  ) +
  scale_color_identity() +
  theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom"),
    legend.key.size = unit(1.2, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

p3



#### FIG. 3E ####
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

# Plot
lr.gene.umap <- DimPlot(lr.pbmc, reduction = "ref.umap", group.by = "celltype.grouped", label = F, repel = T, label.size = 5.5, pt.size = 1.5) + 
  NoLegend() + 
  xlab("UMAP_1") +
  ylab("UMAP_2") + 
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
lr.gene.umap <- LabelClusters(lr.gene.umap, id = "celltype.grouped", fontface = "bold", size = 5.5, repel = TRUE, color = "black")
lr.gene.umap



#### FIG. 3F #####
# Load data
lr.pbmc <- readRDS("LR.PBMC.S3.rds")

# Create a new metadata column to store celltypes
lr.pbmc$celltype.grouped <- lr.pbmc$predicted.celltype.l2
lr.pbmc$celltype.grouped[lr.pbmc$celltype.grouped %in% c("cDC1", "cDC2")] <- "DC"
lr.pbmc$celltype.grouped[lr.pbmc$celltype.grouped %in% c("NK", "NK Proliferating", "NK_CD56bright")] <- "NK"
lr.pbmc$celltype.grouped[lr.pbmc$celltype.grouped %in% c("B memory", "B naive", "B intermediate")] <- "B Cells"
lr.pbmc$celltype.grouped[lr.pbmc$celltype.grouped %in% c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL")] <- "CD4 T Cells"
lr.pbmc$celltype.grouped[lr.pbmc$celltype.grouped %in% c("CD8 Naive", "CD8 TCM", "CD8 TEM")] <- "CD8 T Cells"
lr.pbmc$celltype.grouped[lr.pbmc$celltype.grouped %in% c("Eryth")] <- "Erythroid"
lr.pbmc$celltype.grouped[lr.pbmc$celltype.grouped %in% c("dnT", "ILC")] <- NA
lr.pbmc <- subset(lr.pbmc, !is.na(celltype.grouped))

# Set identities
Idents(lr.pbmc) <- "celltype.grouped"

# Find gene-level markers (positive only)
gene.markers <- FindAllMarkers(lr.pbmc, assay = "RNA", only.pos = TRUE)

# Filter to highly up-regulated genes
gene.markers.filtered <- gene.markers %>%
  group_by(cluster) %>%
  mutate(spec.score = avg_log2FC * (pct.1 - pct.2)) %>%
  dplyr::filter(avg_log2FC > 0.4, pct.1 > 0.25, pct.2 < 0.15)

# Top 3 unique genes per cell type
top3.unique <- data.frame()
used.genes <- c()

# Loop
for (clust in unique(gene.markers.filtered$cluster)) {
  
  df <- gene.markers.filtered %>%
    dplyr::filter(cluster == clust) %>%
    arrange(desc(spec.score))
  
  df <- df[!df$gene %in% used.genes, ]
  top_n <- head(df, 3)
  top3.unique <- rbind(top3.unique, top_n)
  
  used.genes <- c(used.genes, top_n$gene)
}

# Scale data
genes.to.plot <- top3.unique$gene
lr.pbmc <- ScaleData(lr.pbmc, features = genes.to.plot, assay = "RNA")

# Plot
DotPlot(lr.pbmc, features = genes.to.plot, assay = "RNA", dot.scale = 10) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        plot.margin = unit(c(1, 1, 1, 2), "cm"),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +
  scale_color_viridis_c(option = "plasma", direction = -1)