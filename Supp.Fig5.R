# Supplementary Figure 5

# Load packages
rm(list = ls())
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
source("source.R")

{
  ## LOAD SR DATA
  sr.pip <- GetObject("sr.pbmc.sens3")
  sr.tenx <- GetObject("sr.pbmc.10x")
  sr <- Getcpm(sr.tenx)
  # sr.alevin <- GetObject("sr.pbmc.alevin")
  # sr <- Getcpm(sr.alevin)
  
  ## LOAD LR DATA
  lr.masiso <- GetObject("lr.pbmc.masiso")
  lr <- GetObject("lr.pbmc.95")
  lr <- subset(lr, cells = Cells(sr.pip))
  lr <- NormalizeData(lr, normalization.method = "LogNormalize")
  lr <- Gettpm(lr, kind = "lr")
  
  ## LOAD BULK DATA
  bulk <- GetObject("bulk.pbmc")
  bulk <- Gettpm(bulk, kind = "bulk")
}

{ ## GENE-LEVEL ANALYSES
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
  
  # Remove duplicate rows
  sdf <- sdf[!duplicated(sdf), ]
}

{  # txp level analyses
  lr.masiso <- GetObject("lr.pbmc.masiso")
  
  keep.txps <- intersect(names(lr$tpm), names(bulk$tpm))
  keep.txps <- intersect(keep.txps, names(lr.masiso$txp))
  length(keep.txps)
  
  tdf <- data.frame(
    "masiso" = 0,
    "long" = lr$tpm[keep.txps],
    "bulk" = bulk$tpm[keep.txps]
  )
  dim(tdf)
  
  keep.txps.masiso <- intersect(keep.txps, names(lr.masiso$txp))
  tdf[keep.txps.masiso, "masiso"] <- lr.masiso$txp[keep.txps.masiso]
  tdf$masiso <- tdf$masiso*1000000/sum(tdf$masiso)
  cor(tdf, method = "spearman") #0.555239
  
  tdf <- FillTxpMetadata(tdf)
  head(tdf)
  
  cor(tdf[tdf$type == "protein_coding",c(1,2)], 
      method = "spearman") #0.517608
  
  # tdf$fc <- log((tdf$bulk + 0.1) / (tdf$long + 0.1))
  
  # tdf$lfc <- log((tdf$long + 0.1) / (tdf$bulk + 0.1))
  # tdf$mfc <- log((tdf$masiso + 0.1) / (tdf$bulk + 0.1))
  
  # tdf <- tdf[order(tdf$fc),]
  tail(tdf, 20)
  
  {
    sdf <- tdf#[tdf$type == "protein_coding",]
    sdf$mfc <- log10((sdf$masiso + 1) / (sdf$bulk + 1))
    sdf$lfc <- log10((sdf$long + 1) / (sdf$bulk + 1))
    
    radius <- 0.25 
    sdf.p1 <- sdf[sqrt(sdf$lfc^2 + sdf$mfc^2) > radius, ]
    
    n_txps_total   <- nrow(sdf)
    n_txps_plotted <- nrow(sdf.p1)
    n_txps_removed <- n_txps_total - n_txps_plotted
    
    cat("UNIQUE transcripts plotted:", n_txps_plotted, "\n")
    cat("UNIQUE transcripts removed by circle:", n_txps_removed, "\n")
    
    p4 <- ggplot(sdf.p1, aes(y=lfc, x=mfc)) +
      geom_point() +
      geom_hex(bins = 100) +
      scale_fill_viridis_c(option = "plasma", trans = "log10") +
      theme_minimal()
    ggMarginal(p4, yparams = list(fill = "#5C1A8A", color = "black"),
               xparams = list(fill = "#D5722A", color = "black")) # Customize x-axis histogram
  }
}
