# Figure 3 (source)

{
  library(Seurat)
  library(Matrix)
  library(R.utils)
  library(Matrix.utils)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(GenomicAlignments)
  library(tidyr)
  library(glue)
  library(Rsamtools)
  library(ggExtra)
  library(ggrepel)
  library(rtracklayer)
  library(GenomicRanges)
  options(scipen = 999)
}

{
  ReadMat <- function(bpath) {
    mat <- readMM(paste0(bpath, "/matrix.mtx"))
    rownames(mat) <- read.table(paste0(bpath, "barcodes.tsv"))$V1
    colnames(mat) <- read.table(paste0(bpath, "features.tsv"))$V1
    mat <- t(mat)
    CreateSeuratObject(mat)
  }
  
  GetObject <- function(kind) {
    if (kind == "sr.pbmc.sens4") {
      readRDS("SR.PBMC.S3.rds")
    } else if (kind == "sr.pbmc.sens3") {
      readRDS("SR.PBMC.S3.rds")
    } else if (kind == "sr.pbmc.10x") {
      obj <- CreateSeuratObject(Read10X("filtered_feature_bc_matrix"))
      NormalizeData(obj)
    } else if (kind == "sr.pbmc.alevin") {
      apath <- "alevin/"
      mat <- ReadMtx(
        mtx = paste0(apath, "quants_mat.mtx.gz"),
        cells = paste0(apath, "quants_mat_cols.txt"), 
        features = paste0(apath, "quants_mat_rows.txt"), 
        feature.column = 1
      )
      obj <- CreateSeuratObject(t(mat))
      NormalizeData(obj)
    } else if (kind == "sr.pbmc.smartseq") {
      apath <- "HCA.UMIcounts.PBMC.txt"
      mat <- read.table(apath)
      obj <- CreateSeuratObject(mat)
      NormalizeData(obj)
    } else if (kind == "lr.pbmc.masiso") {
      tpath <- "iso"
      tmat <- Read10X(tpath, gene.column = 1)
      
      obj <- NormalizeData(CreateSeuratObject(tmat))
      txp.tpm <- GetAssayData(obj, assay = "RNA", layer = "data")
      txp.tpm <- rowSums(txp.tpm) * 1000000 / sum(rowSums(txp.tpm))
      names(txp.tpm) <- unlist(lapply(names(txp.tpm), function(x) {
        strsplit(x, ".", fixed=T)[[1]][[1]]
      }))
      
      gpath <- "genes"
      gmat <- Read10X(gpath, gene.column = 1)
      cts <- list (
        "gene" = NormalizeData(CreateSeuratObject(gmat)),
        "txp" = txp.tpm
      )
    } else if (kind == "lr.pbmc.95") {
      ReadMat("v3/")
    } else if (kind == "bulk.pbmc") {
      read.table("pbmc.tx.bulk.quant.sf", 
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
    } else {
      print("Can't recognize the Kind")
    }
    
  } #end GetObject function
  
  GroupAtGene <- function(txp.counts) {
    t2g <- read.table("t2gnames.txt", 
                      row.names = 1)
    
    genes <- t2g[names(txp.counts), ]
    tapply(txp.counts, genes, sum)
  }
  
  Getcpm <- function(obj){
    gene.cpm <- rowSums(GetAssayData(obj, assay = "RNA", layer = "data"))
    gene.cpm * 1000000 / sum(gene.cpm)
  }
  
  Gettpm <- function(obj, kind="lr") {
    if (kind=="lr") {
      txp.tpm <- GetAssayData(obj, assay = "RNA", layer = "data")
      txp.tpm <- rowSums(txp.tpm) * 1000000 / sum(rowSums(txp.tpm))
      gene.cpm <- GroupAtGene(txp.tpm)
      cts <- list (
        "tpm" = txp.tpm,
        "cpm" = gene.cpm
      )
    } else if (kind == "bulk") {
      tpm <- setNames(obj$TPM, rownames(obj))
      cts <- list (
        "tpm" = tpm,
        "cpm" = GroupAtGene(tpm)
      )
    }
    cts
  }
  
  FillGeneMetadata <- function(df) {
    t2type <- read.table("t2g2type.txt",
                         row.names = 1)
    t2type[t2type$V3 %in% c("ARMCX5-GPRASP2", "CYB561D2", "GOLGA8M"), "V2"] <- "lncRNA"
    t2type <- unique(t2type)
    rownames(t2type) <- t2type$V3
    
    df$type <- t2type[rownames(df), "V2"]
    df
  }
  
  FillTxpMetadata <- function(df) {
    t2g <- read.table("t2gnames.txt", 
                      row.names = 1)
    t2type <- read.table("t2g2type.txt",
                         row.names = 1)
    t2l <- read.table("tlens.txt", 
                      row.names = 1)
    
    df$gene <- t2g[rownames(df),]
    df$type <- t2type[rownames(df),"V2"]
    df$tlen <- t2l[rownames(df),]
    df
  }
}