# Supplementary Figure 3 and Supplementary Figure 4B

source("Supp.Fig3.Source.R")

{##  LOAD DATA
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
  
  # Remove circle around origin
  radius <- 0.25
  sdf.p1 <- sdf[sqrt(sdf$lfc^2 + sdf$fc^2) > radius, ]
  sdf.p2 <- sdf[sqrt(sdf$sfc^2 + sdf$fc^2) > radius, ]
  sdf.p3 <- sdf[sqrt(sdf$lfc^2 + sdf$sfc^2) > radius, ]
  
  genes_plotted_unique <- bind_rows(
    sdf.p1 %>% dplyr::select(gname),
    sdf.p2 %>% dplyr::select(gname),
    sdf.p3 %>% dplyr::select(gname)
  ) %>%
    distinct(gname)
  
  n_genes_plotted_unique <- nrow(genes_plotted_unique)
  
  # Genes removed by the circle in ALL panels
  genes_removed_unique <- sdf %>%
    filter(
      sqrt(lfc^2 + fc^2) <= radius &
        sqrt(sfc^2 + fc^2) <= radius &
        sqrt(lfc^2 + sfc^2) <= radius
    ) %>%
    distinct(gname)
  
  n_genes_removed_unique <- nrow(genes_removed_unique)
  cat("Unique genes plotted:", n_genes_plotted_unique, "\n")
  cat("Unique genes removed:", n_genes_removed_unique, "\n")
  
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
}



{  # txp level analyses (SUPPLEMENTARY FIGURE 4B)
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
    
    cat("Unique transcripts plotted:", n_txps_plotted, "\n")
    cat("Unique transcripts removed:", n_txps_removed, "\n")
    
    p4 <- ggplot(sdf.p1, aes(y=lfc, x=mfc)) +
      geom_point() +
      geom_hex(bins = 100) +
      scale_fill_viridis_c(option = "plasma", trans = "log10") +
      theme_minimal()
    ggMarginal(p4, yparams = list(fill = "#5C1A8A", color = "black"),
               xparams = list(fill = "#D5722A", color = "black")) # Customize x-axis histogram
  }
}