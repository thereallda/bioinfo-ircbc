#!/usr/bin/env Rscript
# custom function for DE analysis ----
# PCA plot from counts matrix
ggPCA <- function(data, label = colnames(data), palette=NULL) {
    # data: count matrix, can be transformed
    # label: character vector of sample names or labels. Defaults to colnames(data).
    pca <- prcomp(t(data))
    rownames(pca$x) <- label
    pc.var <- round(summary(pca)$importance[2,], 3)
    pca_dat <- as.data.frame(pca$x) %>% mutate(group = str_remove(label, '\\.\\d$'))

    if (is.null(palette)) {
    palette <- paintingr::paint_palette("Spring", length(unique(pca_dat$group)), 'continuous')
    }

    pca_dat %>%
    ggplot(aes(x=PC1, y=PC2)) +
    geom_point(aes(color=group), size=3) +
    ggrepel::geom_text_repel(label=label, max.overlaps = 20) +
    geom_vline(xintercept=0, color='grey80', lty=2) +
    geom_hline(yintercept=0, color='grey80', lty=2) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = 'top', axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    labs(x=paste0('PC1: ', pc.var[1]*100, '%'),
         y=paste0('PC2: ', pc.var[2]*100, '%'))
}
# Correlation heatmap
corHeatmap <- function(data, col.names, row.names) {
  # data: count matrix, can be transformed

  corMat <- cor(data)
  colnames(corMat) <- col.names
  rownames(corMat) <- row.names
  od <- hclust(dist(corMat))$order
  corMat <- corMat[od,od]
  rowDend <-  as.dendrogram(hclust(dist(corMat)))
  corp1 <- ComplexHeatmap::Heatmap(corMat, name='corMat', 
                                  col=colorRampPalette(RColorBrewer::brewer.pal(9,'YlOrRd'))(50),
                                  cluster_columns = F, cluster_rows = rowDend, row_dend_side = 'left',
                                  row_dend_width = unit(2, "cm"),
                                  column_names_rot = 45, row_names_side = 'left', use_raster = T)
  return(corp1)
}
# Wrap-up DESeq2 procedure 
DESeqOne <- function(data, col.data, design, contrast.df, log.fc=1, p.adj=0.05, outdir=NULL) {
  # output directory and prefix (optional)
  if (is.null(outdir)) outdir = 'results/de/'

  dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = col.data,
                              design = design)
  de <- DESeq(dds)
  
  # Noted: DE will be performed as first element vs second element in the vector
  res.ls <- lapply(seq_len(nrow(contrast.df)), function(i) {
        restmp <- results(de, contrast = c('condition',contrast.df[i,1],contrast.df[i,2]), tidy=TRUE) %>% 
            dplyr::rename(GeneID = row)
  })
  names(res.ls) <- paste0('Res_', apply(contrast.df, 1, paste, collapse='_'))

  # save unfiltered DE results
  lapply(seq_along(res.ls), function(x, nm, i) {
    write.csv(x = x[[i]], file = paste0(outdir, nm[[i]], '_DE.csv'), row.names=FALSE)
  }, x=res.ls, nm=names(res.ls))

  # significant DEGs are considered as genes that pass cutoff:
  # |log2 fold change| >= log.fc and p.adj < p.adj
  res.sig.ls <- lapply(seq_along(res.ls), function(x, nm, i) {
    sig <- subset(x[[i]], padj < p.adj & abs(log2FoldChange) >= log.fc)
    # save significant DEGs tables
    write.csv(x = sig, file = paste0(outdir, nm[[i]], '_DE_sig.csv'), row.names=FALSE)
    return(sig)
  }, x=res.ls, nm=names(res.ls))
  names(res.sig.ls) <- names(res.ls)
  de.res <- list("DE.object" = de, "res.ls" = res.ls, "res.sig.ls" = res.sig.ls)
  return(de.res)
}
# volcano plot
ggVolcano <- function(data, title=NULL) {
    # data: DEG results table
    # title: title of plot
  data <- data %>% 
    mutate(stat = if_else(padj < 0.05 & log2FoldChange >= 1, "Up", 
                          if_else(padj < 0.05 & log2FoldChange <= -1, "Down", "NS", missing = "NS"), missing = "NS"))
  
  count.dat <- data %>% dplyr::count(stat) %>% dplyr::mutate(label = paste0(stat, ": ", n))
  
  data %>% 
    ggplot(aes(x=log2FoldChange, y = -log10(padj), color=stat)) +
    geom_point() +
    theme_classic() +
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(color='black')) +
    scale_color_manual(values = c("#4F99B4","#808080","#CBC28D"), labels = count.dat$label) +
    geom_vline(xintercept = c(-1, 1), lty = "dashed", alpha=0.6) +
    geom_hline(yintercept = -log10(0.05), lty = "dashed", alpha=0.6) +
    labs(color='', title = title)
}

# read in meta information
metadata <- read.csv('metadata.csv', header=TRUE, comment.char = "#")
# read in contrast table
contrast_df <- read.csv('contrast.csv', header=TRUE)

# get work directory
wd.curr <- getwd()
# assume run in <project_dir>/src
# set work directory as <project_dir>
setwd(dirname(wd.curr))
# script for standard differential expression analysis based on DESeq2
library(tidyverse)
library(DESeq2)
if (!dir.exists('results/vis')) dir.create('results/vis', recursive = T)
if (!dir.exists('results/de')) dir.create('results/de', recursive = T)

# read in Counts data
counts_df <- read.csv('results/Counts.csv', row.names = 1)
# step 1. filtering low expressed genes
counts_keep <- counts_df[rowSums(counts_df)>10, ]
# step 2. generate QC plots: PCA, correlation, hclust
samples_name <- paste(metadata$condition, metadata$replicate, sep = '.')
# perform vst transformation
vt <- vst(as.matrix(counts_keep))
## PCA
pcaplot <- ggPCA(vt, label = samples_name)
ggsave(plot = pcaplot, filename = 'results/vis/PCA_PC12.pdf', width=10, height=8)
## Correlation heatmap
corp1 <- corHeatmap(vt, col.names = samples_name, row.names = samples_name)
pdf(file='results/vis/CorHeatmap.pdf', bg='white', width=10, height=8)
ComplexHeatmap::draw(corp1)
dev.off()
## Hierarchical clustering 
n.cond <- length(unique(metadata$condition))
colnames(vt) <- samples_name
clust <- hclust(dist(t(vt)))
pdf('results/vis/HClustering.pdf', width=10, height=8)
par(mar = c(8,4.5,2,2)) # tune the margin of plot
as.dendrogram(clust) %>% 
    dendextend::color_branches(k = n.cond) %>% 
    dendextend::color_labels(k = n.cond) %>%
    plot(ylab='Height', main = 'Cluster dendrogram')
dev.off()
# step 3. DE and output de results
rownames(metadata) <- metadata$id
coldatas <- metadata[, 'condition', drop=FALSE]
de.res <- DESeqOne(data = counts_keep, 
                col.data = coldatas, 
                design = ~condition, 
                contrast.df = contrast_df, 
                log.fc=1, p.adj=0.05)

# save DE results as .RData
save(de.res, file = 'results/de/DE.RData')

# step 4. visualization of DEGs
## Volcano plot
de.tab <- de.res$res.ls$Res_WP_L3
vp1 <- ggVolcano(de.tab)
ggsave('results/vis/VolcanoPlot.pdf', vp1, width = 10, height = 8)
## Heatmap
# get significant DEG
de.tab.sig <- de.res$res.sig.ls$Res_WP_L3
sig.degs <- de.tab.sig$GeneID
# get DESeq2 normalized counts
counts_norm <- DESeq2::counts(de.res$DE.object, normalized=TRUE)
# draw DEGs heatmap, scale expression by row
deg_hm1 <- ComplexHeatmap::pheatmap(counts_norm[sig.degs,], scale='row', show_rownames = FALSE)

pdf('results/vis/DEG_Heatmap.pdf', width=5, height=6)
ComplexHeatmap::draw(deg_hm1)
dev.off()

# step 5. GO enrichment
library(clusterProfiler)
library(org.Dm.eg.db)
sig.degs.up <- subset(de.tab.sig, log2FoldChange >= 1)$GeneID
sig.degs.down <- subset(de.tab.sig, log2FoldChange <= -1)$GeneID
ego_up_bp <- enrichGO(sig.degs.up, ont = 'BP', keyType = 'FLYBASE', 
                      universe = rownames(counts_norm), OrgDb = org.Dm.eg.db, 
                      readable = TRUE)
ego_down_bp <- enrichGO(sig.degs.down, ont = 'BP', keyType = 'FLYBASE', 
                      universe = rownames(counts_norm), OrgDb = org.Dm.eg.db,
                      readable = TRUE)
# visualization with barplot
bp1 <- barplot(ego_up_bp) + scale_fill_gradientn(colors=paintingr::paint_palette('Pearlgirl',100,'continuous')) + ggtitle('WP Up-regulated')
bp2 <- barplot(ego_down_bp) + scale_fill_gradientn(colors=paintingr::paint_palette('Pearlgirl',100,'continuous')) + ggtitle('WP Down-regulated')

ggsave('results/vis/GO_BP_WP_UP_DEG.pdf', bp1, width=6, height=4)
ggsave('results/vis/GO_BP_WP_DOWN_DEG.pdf', bp2, width=6, height=4)
