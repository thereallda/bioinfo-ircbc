# helper function 
#' PCA plot from counts matrix
#'
#' @param object A count matrix.
#' @param group Vector of sample groups.
#' @param label Vector of sample names or labels.
#' @param vst.norm if TRUE perform vst transformation.
#' @param palette The color palette for different groups.
#'
#' @return ggplot2 object
#' @export
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom DESeq2 vst
#' @importFrom ggrepel geom_text_repel
#' @importFrom stats prcomp
#' @importFrom paintingr paint_palette
ggPCA <- function(object, group, label=NULL, vst.norm=FALSE, palette=NULL) {
  if (vst.norm) {
    counts_norm <- DESeq2::vst(as.matrix(object))
  } else {
    counts_norm <- object
  }
  
  pca <- prcomp(t(counts_norm))
  pc.var <- round(summary(pca)$importance[2,], 3)
  pca_dat <- as.data.frame(pca$x) %>%
    mutate(group = group)
  
  if (is.null(palette)) {
    palette <- paintingr::paint_palette("Spring", length(unique(pca_dat$group)), 'continuous')
  }
  
  p <- pca_dat %>%
    ggplot(aes(x=PC1, y=PC2)) +
    geom_point(aes(color=group), size=3) +
    geom_vline(xintercept=0, color='grey80', lty=2) +
    geom_hline(yintercept=0, color='grey80', lty=2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'top',
          axis.text = element_text(color='black')) +
    scale_color_manual(values = palette) +
    labs(x=paste0('PC1: ', pc.var[1]*100, '%'),
         y=paste0('PC2: ', pc.var[2]*100, '%'))
  
  if (is.null(label)) {
    label <- colnames(object)
  }
  p <- p + ggrepel::geom_text_repel(label=label, max.overlaps = 20)
  return(p)
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
