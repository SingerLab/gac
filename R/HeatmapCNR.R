#' HeatmapCNR
#'
#' Wrapper function to the ultra-powerful ComplexHeatmap.  It implements Bray-Curtis
#' dissimilarity as the distance metric for clustering cells, and 'ward.D2' from
#' `hclust`.  Bray-Curtis disimilarity
#' is extensively used in Ecology for clustering comunities.
#'
#' If you prefer to use a different method, you can use the native Heatmap
#' function.  E.g. if you prefer genomic bins to cluster, and show a dendrogram
#' when plotting X. By default bins are kept in chromosome order, however, when
#' plotting genes, rows are cells, and these are clustered.
#'
#'
#' @param cnr the CNR bundle
#'
#' @param what wether you want to plot bins or genes
#'
#' @param which.genes  IF you chose genes, you need to specify which ones
#'
#' @param col optional color map, if NULL colors are dependent on the sample ploidy:
#'   For base.ploidy = 2; values are 0 = yellow, 1 = blue, 2 = white, 3-10 reds, >10 greyscale
#'   For base.ploidy = 4; values are a <4 = dark to light blues, 4 = yellow, >4 = light to dark reds 
#'   if cnr$bulk = TRUE, ratios are blue, white, and red.
#' 
#' @param base.ploidy base ploidy, values must be 2 or 4.  Defaults to 2.  Ignored if setting color pallete
#'
#' @param show_row_dend weather to show the row dendrogram, default is FALSE,
#' 
#' @param ... additional parameters from Heatmap
#'
#' @import ComplexHeatmap
#'
#' @return
#' Returns a simple ComplexHeatmap plot clustered using Bray-Disimilarity with
#' vegan::vegdist, and sorted by chromosome location.  
#' 
#' For custumizing your heatmap, please visit the ComplexHeatmap documentation:
#'
#' https://jokergoo.github.io/ComplexHeatmap-reference/book/
#' 
#' @examples
#'
#' ## load example
#' data(cnr)
#'
#' HeatmapCNR(cnr)
#'
#' HeatmapCNR(cnr, what = "genes", which.genes = c("CDK4", "MDM2"))
#'
#' @import ComplexHeatmap
#' @importFrom vegan vegdist 
#' @importFrom circlize colorRamp2
#' 
#' @export
HeatmapCNR <- function(cnr, what = "X", which.genes = NULL,
                       col = NULL, base.ploidy = 2,
                       show_row_dend = FALSE, ...) {

    if(what == "X") {

        use <- as.matrix(cnr$X)

    } else {

        if(what == "genes" & !is.null(which.genes)) {

            use <- cnr[["genes"]][, which.genes]

        }

    }

    if(is.null(col)) {
        
        if(base.ploidy == 2) {
            at <- c(0, 1, 2, 3, 4, 5, 10, 30)
            segment.colors <- circlize::colorRamp2(
                                 breaks = at,
                                 colors = c("#F2D200", "#318CE7",
                                            "#FFFFFF", "#FFA89F",
                                            "#FF523F", "#D40000",
                                            "#8F8F8F", "#141414"))
        } else {
            at = c(0, 1, 2, 3, 4, 8, 16, 32, 64)
            segment.colors <- circlize::colorRamp2(
                                            breaks = at,
                                            colors = c("#313695", "#4575B4", "#74ADD1",
                                                       "#ABD9E9", "#FFFFBF", "#FDAE61",
                                                       "#F46D43", "#D73027", "#A50026"))
        }
        
        if(is.null(col)) {
            bulk.colors <- circlize::colorRamp2(breaks = c(-1,-0.16, 0.16, 1),
                                                colors = c("#2166AC", "#FFFFFF",
                                                           "#FFFFFF", "#D40000"))
        } else {
            bulk.colors <- col
        }
    } else {
        segment.colors <- col
        bulk.colors <- col
    }

    if(what == "X") {

        chrAnnoLeft <- create_chromosome_annotation(cnr)
        
        if(cnr[["bulk"]]) {
            
            Hmap <- ComplexHeatmap::Heatmap(use, name = "X", clustering_distance_columns = function(X) vegan::vegdist(2^X, method = "bray", na.rm = TRUE),
                            cluster_rows = FALSE,
                            left_annotation = chrAnnoLeft,
                            clustering_method_columns = "ward.D2",
                            col = segment.colors,
                            ...)
            
        } else {

            Hmap <- ComplexHeatmap::Heatmap(use, name = "X", clustering_distance_columns = function(X) vegan::vegdist(X, method = "bray", na.rm = TRUE),
                            cluster_rows = FALSE,
                            left_annotation = chrAnnoLeft,
                            clustering_method_columns = "ward.D2",
                            col = segment.colors,
                            ...)
        }
        
    } else {
        
        if(what == "genes" & all(which.genes %in% colnames(cnr$genes))) {
            
            if(cnr[["bulk"]]) {

                Hmap <- ComplexHeatmap::Heatmap(use, name = "genes", clustering_distance_rows = function(genes) vegan::vegdist(2^genes, method = "bray"),
                                clustering_method_rows = "ward.D2",
                                col = bulk.colors,
                                ...)
                
            } else {
                
                Hmap <- ComplexHeatmap::Heatmap(use, name = "genes", clustering_distance_rows = function(genes) vegan::vegdist(genes, method = "bray"),
                                clustering_method_rows = "ward.D2",
                                col = segment.colors,
                                ...)
            }
          
        }
        
    }
    
    return(Hmap)
    
} # HeatmapCNR
