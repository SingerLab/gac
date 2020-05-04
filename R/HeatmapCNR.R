#' HeatmapCNR
#'
#'
#' Wrapper function to the ultra-powerful ComplexHeatmap.  It implements Bray-Curtis
#' dissimilarity as the distance metric for clustering cells, and 'ward.D2' from
#' `hclust`.  Bray-Curtis disimilarity
#' is extensively used in Ecology for clustering comunities.
#'
#' If you prefer to use a different method, you can use the native Heatmap function
#'
#' 
#' @param cnr the CNR bundle
#'
#' @param what wether you want to plot bins or genes
#'
#' @param which.genes  IF you chose genes, you need to specify which ones
#'
#' @param col color function; in this package you have several color palettes,
#'   segCol for segmentation colors: loss and deletions in blue and yellow,
#'   respectively, 2 in white, gains and amplifications in red to gray scale
#' 
#'   lowcol for log2 ratio colors: blue-to-red scale from -2 to +2
#' 
#'   ampCol for amplification : yellow to red scale from 0 to Inf.
#'
#'   fgaCol for Fraction of genome altered: from 0 e.g.. light yellow to 100 % in blue
#'
#'   ploidyCol for ploidy : gradient of purples starting at 1.5
#'
#' @param cluster_rows logical to cluster or not rows (chromosome positions), default FALSE
#' will plot chromosomes in order
#'
#' @param show_row_dend weather to show the row dendrogram, default is FALSE,
#' 
#' @param ... additional parameters from Heatmap
#' For custumizing your heatmap, please visit the ComplexHeatmap documentation:
#'
#' https://jokergoo.github.io/ComplexHeatmap-reference/book/
#'
#'
#'
#' HeatmapCNR(cnr, col = segCol)
#'
#' HeatmapCNR(dna, col = lowCol)
#' 
#' @export
HeatmapCNR <- function(cnr, what = "X", which.genes = NULL,
                       col = segCol, cluster_rows = FALSE,
                       show_row_dend = FALSE, ...) {

    if(what == "X") {

        use <- cnr$X

    } else {

        if(what == "genes" & !is.null(which.genes)) {

            use <- cnr[["genes"]][, which.genes]

        } 

    }
    
    if(what == "X") {

        chrAnno <- rowAnnotation(chr = factor(cnr$chromInfo$chrom))
        
        Hmap <- Heatmap(use, name = "X", clustering_distance_columns = function(X) vegan::vegdist(X, method = "bray"),
                        col = col,
                        cluster_rows = cluster_rows,
                        left_annotation = chrAnno,
                        clustering_method_columns = "ward.D2",
                        show_row_dend = FALSE, ...)
    } else {
        
        if(what == "genes" & all(which.genes %in% colnames(cnr$genes))) {

            Hmap <- Heatmap(use, name = "genes", clustering_distance_rows = function(X) vegan::vegdist(X, method = "bray"),
                            col = col,
                            cluster_rows = cluster_rows,
                            clustering_method_rows = "ward.D2",
                             ...)

        }

    }

    return(Hmap)

} # HeatmapCNR
