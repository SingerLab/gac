#' HeatmapCNR
#'
#'
#' Wrapper function to the ultra-powerful ComplexHeatmap.  It implements Bray-Curtis
#' dissimilarity as the distance metric for clustering cells, and 'ward.D2' from
#' `hclust`.  Bray-Curtis disimilarity
#' is extensively used in Ecology for clustering comunities.
#'
#' If you prefer to use a different method, you can use the native Heatmap function.
#' E.g. if you prefer genomic bins to cluster, and show a dendrogram when plotting X.
#' By default bins are kept in chromosome order, however, when plotting genes, rows
#' are cells, and these are clustered.
#'
#' For convenience, GAC comes with several color palettes as part of the data. To use
#'  any of these run `data(segCol)` for example.  There are the different pallents and
#'  suggested use
#' 
#'   * segCol for segmentation colors: loss and deletions in blue and yellow,
#'   respectively, 2 in white, gains and amplifications in red to gray scale
#' 
#'   * lowCol for log2 ratio colors: blue-to-red scale from -2 to +2
#' 
#'   * ampCol for amplification : yellow to red scale from 0 to Inf.
#'
#'   * fgaCol for Fraction of genome altered: from 0 to 1 light yellow to 1 in blue
#'
#'   * ploidyCol for ploidy : gradient of purples starting at 1.5
#'
#' @param cnr the CNR bundle
#'
#' @param what wether you want to plot bins or genes
#'
#' @param which.genes  IF you chose genes, you need to specify which ones
#'
#' @param show_row_dend weather to show the row dendrogram, default is FALSE,
#' 
#' @param ... additional parameters from Heatmap
#' 
#'
#' 
#'
#' @return
#' Returns a simple ComplexHeatmap plot clustered using Bray-Disimilarity with
#' vegan::vegdist, and sorted by chromosome location.  Several color scales are
#' provided in the package, see segCol, lowCol, ampCol, fgaCol, and ploidyCol.
#' These are useful for setting additional rowAnnotation and HeatmapAnnotations.
#' A default chromosome color palette was left out as different organisms have
#' different chromosome numbers and naming conventions.
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
#' ## load color map for segment colors
#' data(segCol)
#'
#' 
#' HeatmapCNR(cnr, col = segCol)
#'
#' HeatmapCNR(cnr, col = segCol, what = "genes", which.genes = c("CDK4", "MDM2"))
#'
#' @import ComplexHeatmap
#' 
#' @export
HeatmapCNR <- function(cnr, what = "X", which.genes = NULL,
                       show_row_dend = FALSE, ...) {

    if(what == "X") {

        use <- as.matrix(cnr$X)

    } else {

        if(what == "genes" & !is.null(which.genes)) {

            use <- cnr[["genes"]][, which.genes]

        }

    }
    
    if(what == "X") {
        cf <- factor(cnr$chromInfo$chrom)
        grs <-  c("#404040", "#BABABA")
        rp <- ceiling(length(unique(cf))/2)
        chl <- rep(grs, rp)
        chl <- chl[1:length(unique(cf))]
        names(chl) <- unique(cf)
        
        chrAnno <- rowAnnotation(chr = cf, col = list(chr = chl))
        if(cnr[["bulk"]]) {
            
            Hmap <- Heatmap(use, name = "X", clustering_distance_columns = function(X) vegan::vegdist(2^X, method = "bray", na.rm = TRUE),
                            cluster_rows = FALSE,
                            left_annotation = chrAnno,
                            clustering_method_columns = "ward.D2",
                            ...)
            
        } else {

            Hmap <- Heatmap(use, name = "X", clustering_distance_columns = function(X) vegan::vegdist(X, method = "bray", na.rm = TRUE),
                            cluster_rows = FALSE,
                            left_annotation = chrAnno,
                            clustering_method_columns = "ward.D2",
                            ...)
        }
        
    } else {
        
        if(what == "genes" & all(which.genes %in% colnames(cnr$genes))) {
            
            if(cnr[["bulk"]]) {
                Hmap <- Heatmap(use, name = "genes", clustering_distance_rows = function(genes) vegan::vegdist(2^genes, method = "bray"),
                                clustering_method_rows = "ward.D2",
                                ...)
                
            } else {
                
                Hmap <- Heatmap(use, name = "genes", clustering_distance_rows = function(genes) vegan::vegdist(genes, method = "bray"),
                                clustering_method_rows = "ward.D2",
                                ...)
            }
          
        }
        
    }
    
    return(Hmap)
    
} # HeatmapCNR
