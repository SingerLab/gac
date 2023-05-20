#' plot consensus clustering
#'
#' @param cnr a cnr bundle
#'
#' @param k parameter K to plot
#'
#' @param col color map, defaults to Color Brewer `YlGnBu`
#' 
#' @param ... additional argoments passed to Heatmap
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' 
#' @export
plot_ccp <- function(cnr, k = NULL, col = NULL, ...) {

    if(is.null(col)) {
        cluster.block.colors <- circlize::colorRamp2(
                                              breaks = c(0, 0.5, 1),
                                              col = c("#FFFFD9", "#41B6C4", "#081D58"))
    } else {
        cluster.block.colors <- col
    }

    if(is.null(k)) {
        k <- cnr$kCC
    }
    
    tmpH = ComplexHeatmap::Heatmap(round(
                               cnr[["ccp"]][[k]][["consensusMatrix"]], digits = 2),
                               cluster_rows = cnr[["ccp"]][[k]][["consensusTree"]],
                               cluster_columns = cnr[["ccp"]][[k]][["consensusTree"]],
                               column_title = paste("K", "=", k),
                               col = cluster.block.colors,
                               ...)
    tmpH
}
