#' plot consensus clustering
#'
#' @param cnr a cnr bundle
#'
#' @param k parameter K to plot
#'
#' @param ... additional argoments passed to Heatmap
#'
#' @importFrom ComplexHeatmap Heatmap
#' 
#' @export
plot_ccp <- function(cnr, k, ...) {
    
    tmpH = ComplexHeatmap::Heatmap(round(
                               cnr[["ccp"]][[k]][["consensusMatrix"]], digits = 2),
                               cluster_rows = cnr[["ccp"]][[k]][["consensusTree"]],
                               cluster_columns = cnr[["ccp"]][[k]][["consensusTree"]],
                               column_title = paste("K", "=", k),
                   ...)
    tmpH
}
