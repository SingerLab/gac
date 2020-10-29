#' plot consensus clustering
#'
#' @param cnr a cnr bundle
#'
#' @param k parameter K to plot
#'
#' @param ... additional argoments passed to Heatmap
#'
#' @import ComplexHeatmap
#'
#' 
#' @export
plot_ccp <- function(cnr, k, ...) {
    
    tmpH = ComplexHeatmap::Heatmap(round(
                               cnr[["ccp"]][[k]][["consensusMatrix"]], digits = 2),
                               cluster_rows = cnr[["ccp"]][[k]][["consensusTree"]],
                               cluster_columns = cnr[["ccp"]][[k]][["consensusTree"]],
                               column_dend_reorder = TRUE, row_dend_reorder = TRUE,
                               column_title = paste("K", "=", k),
                   show_row_names = FALSE, show_column_names = FALSE,
                   ...)
    tmpH
}
