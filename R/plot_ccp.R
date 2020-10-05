#' plot consensus clustering
#'
#' @param cnr a cnr bundle
#'
#' @param k parameter K to plot
#'
#' @param col colour function for the heatmap
#' 
plot_ccp <- function(cnr, k, col = ccClustCol, ...) {
    assertthat::assert_that("cluster" %in% names(cnr$Y))
    set.seed(81)

    clColors <- sample(colors(), length(unique(cnr$Y$cluster)))
    names(clColors) <- unique(cnr$Y$cluster)
    ccColors <- sample(colors(), k)
    names(ccColors) <- k
    
    tmpH = Heatmap(round(cnr[["ccp"]][[k]][["consensusMatrix"]], digits = 2),
                   cluster_rows = cnr[["ccp"]][[k]][["consensusTree"]],
                   cluster_columns = cnr[["ccp"]][[k]][["consensusTree"]],
                   column_dend_reorder = TRUE, row_dend_reorder = TRUE,
                   column_title = paste("K", "=", k),
                   top_annotation =
                       HeatmapAnnotation(cc = cnr[["ccp"]][[k]][["consensusClass"]],
                                         BrayC = cnr$Y$cluster,
                                         col = list(cc = ccColors,
                                                    BrayC = clColors),
                                         show_legend = FALSE),
                   left_annotation =
                       rowAnnotation(cc = cnr[["ccp"]][[k]][["consensusClass"]],
                                     BrayC = cnr$Y$cluster,
                                     col = list(cc = ccColors,
                                                BrayC = clColors),
                                     show_legend = FALSE),
                   show_row_names = FALSE, show_column_names = FALSE)
}
