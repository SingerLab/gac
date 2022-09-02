#' plot genome wide correlations
#'
#' @param cnr a cnr bundle
#'
#' @param ... additional parameters passed to \code{Heatmap}
#' 
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom grid gpar unit
#' @importFrom circlize colorRamp2
#' @importFrom assertthat assert_that
#' 
#' @export
plot_cn_correlations  <- function(cnr, ...) {

    assertthat::assert_that(!is.null(cnr[["gw_corr"]]))
    
    chrAL <- create_chromosome_annotation(cnr, side = "left")
    chrAT <- create_chromosome_annotation(cnr, side = "top")

    pwColors <- circlize::colorRamp2(breaks = c(0, 0.1, 0.5, .95),
                                     colors = c("#FFFFFF", "#EFF3FF",
                                                "#6BAED6", "#08519C"))
    
    gwH <- Heatmap(cnr[["gw_corr"]],
                   col = pwColors,
                   border = "#FFFFFF",
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_split = cnr$chromInfo$bin.chrom,
                   column_split = cnr$chromInfo$bin.chrom,
                   cluster_column_slices = FALSE,
                   cluster_row_slices = FALSE,
                   row_gap = grid::unit(0, "mm"),
                   column_gap = grid::unit(0, "mm"),
                   top_annotation = chrAT,
                   left_annotation = chrAL,
                   show_heatmap_legend = FALSE, ...
                   )
} # end plot_cn_correlations
