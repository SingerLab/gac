#' Read DepMap D2_combined file
#'
#' This function is file specific to the Depmap project D2_combined file
#'
#' It will read, and process the row.names to preserve the hgnc.symbols and entrezID 
#'
#' @param depmap Genetic Dependency Combined RNAi; Depmap D2_conbined object  (see: https://depmap.org/portal/)
#'
#' @param keep which cell lines you want to keep
#'
#' @return
#' Returns the depscores for a set of cells of interest
#' 
#' @examples
#'
#' \dontrun{
#'
#' d2_map <- read.csv("inst/extdata/D2_combined_gene_dep_scores.csv", header = TRUE, row.names = 1)
#'
#' ## select cells to use
#' keep <- grep("BONE", names(depmap))
#' 
#' d2_bone <- DepScores(depmap = d2_map, keep = keep)
#'
#'
#' depMeans <- rowMeans(d2_bone[,keep])
#'
#' }
#' 
#' @export
DepScores <- function(depmap, keep) {
    depScores <- depScores[, keep]
    
    attr(depScores, "hgnc.symbol") <- unlist(gsub("(.*) \\(([0-9]+)\\)", "\\1", rownames(depScores)))
    attr(depScores, "entrezID") <- unlist(gsub("(.*) \\(([0-9]+)\\)", "\\2", rownames(depScores)))
    rownames(depScores) <- attr(depScores, "hgnc.symbol")
    
    depScores
    
} ## DepScores
