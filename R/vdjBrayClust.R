#' VDJ cells bray clustering 
#'
#' @param cnr a cnr bundle, previously genotyped for VDJ
#'
#' @param vdj.genes list of VDJ genes for clustering. All must be in the
#'    gene.index `hgnc.symbol` column
#' 
#' @importFrom ape as.phylo
#' @importFrom assertthat assert_that
#' 
#' @export
vdjBrayClust <- function(cnr, vdj.genes = NULL) {

    assertthat::assert_that(!is.null(vdj.genes))
    assertthat::assert_that(all(vdj.genes %in% cnr$gene.index$hgnc.symbol))
    
    dd <- vegan::vegdist(cnr$genes[names(cnr$vdj.cells), vdj.genes])
    
    hc <- hclust(dd, method = "ward.D2")

    cnr[["vdjDist"]] <- dd
    cnr[["vdjClust"]] <- hc
    cnr[["vdjPhylo"]] <- ape::as.phylo(hc)
    
    return(cnr)
}

#' plot VDJ clustering
#'
#' @param cnr a cnr bundle with clustered vdj cells
#'
#' @param type type of tree, from \code{\link[ape]{plot.phylo}}, default is "fan", 
#' 
#' @param ... additional parameters to \code{\link[ape]{plot.phylo}}
#'
#' @importFrom assertthat assert_that
#' 
#' @export
vdjPlotClust <- function(cnr, type = "fan",  ...) {
    assertthat::assert_that(!is.null(cnr[["vdjPhylo"]]))
    plot(cnr[["vdjPhylo"]], type = type,  ...)
}
