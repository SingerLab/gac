#' Mark genes
#'
#' Produces a named list with gene name and bin.id to use in ComplexHeatmap::anno_mark
#'
#' @param cnr  the cnr bundle.  Only gene index will be used
#'
#' @param gene.list a gene list
#'
#'
#' @return
#' returns a vector of named gene.list (gg) w/their rownames/bin.id
#' 
#' @examples
#'
#' data(cnr)
#'
#' aa <- mark.genes(cnr$gene.index, gene.list = c("CDK4", "MDM2"))
#'
#' geneAnno <- rowAnnotation(genes = anno_mark(at = aa, labels = names(aa)))
#'
#' 
#' 
#' @export
mark.genes <- function(cnr, gene.list) {
    ## genes to mark in Heatmap
    gg <- cnr$gene.index[gene.list, "bin.id"]
    names(gg) <- gene.list
    return(gg)
} ## mark.genes
