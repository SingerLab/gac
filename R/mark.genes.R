#' Mark genes
#'
#' Produces a named list with gene name and bin.id to use in ComplexHeatmap::anno_mark
#'
#' @param cnr  the cnr bundle
#'
#' @param gene.list the list of genes you wish to mark
#'
#'
#' @return
#' returns a vector of named gene.list (gg) w/their rownames/bin.id
#' 
#' @examples
#'
#' data(cnr)
#'
#' aa <- mark.genes(cnr, gene.list = c("CDK4", "MDM2"))
#'
#' geneAnno <- rowAnnotation(genes = anno_mark(at = aa, labels = names(aa)))
#'
#' 
#' 
#' @export
mark.genes <- function(cnr, gene.list) {
    ## genes to mark in Heatmap

    ggL <- cnr$gene.index$hgnc.symbol %in% gene.list

    na <- setdiff(gene.list, cnr$gene.index$hgnc.symbol)

    if(length(na) > 0) {
        warning(na)
    }
    
    gg <- cnr$gene.index[ggL, "bin.id"]
    names(gg) <- cnr$gene.index[ggL, "hgnc.symbol"]
    return(gg)
    
} ## mark.genes
