#' summary of cnr bundle
#'
#' @param cnr a cnr bundle
#'
#' @param detailed logical, to provide additional detail on the bundle. Default (FALSE)
#'
#' @param display  logical, weather to display data and other results. Default (FALSE)
#'
#' @return
#' A summary of the contents in the cnr
#'
#' @examples
#'
#' data(cnr)
#'
#' summaryCNR(cnr)
#'
#' @importFrom dplyr tibble
#' 
#' @export
summaryCNR <- function(cnr, detailed = FALSE, display = FALSE) {

    ncells(cnr, detailed = detailed, display = FALSE)
    npheno(cnr, detailed = detailed, display = FALSE)
    nbins(cnr, detailed = detailed, display = FALSE)
    ngenes(cnr, detailed = detailed, display = FALSE)
    message()
    message("Bray-Curtis Diss.: ", !is.null(cnr$cdb))
    message("Phylo Tree: ", !is.null(cnr$phylo))
    message("Consensus Clustering: ", !is.null(cnr$ccp))
    message("K Spectral: ", !is.null(cnr$kStats))
    if(!is.null(cnr$kStats)) {
        cnr$optK
    }
    message("Alteration Frequencies: ", any(grepl("AmpFQ", names(cnr$gene.index))))
    message("GISTIC2: ", !is.null(cnr$gistic))
    message("VDJ Genotyping: ", !is.null(cnr$vdj.cells))
    message("Cluster Heterogeneity: ", !is.null(cnr$cluster_heterogeneity))
    message("Cluster Profiles: ", !is.null(cnr$DDRC.df))
    
}

#' number of cells in cnr
#'
#' @param cnr a cnr bundle
#'
#' @param detailed logical, detailed summary
#'
#' @param display display cell ID's
#' 
#' @importFrom assertthat assert_that
#' 
#' @export
ncells <- function(cnr, detailed = FALSE, display = FALSE) {

    nxc <- ncol(cnr$X)
    nyc <- nrow(cnr$Y)
    nqc <- nrow(cnr$qc)

    ncc <- length(cnr$cells)

    assertthat::assert_that(nxc == nyc)
    assertthat::assert_that(nxc == nqc)
    assertthat::assert_that(nxc == ncc)
    
    if(is.null(cnr$vdj.cells)) {
        nvc <- 0
    } else {
        nvc <- length(cnr$vdj.cells)
    }
    
    message("total cells: ", nyc)
    message("vdj cells:   ", nvc)

    if(detailed & "cell.type" %in% names(cnr$Y)) {
        ntc <- table(cnr$Y$cell.type)
        message("detailed cells")
        print(ntc)
        message()
    }

    if(display) print(cnr$cells)

}

#' number of phenotypes
#' @param cnr a cnr bundle
#'
#' @param detailed detailed summary
#'
#' @param display logical display example
#'
#' @importFrom dplyr tibble
#'
#' @export
npheno <- function(cnr, detailed = FALSE, display = FALSE) {
    np <- ncol(cnr$Y)
    message("number of phenotypes:", np)

    if(detailed) {
        message()
        message("Detailed Phenotype Summary:")
        print(summary(cnr$Y))
        message()

    }

    if(display) dplyr::tibble(cnr$Y)

}

#' number of bins
#' @param cnr a cnr bundle
#'
#' @param detailed detailed summary
#'
#' @param display  logical, display chromInfo
#'
#' @importFrom dplyr tibble
#' 
#' @export
nbins <- function(cnr, detailed = FALSE, display = FALSE) {

    nb <- nrow(cnr$X)
    
    message("number of bins: ", nb)

    if(detailed) {
        message()
        message("Detailed bin summary:")
        print(table(cnr$chromInfo$bin.chrom))
        message()
    }

    if(display) dplyr::tibble(cnr$chromInfo)
}
    
    
#' number of genes in gene.index
#' @param cnr a cnr bundle
#'
#' @param detailed detailed summary
#'
#' @param display  logical, display gene.index
#'
#' @importFrom dplyr tibble
#'
#' @export
ngenes <- function(cnr, detailed = FALSE, display = FALSE) {

    ng <- nrow(cnr$gene.index)
    ngc <- ncol(cnr$genes)
    
    message("number of genes in index: ", ng)
    message("number of genes in use: ", ngc)

    if(detailed) {
        message()
        message("Detailed Summary of Genes:")
        message("number of genes per chromosome:")
        print(table(cnr$gene.index$seqnames))
        message()
        message("summary of genes per bin:")
        print(summary(as.numeric(table(cnr$gene.index$bin.id))))
        message()
    }

    if(display) dplyr::tibble(cnr$gene.index)
}
    
    
