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
#' summary_cnr(cnr)
#'
#' @importFrom dplyr tibble
#' 
#' @export
summary_cnr<- function(cnr, detailed = FALSE, display = FALSE) {

    message("Data Type:  ")
    print(cnr_data_type(cnr, detailed = detailed))
    message()
    message("number of cells/samples:")
    print(ncells(cnr, detailed = detailed, display = display))
    message("number of phenotypes:")
    print(npheno(cnr, detailed = detailed, display = display))
    message("number of bins:")
    print(nbins(cnr, detailed = detailed, display = display))
    message("number of genes:")
    print(ngenes(cnr, detailed = detailed, display = display))
    message()
    message("Tasks performed:")
    message("*   Bray-Curtis Diss.: ", !is.null(cnr$cdb))
    message("*   Phylo Tree: ", !is.null(cnr$phylo))
    message("*   Consensus Clustering: ", !is.null(cnr$ccp))
    message("*   K Spectral: ", !is.null(cnr$kStats))
    if(!is.null(cnr$kStats)) {
        cnr$optK
    }
    message("*   Alteration Frequencies: ", any(grepl("AmpFQ", names(cnr$gene.index))))
    ## message("GISTIC2: ", !is.null(cnr$gistic))
    message("*   VDJ Genotyping: ", !is.null(cnr$vdj.cells))
    message("*   Cluster Heterogeneity: ", !is.null(cnr$cluster_heterogeneity))
    message("*   Cluster Profiles: ", !is.null(cnr$DDRC.df))
    
    return(invisible(NULL))
    
}

#' data type in CNR
#'
#' @param cnr a cnr bundle
#' 
#' @param detailed logical, if TRUE print data type with data range, default FALSE
#'
#' @examples
#' data(cnr)
#' cnr_data_type(cnr)
#'
#' cnr_data_type(cnr, detailed = TRUE)
#' 
#' @export
cnr_data_type <- function(cnr, detailed = FALSE) {

    out <- ifelse(cnr$bulk, "log2 Ratio", "Integer Copy Number")

    if(detailed) {
        out = list("data.type" = out, "range" = range(cnr$X))
    }
    
    return(out)
}


#' number of cells in cnr
#'
#' @param cnr a cnr bundle
#'
#' @param detailed logical, detailed summary e.g print total cells, n vdj cells, & cells by cell type
#'
#' @param detailed.column character, name of column to partition cells, default is NULL
#'
#' @param display display cell ID's
#'
#' @importFrom assertthat assert_that
#'
#' @examples
#' data(cnr)
#'
#' ncells(cnr)
#'
#' ncells(cnr, detailed = TRUE)
#'
#' ncells(cnr, detailed = TRUE,
#'        detailed.column = "category1")
#' 
#' @export
ncells <- function(cnr, detailed = FALSE,
                   detailed.column = NULL,
                   display = FALSE) {

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

    if(detailed && !is.null(detailed.column)) {
        
        assertthat::assert_that(detailed.column %in% names(cnr$Y))
        
        ntc <- table(cnr$Y[, detailed.column])
        
        nxc = c("n.total.cells" = nyc,
                "n.vdj.cells" = nvc,
                ntc)
    } else {
        if(detailed & is.null(detailed.column)) {
            nxc <- c("n.total.cells" = nyc,
                     "n.vdj.cells" = nvc)
        }
    }
    if(display) print(cnr$cells)
    
    return(nxc)
    
} ## end of ncells

#' number of phenotypes
#' @param cnr a cnr bundle
#'
#' @param detailed detailed summary
#'
#' @param display logical display example
#'
#' @examples
#' data(cnr)
#'
#' npheno(cnr)
#'
#' npheno(cnr, detailed = TRUE)
#'
#'
#' @export
npheno <- function(cnr, detailed = FALSE,
                   display = FALSE) {
    np <- ncol(cnr$Y)
    ## message("number of phenotypes:", np)

    if(detailed) {
        message()
        message("Detailed Phenotype Summary:")
        print(summary(cnr$Y))
        message()

    }

    if(display) cnr$Y

    return(np)
}

#' number of bins
#' @param cnr a cnr bundle
#'
#' @param detailed detailed summary
#'
#' @param display  logical, display chromInfo
#'
#' @examples
#'
#' data(cnr)
#' 
#' nbins(cnr)
#'
#' nbins(cnr, detailed = TRUE)
#'
#' @importFrom dplyr tibble
#' 
#' @export
nbins <- function(cnr, detailed = FALSE, display = FALSE) {

    nb <- nrow(cnr$X)
    
    ## message("number of bins: ", nb)

    if(detailed) {
        message()
        message("Detailed bin summary:")
        print(table(cnr$chromInfo$bin.chrom))
        message()
    }

    if(display) dplyr::tibble(cnr$chromInfo)

    return(nb)
}
    
    
#' number of genes in gene.index
#' @param cnr a cnr bundle
#'
#' @param detailed logical, detailed summary, if TRUE prin number of genes by chromosome,
#'   default FALSE
#'
#' @param display  logical, display gene.index
#'
#' @importFrom dplyr tibble
#'
#' @examples
#' data(cnr)
#'
#' ngenes(cnr)
#'
#' ngenes(cnr, detailed = TRUE)
#' 
#' 
#' @export
ngenes <- function(cnr, detailed = FALSE, display = FALSE) {

    ng <- nrow(cnr$gene.index)
    ngc <- ncol(cnr$genes)
    

    if(detailed) {
        message()
        message("number of genes in index: ", ng)
        message("number of genes in use: ", ngc)
        message()
        message("Detailed Summary of Genes:")
        message("number of genes per chromosome:")
        print(table(cnr$gene.index$chrom))
        message()
        message("summary of genes per bin:")
        print(summary(as.numeric(table(cnr$gene.index$bin.id))))
        message()
    }
    
    if(display) dplyr::tibble(cnr$gene.index)

    if(ng != ngc) {
    out <- c("n.genes.data" = ngc,
             "n.genes.index" = ng)
    } else {
        out = ngc
    }
    
    return(out)
}
    
    
