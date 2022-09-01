#' plot pathview of specific clones
#'
#' this function is a wrapper to plot CNA detected in clones onto KEGG pathways.
#' Requires Bioconductor/Pathview to be installed
#' 
#' @param cnr a cnr bundle after running \code{\link{cluster_heterogeneity}}
#'
#' @param clones a list of clones (names from DDRC.g)
#'
#' @param pathway.id a pathview pathway ID use \code{link[pathview]{paths.hsa}}.
#'  Default "hsa05200"
#'
#' @param species use species from pathview. Default "hsa" for humans
#'
#' @param ploidy.norm ploidy to normalize (e.g. if 2, then CN = 2 is now 0)
#'
#' @param gene.idtype type of gene id in data matrix. GAC uses hgnc.symbol
#'
#' @param out.suffix filename output suffix, will be appendd after the clones
#'
#' @param node.sum how to summarise node data when gene has multiple mappings.
#' Possible values are: "sum","mean", "median", "max", "max.abs" and "random".
#' Passed to \code{\link[pathview]{node.map}}.  default is 'max.abs'.
#'
#' @param ... additional parameters passed to pathview
#'
#' @return
#'
#' Returns a normalized CN representation of a KEGG pathway.  By default, the scale
#' shows the copy number deviation from the expected ploidy of the tumor.
#'
#' Output is written to the working directory as a PNG, among other files
#'
#' @examples
#' data(cnr)
#' data(paths.hsa, package = "pathview")
#' 
#' noisy.cells <- cnr$qc$cellID[cnr$qc$qc.status == "FAIL"]
#'
#' ## reduced pipeline to genrate DDRC clone profiles
#' cnr <- excludeCells(cnr, excl = noisy.cells)
#' cnr <- phyloCNR(cnr, root.cell = "cell0")
#' cnr <- setBrayClusters(cnr)
#' cnr <- consensusClusterCNR(cnr, iters = 20, maxK = 40)
#' cnr <- doKSpectral(cnr)
#' cnr <- setKcc(cnr)
#' cnr <- cluster_heterogeneity(cnr, by = "category1",
#'           cluster_column = "ConsensusC")
#' cnr <- get_cluster_profiles(cnr)
#' 
#' clone.1 <- colnames(cnr$DDRC.g)[1]
#'
#' top.3.clones <- names(rev(sort(table(cnr$Y$final_cluster))))[1:3]
#' 
#' pathview_clones(cnr, clones = clone.1)
#' pathview_clones(cnr, clones = top.3.clones)
#' 
#' @references
#' If you use this module, please cite:
#' 
#' Weijun Luo and Cory Brouwer. Pathview: an R/Bioconductor package for
#' pathway-based data integration and visu- alization. Bioinformatics,
#' 29(14):1830-1831, 2013. doi: 10.1093/bioinformatics/btt285.
#' 
#' @import pathview
#' @importFrom assertthat assert_that
#' @importFrom utils data
#' 
#' @export
pathview_clones <- function(cnr, clones,
                            pathway.id = "hsa05200",
                            species = "hsa",
                            ploidy.norm = 2,
                            gene.idtype = "SYMBOL",
                            out.suffix = "clones",
                            node.sum = "max", ...) {

    ## check required variables
    assertthat::assert_that(!is.null(clones))
    assertthat::assert_that(is.character(clones))

    ## check required variables
    assertthat::assert_that(!is.null(pathway.id))
    assertthat::assert_that(is.character(pathway.id))

    ## check required variables
    assertthat::assert_that(!is.null(species))
    assertthat::assert_that(is.character(species))

    ## data("bods", package = "pathview")    

    ## checking DDRC.g is present
    if(is.null(cnr[["DDRC.g"]])) {
        stop("No clone profiles available, please run cluster_heterogeneity")
    }

    ## check all clones are in the data
    assertthat::assert_that(all(clones %in% colnames(cnr[["DDRC.g"]])))
    
    ## set up clones string
    clones.txt <- paste(clones, collapse = "_")

    tmp1 <- cnr[["DDRC.g"]][, clones]
    tmp1 <- tmp1 - ploidy.norm
    tmp1[tmp1 == -2] <- (ploidy.norm *2) * -1
    tmp1[tmp1 == -1] <- ploidy.norm * -1

    ## fx
    pwv <- pathview::pathview(gene.data = tmp1,
             pathway.id = pathway.id,
             species = species,
             out.suffix = paste(clones.txt, out.suffix,
                                sep = "_"),
             kegg.native = TRUE,
             same.layer = FALSE, gene.idtype = "SYMBOL", limit = c(-4, 4),
             high = "#D40000FF", mid = "#FFFFFFFF", low = "#005EB8FF",
             node.sum = node.sum,
             ...)
    return(pwv)
} # end pathview_clones
