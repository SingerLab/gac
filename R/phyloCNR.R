#' Calculating cell-to-cell distances, heirarchical clustering, and generating a tree class `phylo`
#'
#' @param cnr a cnr bundle
#'
#' @param dist.method method for calculating cell-to-cell distance
#' (see vegan::vegdist)
#'
#' @param hclust.method method for heirarchical clustering (see hclust)
#'
#' @param root.cell a cellID to root the three.
#'
#' @param ... other parameters
#'
#' @return
#'
#' Creates a cell-to-cell distance matrix, runs heirarchical clustering,
#'  and converts the `hclust` object to an `ape` class `phylo` object.
#'
#' \itemize{
#'   \item cdb cell to cell Bray-Curtis dissimiarly 
#'   \item hcdb heirarchical clustering of distance matrix
#'   \item phylo ape class `phylo` object
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom ape as.phylo root
#' 
#' @export
phyloCNR <- function(cnr, dist.method = "bray", hclust.method = "ward.D2", root.cell = NULL, ...) {

    if(!is.null(root.cell)) {
        assertthat::assert_that(length(root.cell) == 1)
        assertthat::assert_that(any(root.cell %in% names(cnr$X)),
                                msg = "root.cell not found. the root.cell must be contained within the data")
    }
    
    cnr[["cdb"]] <- distCNR(cnr, method = dist.method, ...)

    cnr[["hcdb"]] <- hclustCNR(cnr, method = hclust.method, ...)

    cnr[["phylo"]] <- ape::as.phylo(cnr[["hcdb"]])
    
    if(!is.null(root.cell)) {
        cnr[["phylo"]] <- ape::root(cnr[["phylo"]], outgroup = root.cell,
                                    resolve.root = TRUE, ...)
    }
    
    return(cnr)
    
}



#' Calculating clone-to-clone distances, and estimating a phylogenetic tree
#' 
#' @param cnr a cnr bundle
#'
#' @param exclude.clones vector of clones to exclude
#'
#' @param root.clone a single clone to root the tree
#'
#' @param dist.method method for calculating cell-to-cell distance
#' (see vegan::vegdist)
#'
#' @param hclust.method method for heirarchical clustering (see hclust)
#'
#' @param ... other parameters
#'
#' @return
#'
#' Creates a cell-to-cell distance matrix, runs heirarchical clustering,
#'  and converts the `hclust` object to an `ape` class `phylo` object.
#'
#' \itemize{
#'   \item cdb cell to cell Bray-Curtis dissimiarly 
#'   \item hcdb heirarchical clustering of distance matrix
#'   \item phylo ape class `phylo` object
#' }
#'
#' @importFrom assertthat assert_that
#' @importFrom ape as.phylo root nj
#' @importFrom vegan vegdist
#' 
#' @export
phyloDDRC <- function(cnr, exclude.clones = NULL, root.clone = NULL,
                      dist.method = "manhattan", hclust.method = "ward.D2",  ...) {
    
    assertthat::assert_that(!is.null(cnr$DDRC.df),
                            msg = "DDRC.df not found, please run get_cluster_profiles")
    
    if(!is.null(root.clone)) {
        assertthat::assert_that(length(root.clone) == 1)
        assertthat::assert_that(any(root.clone %in% colnames(cnr$DDRC.df)),
                                msg = "root.cell not found. the root.cell must be contained within the data")
    }
    
    if(!is.null(exclude.clones)) {
        assertthat::assert_that(all(exclude.clones %in% colnames(cnr$DDRC.df)),
                                msg = "not all exclude.clones are present in DDRC.df")
        keep.clones <- setdiff(colnames(cnr$DDRC.df), exclude.clones)
    } else {
        keep.clones <- colnames(cnr$DDRC.df)
    }
    
    cnr[["DDRC.dist"]] <- vegan::vegdist(cnr$DDRC.df[, keep.clones], method = dist.method, ...)
    
    cnr[["DDRC.nj"]] <- ape::nj(cnr[["DDRC.dist"]])
    
    if(!is.null(root.clone)) {
        assertthat::assert_that(length(root.clone) == 1,
                                msg = "there is more than one root clone")
        
        cnr[["DDRC.nj"]] <- ape::root(cnr[["DDRC.nj"]], outgroup = root.clone,
                                      resolve.root = TRUE, ...)
    }
    
    return(cnr)
    
} # end phyloDDRC

