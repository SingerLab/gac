#' Calculating cell-to-cell distances, heirarchical clustering, and generating a tree class `phylo`
#'
#' @param cnr a cnr bundle
#'
#' @param root.cell a cellID to root the three.
#'
#' @param dist.method method for calculating cell-to-cell distance
#' (see \link[vegan]{vegdist})
#'
#' @param hclust.method method for heirarchical clustering (see hclust)
#'
#' @param tree.method minimum evolution phylogenetics method, can be
#' `bal`, `ols` or NULL. Default is `bal`.
#' 
#' @param ... other parameters passed to \link[vegan]{vegdist}
#'
#' @return
#'
#' Creates a cell-to-cell distance matrix, runs heirarchical clustering.
#'  By defaul cell phylogenetics is infered by \link[ape]{fastme.bal},
#'  and alternatively by `fastme.ols`. When `tree.method = NULL, the
#' `hclust` object is converted  to an `ape` class `phylo` object to
#'  represent the cell phylogeny.
#' 
#' \itemize{
#'   \item cdb cell to cell Bray-Curtis dissimiarly 
#'   \item hcdb heirarchical clustering of distance matrix
#'   \item phylo ape class `phylo` object
#' }
#'
#' @examples
#'
#' data(cnr)
#'
#' ## unrooted cell phylogenetic tree
#' cnr <- phylo_cnr(cnr)
#'
#' cnr$phylo
#' plot(cnr$phylo)
#' 
#'
#' ## rooted cell phylogenetic tree
#' cnr <- phylo_cnr(cnr, root.cell = "cell0")
#'
#' cnr$phylo
#' plot(cnr$phylo)
#'
#' @importFrom assertthat assert_that
#' @importFrom ape as.phylo fastme.bal fastme.ols root 
#'
#' @references
#' Vincent Lefort,  Richard Desper,  Olivier Gascue. 2015. "FastME 2.0:
#'   A Comprehensive, Accurate, and Fast Distance-Based Phylogeny Inference
#'   Program". Molecular Biology and Evolution, Volume 32, Issue 10, October
#'   2015, Pages 2798–2800. <https://doi.org/10.1093/molbev/msv150?>
#'
#' Paradis E, Schliep K (2019). “ape 5.0: an environment for modern
#'   phylogenetics and evolutionary analyses in R.” Bioinformatics,
#'   *35*, 526-528.  <https://doi.org/10.1093/bioinformatics/bty633>.
#' 
#' Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P,
#'   O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M,
#'   Bedward M, Bolker B, Borcard D, Carvalho G, Chirico M, De Caceres M,
#'   Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan
#'   G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T,
#'   Stier A, Ter Braak C, Weedon J (2022). _vegan: Community Ecology
#'   Package_. R package version 2.6-4,
#'   <https://CRAN.R-project.org/package=vegan>.
#'
#' 
#' @export
phylo_cnr <- function(cnr, root.cell = NULL, dist.method = "bray",
                     hclust.method = "ward.D2",  tree.method = "bal",
                     ...) {

    if(!is.null(tree.method)) {
        assertthat::assert_that(tree.method %in% c("bal", "ols"))
    }
    
    if(!is.null(root.cell)) {
        assertthat::assert_that(length(root.cell) == 1)
        assertthat::assert_that(any(root.cell %in% names(cnr$X)),
                                msg = "root.cell not found. the root.cell must be contained within the data")
    }
    
    cnr[["cdb"]] <- distCNR(cnr, method = dist.method, ...)

    cnr[["hcdb"]] <- hclustCNR(cnr, method = hclust.method, ...)

    if(is.null(tree.method)) {
        cnr[["phylo"]] <- ape::as.phylo(cnr[["hcdb"]])
    }
    if(tree.method == "bal") {
        cnr[["phylo"]] <- ape::fastme.bal(cnr[["cdb"]])
    }
    if(tree.method == "ols") {
        cnr[["phylo"]] <- ape::fastme.ols(cnr[["cdb"]])
    }
    

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
#' @param tree.method method for phylogenetic tree construction
#'  either "bal", "ols", or "nj". default is "bal"
#'
#' @param spr parameter specific to \link[ape]{fastme.bal}. Performs
#' Subtree Pruning and Regrafting.  Ignored in all other methods.
#'  Default is TRUE
#' 
#'
#' @param ... other parameters passed to root
#'
#' @return
#'
#' Creates a cell-to-cell distance matrix, runs heirarchical clustering.
#'  By default a phylogenetic tree is constructed using \link[ape]{fastme.bal},
#' alternative algorithms are \link[ape]{fastme.ols} or \link[ape]{nj}.
#'
#' \itemize{
#'   \item DDRC.dist clone to clone Bray-Curtis dissimiarly 
#'   \item DDRC.phylo phylogenetic tree of class `phylo` from R/ape
#' }
#'
#' @examples
#'
#' data(cnr)
#'
#' cnr <- phylo_cnr(cnr)
#' cnr <- setBrayClusters(cnr)
#' cnr <- get_cluster_profiles(cnr, cluster.column = "BrayC")
#' 
#' ## unrooted cell phylogenetic tree
#' cnr <- phylo_ddrc(cnr)
#'
#' cnr$DDRC.phylo
#' plot(cnr$DDRC.phylo)
#' 
#'
#' ## rooted cell phylogenetic tree
#' cnr <- phylo_ddrc(cnr, root.clone = "C2")
#'
#' cnr$DDRC.phylo
#' plot(cnr$DDRC.phylo)
#'
#' ## change method to ols
#' cnr <- phylo_ddrc(cnr, root.clone = "C2", tree.method = "ols")
#'
#' cnr$DDRC.phylo
#' plot(cnr$DDRC.phylo)
#'
#' @references
#' Vincent Lefort,  Richard Desper,  Olivier Gascue. 2015. "FastME 2.0:
#'   A Comprehensive, Accurate, and Fast Distance-Based Phylogeny Inference
#'   Program". Molecular Biology and Evolution, Volume 32, Issue 10, October
#'   2015, Pages 2798–2800. <https://doi.org/10.1093/molbev/msv150?>
#' 
#' Paradis E, Schliep K (2019). “ape 5.0: an environment for modern
#'   phylogenetics and evolutionary analyses in R.” Bioinformatics,
#'   *35*, 526-528.  <https://doi.org/10.1093/bioinformatics/bty633>.
#' 
#' Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P,
#'   O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M,
#'   Bedward M, Bolker B, Borcard D, Carvalho G, Chirico M, De Caceres M,
#'   Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan
#'   G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T,
#'   Stier A, Ter Braak C, Weedon J (2022). _vegan: Community Ecology
#'   Package_. R package version 2.6-4,
#'   <https://CRAN.R-project.org/package=vegan>.
#'
#' @importFrom assertthat assert_that
#' @importFrom ape as.phylo root nj fastme.bal fastme.ols
#' @importFrom vegan vegdist
#' 
#' @export
phylo_ddrc <- function(cnr, root.clone = NULL,
                      exclude.clones = NULL, 
                      dist.method = "bray",
                      tree.method = "bal",
                      spr = TRUE, ...) {
    
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
    
    cnr[["DDRC.dist"]] <- vegan::vegdist(t(cnr$DDRC.df[, keep.clones]), method = dist.method, ...)
    
    if(tree.method == "bal") {
        cnr[["DDRC.phylo"]] <- ape::fastme.bal(cnr[["DDRC.dist"]], spr = spr)
    }
    if(tree.method == "ols") {
        cnr[["DDRC.phylo"]] <- ape::fastme.ols(cnr[["DDRC.dist"]])
    }
    if(tree.method == "nj") {
        cnr[["DDRC.phylo"]] <- ape::nj(cnr[["DDRC.dist"]])
    }
    
    if(!is.null(root.clone)) {
        assertthat::assert_that(length(root.clone) == 1,
                                msg = "there is more than one root clone")
        
        cnr[["DDRC.phylo"]] <- ape::root(cnr[["DDRC.phylo"]], outgroup = root.clone,
                                      resolve.root = TRUE, ...)
    }
    
    return(cnr)
    
} # end phylo_ddrc


#' calculating cell-to-cell distances for clustering
#'
#' @param cnr a cnr object
#'
#' @param method method for calculating distances, defaults to Bray-Curtis
#' dissimilarity
#'
#' @param ... other parameters passed to vegan::vegdist
#'
#' @importFrom vegan vegdist
#' 
#' @return
#' Returns a cell-to-cell distance matrix of class `dist`
#' 
#' @keywords internal
#' @noRd
distCNR <- function(cnr, method = "bray", ...) {
    
    if(cnr$bulk) {
        
        ## Bray-Curtis dissimilarity doesn't work well on log2 ratio,
        ## transforming back to ratio
        cdb <- vegan::vegdist(t(2^cnr[["X"]]), method = method, ...)

    } else {

        cdb <- vegan::vegdist(t(cnr[["X"]]), method = method, ...)

        return(cdb)

    }
} ## distCNR


#' heirarchical clustering
#'
#' @param cnr a cnr object
#'
#' @param method method for heirarchical clustering, defaults to "ward.D2"
#'
#' @param ... other parameters passed to hclust
#'
#' @return
#' Returns a heirarchical clustering object
#'
#' @importFrom stats hclust
#' 
#' @keywords internal
#' @noRd
hclustCNR <- function(cnr, method = "ward.D2", ...) {

    hcdb <- stats::hclust(cnr[["cdb"]], method = method, ...)
    
    return(hcdb)

} ## hclustCNR


