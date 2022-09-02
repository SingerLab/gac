#' implementation of SCclust pipeline
#'
#' @param cnr a cnr bundle
#'
#' @param cytobands UCSC cytobands file
#'
#' @param hc.method clustering method for hclust_tree
#'
#' @param tree.method clustering method for tree_py
#' 
#' @param run.fisher logical, weather you want to run sim_fisherCNR. defaults TRUE.
#'    If FALSE, it assumes you have ran sim_fisherCNR before, and the cnr already
#'    contains cnr[["pins"]] and cnr[["fisher"]]
#'
#' @param ... additional arguments passed
#'
#' @source \url{https://github.com/KrasnitzLab/SCclust}
#'
#' @import SCclust
#' 
#' @keywords internal
#' @noRd
scClustCNR <- function(cnr, cytobands, hc.method = "average",
                       tree.method = "average", run.fisher = TRUE, ...) {
    if(run.fisher) {
        cnr <- sim_fisherCNR(cnr, cytobands = cytobands, ...)
    }
    
    cnr <- fisher_tree(cnr, hc.method = hc.method, tree.method = tree.method, ...)
    
    cnr <- getSubclonesCNR(cnr, ...)
    
    return(cnr)

} ## end scClustCNR

    
#' implementation of SCclust sim_fisher_wrapper
#'
#' scClust uses Fisher distance
#'  
#' @param cnr a cnr bundle
#'
#' @param cytobands UCSC cytobands file
#'
#' @param centromere cytoband boundaries of the centromere (will be excluded)
#' 
#' @param nsim number of Fisher wrapper simulations
#' 
#' @param ... additional parameters pased to SCclust::sim_fisher_wrapper
#'
#' @return
#'
#' SCclust pipeline and objects
#' 
#' @import SCclust
#' 
#' @source \url{https://github.com/KrasnitzLab/SCclust}
#'
#' @examples
#'
#' \dontrun{
#' ## the simulated data does not run well
#' data(cnr, hg19_cytoBand)
#' cnr <- sim_fisherCNR(cnr, cytobands = hg19_cytoBand)
#' }
#' 
#' @keywords internal
#' @noRd
sim_fisherCNR <- function(cnr, cytobands, centromere = c("p11", "q11"),
                          nsim = 200, ...) {

    cnr[["cytobands"]] <- cytobands

    cnr[["centroareas"]] <- SCclust::calc_centroareas(cytobands,
                                                      centromere = centromere)
    
    cnr[["centrobins"]] <- SCclust::calc_regions2bins(cnr$chromInfo,
                                                      cnr$centroareas)

    pinX <- cbind(cnr$chromInfo[, c("chrom", "chrompos", "abspos")], cnr$X)
    
    cnr[["pins"]] <- SCclust::calc_pinmat(cnr$chromInfo, pinX, 
                                          dropareas = cnr$centroareas)
                                          
    
    cnr[["fisher"]] <- SCclust::sim_fisher_wrapper(cnr$pins$pinmat,
                                                   cnr$pins$pins,
                                                   nsim = nsim,
                                                   ...)
    
    return(cnr)
}


#' implementation of fisher_dist, hc, and tree
#'
#' @param cnr a cnr bundle
#'
#' @param hc.method method on hcclust_tree, defaults to average
#'
#' @param tree.method method on tree_py
#'
#' @param ... additional parameters pased to SCclust::sim_fisher_wrapper
#'
#' @import SCclust
#'
#' @source \url{https://github.com/KrasnitzLab/SCclust}
#' 
#' @keywords internal
#' @noRd
fisher_tree <- function(cnr, hc.method = "average", tree.method = "average", ...) {

    cnr[["mfdr"]] <- SCclust::fisher_fdr(cnr$fisher$true, cnr$fisher$sim,
                                         cnr$cells, ...)
    cnr[["mdist"]] <- SCclust::fisher_dist(cnr$fisher$true, cnr$cells)

    cnr[["hc"]] <- SCclust::hclust_tree(cnr$pins$pinmat, cnr$mfdr, cnr$mdist,
                               hcmethod = hc.method)
    cnr[["tree_df"]] <- SCclust::tree_py(cnr$mdist, method = tree.method)

    return(cnr)
    
}


#' implementation of finding clones and subclones from SCclust
#'
#' @param cnr a cnr bundle
#'
#' @param ... additional parameters for find_subclones
#'
#' @source \url{https://github.com/KrasnitzLab/SCclust}
#' 
#' @import SCclust
#'
#' @keywords internal
#' @noRd
getSubclonesCNR <- function(cnr,  ...) {

    cnr[["hc"]] <- SCclust::find_clones(cnr[["hc"]])
    cnr[["subclones"]] <- SCclust::find_subclones(cnr[["hc"]], cnr$pins$pinmat,
                                         cnr$pins$pins, ...)

    cnr <- addPheno(cnr, cnr[["subclones"]])
    
    return(cnr)
}

