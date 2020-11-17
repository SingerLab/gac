#' Build a consensus copy number profile for each cluster or clone
#'
#' @param cnr a cnr bundle containing a `final_cluster` column in Y.
#'
#' @param minimum_cells minimum number of cells in a cluster, best if greater than 3
#' to estimate a median
#'
#' @return
#'
#' Function returns the cnr with three additional tables. 
#' * uclust : number of cells in each final_cluster, only clusters greater than the specified minimum number of cells is shown
#' * DDRC.df: a matrix containing the representative profile for each `final_cluster` at the bin level.
#' * DDRC.g : interpolation of the DDRC.df at the gene level for each `final_cluster`.
#' 
#' @export
get_cluster_profiles <- function(cnr, minimum_cells = 3) {

    uclust <- use_clusters(cnr, minimum_cells = minimum_cells)
    
    ## generate consensus profile for each cluster and bind them
    ## this also removes clusteres with a low number of cells
    DDRC.df <- sapply(names(uclust), function(clust) {
    ## getting consensus copy number profile for each cluster
        udf <- cnr$X[, cnr$Y$cellID[cnr$Y$final_cluster == clust]]
        if ( ncol(udf) != 0 && ! is.null(ncol(udf)) ) {
            apply(udf, 1, median, na.rm = TRUE)
        }
    })
    
    DDRC.df[DDRC.df >=2] <- ceiling(DDRC.df[DDRC.df >=2])
    DDRC.df[DDRC.df <2] <- floor(DDRC.df[DDRC.df <2])
    
    DDRC.g <- round(t(expand2genes(DDRC.df, gene.index = cnr$gene.index)))
    
    cnr[["uclust"]] <- uclust
    cnr[["DDRC.df"]] <- DDRC.df
    cnr[["DDRC.g"]] <- DDRC.g

    return(cnr)
    
}


#' select clusters to use based on a minimum number of cells
#'
#' @param cnr a cnr bundle containing the output of setCluster
#'
#' @param minimum_cells minimum number of cells in a cluster, must be greater than 3
#' to estimate a median
#'
#' @export
use_clusters <- function(cnr, minimum_cells = 3) {

    ucd <- table(cnr[["Y"]]$final_cluster)
    ucd <- ucd[ucd >= minimum_cells]
    
    return(ucd)

} # use_clusters

