#' Set clusters using one-cell and multi-cell height intersection
#'
#' @param cnr a cnr bundle
#'
#' @param tree.height height of the tree \link{optClust}
#'
#' @param prefix a prefix to append before the cluster number
#'
#' @param ... additional parameters passed to \link{optClust}
#' 
#' @return
#'
#' Returns a cnr object with cluster membership based on Bray-Curtis
#'  Ward.D2 heirarchical clustering. If hieght is NULL, then the
#'  the tree height is calculated using \link{minimum.intersect}. Where
#'  optimal height is the value at the intersect between the number of
#'  one-cell and multi-cell clusters.
#' 
#' 
#' @importFrom assertthat assert_that
#' @importFrom stats cutree
#' 
#' @export
setBrayClusters <- function(cnr, tree.height = NULL, prefix = "C", ...) {
    
##    assertthat::assert_that("BrayC" %in% colnames(cnr[["Y"]]) == FALSE,
##                            msg = "a cluster column already exists")

    if(is.null(tree.height)) {
        mocp <- optClust(cnr, ...)
        tree.height <- minimum.intersect(mocp)
        message("tree.height not set, using minum intersect point of ",
                tree.height, " as tree.height")
    } else {
        tree.height <- tree.height
    }
    
    BrayC <-  cutree(cnr[["hcdb"]], h = tree.height)
    
    if(!is.null(prefix)) {
        BrayC <- paste0(prefix, BrayC)
    }
    
    cnr[["Y"]]$BrayC <- BrayC
    cnr[["tree.height"]] <- tree.height
    
    return(cnr)
}

#' Set cluster membership for K clusters
#'
#'
#' @param cnr a cnr object
#'
#' @param kCC number of K clusters to use, 
#'  select optimum or maximum k from kSpectral.  Default is NULL
#'  which will pull the optK[kCC] from kStats
#'
#' @param prefix prefix charcter to append to Consensus Clusters
#' 
#' @importFrom assertthat assert_that
#' 
#' @return
#'
#' returns cluster membership based on consensus clustering for a
#'  specified kCC
#' 
#' @export
setKcc <- function(cnr, kCC = NULL, prefix = "X") {

    ## assertthat::assert_that("ConsensusC" %in% colnames(cnr[["Y"]]),
    ##                        msg = "a consensus cluster column already exists")
    
    if(is.null(kCC) & !is.null(cnr[["optK"]]["kCC"])) {
        kCC <- cnr[["optK"]]["kCC"]
        message("Using recommended kCC = ", kCC, " by kSpectral")
    } else {
        assertthat::assert_that(kCC < length(cnr[["ccp"]]))
    }
    
    cc_membership <- as.data.frame(factor(cnr[["ccp"]][[kCC]]$consensusClass))
    colnames(cc_membership) <- "ConsensusC"

  if(!is.null(prefix)) {
        cc_membership[,1] <- paste0(prefix, cc_membership[,1])
    }

    cnr <- addPheno(cnr, df = cc_membership, by.x = "cellID", by.y = 0,
                    sort = FALSE)

  cnr[["kCC"]] <- kCC

  return(cnr)
}

#' rank clones
#'
#' @param cnr a cnr object
#'
#' @param cluster_column which cluster column to use; must be one of
#'   "BrayC", "ConsensusC", or "final_cluster"
#' 
#' @param rank_by option to rank, 'fga' to rank by genomic complexity based
#'   on the fraction of genome altered. Or 'frequency', ranked from most to
#'   least least prevalent.
#'
#' @param fga.method option to estimate cluster FGA, either 'genomic' or 'proxy'.
#'   In method 'genomic' the altered bins are weighed by their genomic length in
#'   base pairs (bp), and the fraction is estimated based on the total genome
#'   length.  Method 'proxy' does not weigh by genomic lenght, instead it
#'   approximates FGA by estimating the proportion of altered bins.  Alterations
#'   are deviations from a diploid state e.g. != 2. 
#' 
#' @param prefix prefix charcter to append to Consensus Clusters
#'
#' @return
#'
#' rank_clones returns a ranked_cluster vector to `Y`.  The ranking of the clone is
#' user specified by frequency or genomic complexity.  When ranking clones by frequency,
#' the specified clusters are ranked in descending order, where the most frequent clone
#' is assigned clone `1`, and all less frequent clones are assigned a numerical label
#' in sequential order.  Ties are broke at random.  When ranking by genomic complexity,
#' then the fraction of the genome altered (FGA) is computed based on a diploid genome (2n).
#' And ranked in ascending order, where the least altered cell is clone `1` and more altered
#' clones are assigned a numerical label in sequential order.  Ties are broken at random.
#'
#' Clusters with 3 or more cells are sumarized using the same parameters as in \link{get_cluster_profiles}.
#' However, for clusters with 2 cells, are sumarized by taking the `max` at each bin, and those with one cell, then
#' the cell was used.
#'
#' @importFrom assertthat assert_that
#' 
#' @export
rank_clones <- function(cnr, cluster_column, rank_by, fga.method, prefix = NULL) {
    
    assertthat::assert_that(any(cluster_column %in% c("BrayC", "ConsensusC", "final_cluster")),
                            msg = "cluster_column should be either 'BrayC' 'ConsensusC' _or_ 'final_cluster'")
    
    assertthat::assert_that(any(cluster_column %in% colnames(cnr[["Y"]])),
                            msg = "cluster_column does not exist in cnr[['Y']]")
    
    assertthat::assert_that(!is.null(rank_by),
                            msg = "rank_by is NULL, please select rank option: 'frequency' or 'fga'")

    assertthat::assert_that(!any("ranked_cluster" %in% names(cnr$Y)),
                            msg = "Clusters have already been ranked.  Please remove or rename column `ranked_cluster`")
    
    if(rank_by == "frequency") {
        
        ranked_cluster <- rev(sort(table(cnr[["Y"]][, cluster_column])))
        cnr[["clone_ranking_method"]] <- rank_by

    } else {
        
        if(rank_by == "fga") {
            
            assertthat::assert_that(!is.null(fga.method))
                        
            ucd <- table(cnr[["Y"]][, cluster_column])
            uclust <- ucd[ucd >= 3]
            
            if(length(uclust) != 0) {
                ddrc.df <- sapply(names(uclust), function(clust) {
                    udf <- cnr$X[, cnr$Y$cellID[cnr$Y[, cluster_column] == clust]]
                    if ( ncol(udf) != 0 && !is.null(ncol(udf)) ) {
                        apply(udf, 1, median, na.rm = TRUE)
                    }
                })
            } else {
                ddrc.df <- NULL
            }
            
            uclust <- ucd[ucd == 2]
            
            if(length(uclust) != 0) {
                ndrc.df <- sapply(names(uclust), function(clust) {
                    udf <- cnr$X[, cnr$Y$cellID[cnr$Y[, cluster_column] == clust]]
                    if( ncol(udf) != 0 && !is.null(ncol(udf)) ) {
                        apply(udf, 1, max, na.rm = TRUE)
                    }
                })
            } else {
                ndrc.df <- NULL
            }
            
            uclust <- ucd[ucd == 1]

            if(length(uclust) != 0) {
                udf.cells <- cnr$Y$cellID[cnr$Y[, cluster_column] %in% names(uclust)]
                udrc.df <- cnr$X[, udf.cells]
                colnames(udrc.df) <- names(uclust)
            } else {
                udrc.df <- NULL
            }

            ddrc <- list(ddrc.df, ndrc.df, udrc.df)
            ddrc <- ddrc[-which(sapply(ddrc, is.null))]
            
            ddrc.df <- do.call(cbind, ddrc)

            
            if(fga.method == "genomic") {
                
                assertthat::assert_that(any("bin.length" %in% names(cnr$chromInfo)))
                
                ddrc.fga <- colSums(binary.X(ddrc.df) * cnr$chromInfo$bin.length) / sum(cnr$chromInfo$bin.length)

            } else {
            
                if(fga.method == "proxy") {
                    ddrc.fga <- colSums(binary.X(ddrc.df)) / nrow(cnr$X)
                }
            }
            
            ranked_cluster <- sort(ddrc.fga)
            cnr[["clone_ranking_method"]] <- paste(rank_by, fga.method, sep = ":")
            
        }
        
        ranked_clones <- data.frame(clone = names(ranked_cluster))
        ranked_clones$ranked_cluster <- 1:nrow(ranked_clones)

        if(!is.null(prefix)) {
            ranked_clones$ranked_cluster <- paste(prefix, ranked_clones$ranked_cluster, sep = "")
        }
        
        cellsTmp <- cnr[["Y"]][, c("cellID", cluster_column)]
        cellsTmp <- merge(cellsTmp, ranked_clones, by.x = cluster_column,
                          by.y = "clone", sort = FALSE)
        
        cnr <- addPheno(cnr, df = cellsTmp, sort = FALSE)
        
    }
    
    return(cnr)
    
} ## end ranked_clones
