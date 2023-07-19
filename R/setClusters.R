#' Set Bray clusters 
#'
#' @param cnr a cnr bundle
#'
#' @param tree.height height of the tree \code{\link{optClust}}
#'
#' @param prefix a prefix to append before the cluster number
#'
#' @param opt.method optimization method, must be either "mi" or
#'  'maxp' ; default is 'mi'
#'
#' @param ... additional parameters passed to \link{optClust}
#' 
#' @return
#'
#' Returns a cnr object with cluster membership based on Bray-Curtis
#'  Ward.D2 heirarchical clustering. If method is \code{"mi"} the
#'  the tree height is calculated using \code{minimum.intersect}.
#'  Where optimal height is the value at the intersect between the number of
#'  one-cell and multi-cell clusters.  If method is \code{"maxp"}, 
#'  %CMC and tree height are maximized.
#'
#' @examples
#' data(cnr)
#' noisy.cells <- cnr$qc$cellID[cnr$qc$qc.status == "FAIL"]
#' cnr <- excludeCells(cnr, excl = noisy.cells)
#' 
#' cnr <- phylo_cnr(cnr, root.cell = "cell0")
#' 
#' cnr <- setBrayClusters(cnr)
#'
#' cnr <- setBrayClusters(cnr, opt.method = "maxp")
#' 
#' cnr <- setBrayClusters(cnr, tree.height = 0.16)
#' 
#' @importFrom assertthat assert_that
#' @importFrom stats cutree
#' 
#' @export
setBrayClusters <- function(cnr, tree.height = NULL, prefix = "C",
                            opt.method = "maxp", ...) {
    
    if("BrayC" %in% colnames(cnr[["Y"]])) {
        warning("existing BrayC column will be ovewritten")
    }

    if(is.null(tree.height)) {
        mocp <- optClust(cnr, ...)
        
        if(opt.method == "mi") {
            tree.height <- minimum.intersect(mocp)
        } else {
            tree.height <- maximum.percentage(mocp)
        }
        
        message("tree.height not set, using minimum intersect point of ",
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
#' @param overwrite whether to overwrite the ConsensusC column,
#'   default is TRUE
#' 
#' @importFrom assertthat assert_that
#' 
#' @return
#'
#' Appends cluster membership based on consensus clustering for a
#'  specified kCC to the cnr$Y phenotype matrix.  Column name is
#'  ConsensusC.
#'
#' If `overwrite = TRUE` this column will be overwritten in subsequent
#' analyses.  To prevent this, change the name of the column to a
#' differnt name e.g. `kCC.10` prior to the next run. Or set
#' `overwrite = FALSE`, this will issue a warning, and
#' will rename append `.x` and `.y` suffix from `merge`.
#'
#' @examples
#'
#' data(cnr)
#' noisy.cells <- cnr$qc$cellID[cnr$qc$qc.status == "FAIL"]
#'
#' cnr <- excludeCells(cnr, excl = noisy.cells)
#'
#' cnr <- phylo_cnr(cnr, root.cell = "cell0")
#'
#' cnr <- run_consensus_clustering(cnr, iters = 20, maxK = 40,
#'        verbose = TRUE)
#'
#' cnr <- doKSpectral(cnr)
#'
#' cnr <- setKcc(cnr)
#' 
#' 
#' @export
setKcc <- function(cnr, kCC = NULL, prefix = "X", overwrite = TRUE) {

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
    cc_membership$cellID <- rownames(cc_membership)
    
  if(!is.null(prefix)) {
        cc_membership[,1] <- paste0(prefix, cc_membership[,1])
    }

    
    if("ConsensusC" %in% names(cnr$Y)) {
        if(overwrite) {
            warning("Overwriting ConsensusC")
            ccx <- which(names(cnr$Y) == "ConsensusC")
            cnr$Y <- cnr$Y[, -ccx]
        } else {
            warning("ConsensusC exists, check names in Y matrix before continuing, and modify accordingly")
        }
    } 
    
    cnr <- addPheno(cnr, df = cc_membership, by = "cellID",
                    sort = FALSE)
    
    cnr[["kCC"]] <- kCC
    
    return(cnr)
}


#' Optimizing clustering for single-cell copy number
#'
#' The optClust is a simple attempt to set a threshold on where to cut the an
#' hclust tree.  The aim is to maximize the number of clusters with multiple-cells
#' while minimizing the number of one-cell clusters.
#'
#' @return
#'
#' Returns a matrix where rows are the heights in the opt.range and three columns.
#'
#' One-cell specifies the number of clusters with only one cell, 
#' Multi-cell specifies the number of clusters with multiple cells (>=2), and
#' the percentage of cells in multi-cell clusters (%CMC)
#'
#' A suitable threshold is one that minimizes the one-cell cluster and maximizes
#' the %CMC
#' 
#' @param cnr a cnr bundle
#'
#' @param opt.range  range of tree heights that need to be optimized
#'
#' @examples
#'
#' data(cnr)
#'
#' cnr <- phylo_cnr(cnr)
#'
#' ## macro optimization
#' optClust(cnr, opt.range = seq(0, 0.3, by = 0.05))
#'
#' ## for micro-optimization
#' optClust(cnr, opt.range = seq(0, 0.2, by = 0.001))
#'
#' @export
optClust <- function(cnr, opt.range = seq(0.005, 0.6, by = 0.005)) {

    clL <- sapply(opt.range, function(h) {
        ctbl <- sort(table(cutree(cnr[["hcdb"]], h = h)),
                     decreasing = TRUE)
    })
    names(clL) <- opt.range

    omt <- matrix(t(sapply(clL, function(i) {
        cbind(sum(i == 1), sum(i != 1), sum(i != 1)/length(i))
    })),
    ncol = 3,
    dimnames = list(opt.range, c("One-cell", "Multi-cell", "%CMC")))
    
    return(omt)
    
} # optClust


#' Minimum tree.height at the intersect of one-cell and multi-cell clusters
#'
#' @param mocp output from \link{optClust}
#'
#' @export
minimum.intersect <- function(mocp) {

    mi <- as.numeric(names(which(mocp[,1] < mocp[,2])[1]))
    return(mi)
    
} #end minimum.intersect


#' Max percent
#'
#' @param mocp output from \code{\link{optClust}}
#' 
#' @param low.boundary lower percent of clones with multiple cells,
#'  default is .70
#'
#' @param high.boundary highest percent of clones with multiple cells,
#'  default is .90
#'
#' @return
#'
#' The maximum tree height located at the maximum percent of clones with
#'  multiple cells between the low and high boundaries.
#'
#' @importFrom assertthat assert_that
#' 
#' @export
maximum.percentage <- function(mocp, low.boundary = .70, high.boundary = .90) {

    assertthat::assert_that("%CMC" %in% colnames(mocp))
    
    maxp <- mocp[mocp[, "%CMC"] >= low.boundary , ]
    maxp <- maxp[maxp[, "%CMC"] <= high.boundary, ]

    maxx <- max(maxp[, "%CMC"])
    maxp <- maxp[maxp[, "%CMC"] == maxx,]
    
    max.tree.height <- max(as.numeric(rownames(maxp)))
    
    return(max.tree.height)
    
} #end maximum.percentage


#' rank clones
#' 
#' STILL IN DEVELOPMENT: UNTESTED
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
