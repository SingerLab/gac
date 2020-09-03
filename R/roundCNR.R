#' rounded quantal matrix
#'
#' Rounds an X matrix to the nearest integer with special
#' threasholds specific to single-cell quantal data from the Varbin algorithm
#' for deletions, losses, and neutral (2 copies)
#' 
#'
#' @param X a numerical matrix composed of bins with copy numbers estimates integer or numerical
#' @param neut maximum value that will be floored to a neutral copy number for diploid species i.e. 2
#' @param loss lower bound for neutral, and upper bound threshold to be considered a loss of one copy
#' @param del maximum threshold to consider a total deletion of the bin i.e. 0 copies
#'
#' @return
#' Rounds X based on quantal
#' 
#' @examples
#'
#' data(copynumbers)
#'
#' cni <- roundCNR(copynumbers)
#' 
#' Heatmap(cni)
#' 
#' 
#' @export
roundCNR <- function(X, neut = 2.5, loss = 1.2, del = 0.2) {
    X[X <= neut & X > loss] <- 2
    X[X <= loss & X > del] <- 1
    X[X <= del] <- 0
    X <- round(X)
    X
}

