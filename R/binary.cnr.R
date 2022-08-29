#' build binary CNR from integer copy number data
#'
#' This function builds a binary, incidence, matrix from the cnr$X.  It was designed as a helper function to generate the input for infSCITE.
#'
#' Because infSCITE was developed for mutation data, this function creates a precense/absence matrix of the data
#' 
#' By default, anything not diploid is 1.
#'
#'
#' @param cnr a copy number matrix to convert to an incidence matrix
#' 
#' @param base.ploidy expected cell ploidy, e.g. 2N = 2, 4N = 4.
#' 
#' @return
#'
#' Returns an incidence matrix for X as part of the cnr object
#'
#' \itemize{
#'   \item Z incidence matrix from X, copy number != 2 is 1, all 2 are 0
#' }
#'
#' 
#' @examples
#'
#' data(cnr)
#'
#' Z <- binary.cnr(cnr)
#' 
#' 
#' @export
binary.cnr <- function(cnr, base.ploidy = 2) {

    Z <- cnr[["X"]]
    Z[Z != base.ploidy] <- 1
    Z[Z == base.ploidy] <- 0
    
    cnr[["Z"]] <- Z
    
    return(cnr)
    
} ## binary.cnr



#' build binary matrix from integer copy number data
#'
#' This function builds a binary, incidence, matrix from the cnr$X.  It was designed as a helper function to generate the input for infSCITE.
#'
#' Because infSCITE was developed for mutation data, this function creates a precense/absence matrix of the data
#' 
#' By default, anything not diploid is 1.
#'
#'
#' @param X a copy number matrix to convert to an incidence matrix
#'
#' @param base.ploidy expected cell ploidy, e.g. 2N = 2, 4N = 4
#'
#' @return
#'
#' Returns an incidence matrix for X, when X != base.ploidy
#' 
#' @examples
#'
#' data(cnr)
#'
#' Z <- binary.X(cnr$genes[, c("CDK4", "MDM2")])
#' 
#' 
#' @export
binary.X <- function(X, base.ploidy = 2) {
    Z <- X
    Z[Z != base.ploidy] <- 1
    Z[Z == base.ploidy] <- 0
    Z
} ## binary.X

