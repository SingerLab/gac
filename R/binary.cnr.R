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
#' 
#' @examples
#'
#' data(cnr)
#'
#' Z <- binary.cnr(cnr$genes[, c("CDK4", "MDM2")])
#' 
#' 
#' @export
binary.cnr <- function(X) {
    Z <- X
    Z[Z == 2] <- 0
    Z[Z != 0] <- 1
    Z
} ## binary.cnr
