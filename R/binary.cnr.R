#' build binary matrix from integer copy number data
#'
#' This function builds a binary, incidence, matrix from the cnr$X.  It was designed as a helper function to generate the input for infSCITE.
#'
#' Because infSCITE was developed for mutation data, this function creates a precense/absence matrix of the data
#' 
#' By default, anything not diploid is 1.
#'
#'
#' @param cnr a cnr matrix either the bins or genes
#'
#' 
#' @examples
#'
#' data(cnr)
#'
#' bimat <- binary.cnr(cnr$genes[, c("CDK4", "MDM2")])
#' 
#' 
#' @export
binary.cnr <- function(cnr) {
    Z <- cnr$X
    Z[Z == 2] <- 0
    Z[Z != 0] <- 1
    Z
} ## binary.cnr
