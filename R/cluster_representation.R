#' cluster representation across the sample
#'
#' @param cnr a cnr bundle
#'
#' @export
consensus_representation <- function() {
    
}


#' binary matrix 
#'
#' @param mat a matrix
#'
#' @export
binary.mat <- function(mat) {
    mat[mat >= 1] <- 1
    mat
}

