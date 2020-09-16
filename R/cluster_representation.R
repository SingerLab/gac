#' cluster representation across the sample
#'
#' @param cnr a cnr bundle
#'
#' @param by a categorical variable used to stratify the cell population
#' 
#' @import entropy
#' 
#' @export
cluster_representation <- function(cnr, by = NULL) {

    assertthat::is.string(by)
    assertthat::assert_that(by %in% names(cnr$Y))

    cc <- table(cnr$Y[, c("cluster", by)])
    cc <- cbind(cc, rowSums(cc))
    colnames(cc)[ncol(cc)] <- "n.cells"
    cc <- data.frame(cc)
    
    cc$n.samples <- rowSums(binary.mat(cc[,1:(ncol(cc)-1)]))
    cc$entropy <- apply(cc[,1:(ncol(cc)-2)], 1, entropy::entropy)

    cnr[["consensus_counts"]] <- cc
    
    return(cnr)
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

