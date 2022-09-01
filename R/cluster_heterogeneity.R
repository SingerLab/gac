#' cluster heterogeneity across the sample
#'
#' @param cnr a cnr bundle
#'
#' @param by a categorical variable used to stratify the cell population.  If
#' NULL i.e. no stratification, the representation will be done overall
#'
#' @param cluster_column column to use as cluster
#' 
#' @examples
#'
#' data(cnr)
#'
#' cnr <- phyloCNR(cnr)
#'
#' cnr <- setBrayClusters(cnr, tree.height = 0.065)
#'
#' cnr <- cluster_heterogeneity(cnr, cluster_column = "BrayC")
#'
#' @importFrom assertthat assert_that is.string
#' 
#' @export
cluster_heterogeneity <- function(cnr, by = NULL, cluster_column = NULL) {

    assertthat::assert_that(!is.null(cluster_column))
    assertthat::assert_that(any(cluster_column %in% names(cnr$Y)))
    
    if(is.null(by)) {

        occ <- table(cnr$Y[, cluster_column])
        ns <- rep(1, length(occ))
        
    } else {
        
        assertthat::is.string(by)
        assertthat::assert_that(by %in% names(cnr$Y))
        
        occ <- table(cnr$Y[, c(cluster_column, by)])
        ns <- rowSums(binary.mat(occ))
        
    }
    
    occ <- cluster_representation(occ)

    cnr$Y$final_cluster <- cnr$Y[, cluster_column]

    cnr[["cluster_heterogeneity"]] <- occ
    
    return(cnr)
} # cluster_heterogeneity


#' binary matrix 
#'
#' @param mat a matrix
#'
#' @export
binary.mat <- function(mat) {
    mat[mat >= 1] <- 1
    mat
} # end binary.mat


#' cluster representation
#'
#' @param cc cluster counts
#'
#' @examples
#' 
#' ## cluster table, no subsamples
#' cc <- table(sample(paste0("C", 1:12), size = 1000, replace = TRUE))
#' 
#' ( cc <- cluster_representation(cc) )
#'
#' cc <- data.frame(cluster = sample(paste0("C", 1:12), size = 1000,
#'                                   replace = TRUE),
#'             sample = sample(paste0("S", 1:6), size = 1000, replace = TRUE))
#' 
#' cc <- table(cc[, c("cluster", "sample")])
#' cc[cc <= 10] <- 0
#'
#' ns <- rowSums(binary.mat(cc))
#' 
#' ( cc <- cluster_representation(cc) )
#'
#' @importFrom entropy entropy
#' 
#' @keywords internal
#' 
#' @export
cluster_representation <- function(cc) {

    if(is.matrix(cc)) {
        n.cells.clone <- rowSums(cc)
        freqs <- cc / n.cells.clone
        colnames(freqs) <- paste0(colnames(cc), ".fq")
        n.regions <- rowSums(binary.mat(cc))
        enropy <- apply(cc, 1, entropy::entropy)
        
    } else {
        n.cells.clone <- sum(cc)
        freqs <- cc / n.cells.clone
        n.regions <- 1
        entropy <- entropy::entropy(cc)
    }
    
    ccx <- data.frame(
        cbind(cc,
              n.cells = n.cells.clone,
              n.regions = n.regions,
              freqs,
              entropy = entropy)
    )
    
    return(ccx)
    
} # end heterogeneity

