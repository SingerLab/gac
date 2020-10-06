#' cluster heterogeneity across the sample
#'
#' @param cnr a cnr bundle
#'
#' @param by a categorical variable used to stratify the cell population.  If
#' NULL i.e. no stratification, the representation will be done overall
#' 
#' @import entropy
#'
#' @examples
#'
#' data(cnr)
#'
#' cnr <- phyloCNR(cnr)
#'
#' cnr <- setClusters(cnr, tree.height = 0.065)
#'
#' cnr <- cluster_heterogeneity(cnr)
#'
#' 
#' @export
cluster_heterogeneity <- function(cnr, by = NULL) {

    if(is.null(by)) {

        occ <- table(cnr$Y$cluster)
        ns <- rep(1, length(occ))
        
    } else {
        
        assertthat::is.string(by)
        assertthat::assert_that(by %in% names(cnr$Y))
        
        occ <- table(cnr$Y[, c("cluster", by)])
        ns <- rowSums(binary.mat(occ))
        
    }
    
    occ <- cluster_representation(occ, n.samples = ns)
    
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
#' @param n.samples vector containing the number of samples. Must
#' mach number of rows in cc
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
#' ( cc <- cluster_representation(cc, n.samples = ns) )
#'
#'@export
cluster_representation <- function(cc, n.samples = 1) {

    acc <- ifelse(is.na(ncol(cc)), cc, rowSums(cc))
    
    cc <- data.frame(
        cbind(cc,
              n.cells = cc,
              overall.freq = acc / sum(cc),
              n.samples = n.samples,
              entropy = entropy::entropy(cc))
    )
    
    return(cc)
    
} # end heterogeneity

