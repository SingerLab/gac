#' Run Specral K to select optimum and maximum K from Consensus Clustering
#'
#' @param cnr a cnr bundle, with `ccp`
#'
#' @param ... additional parameters passed 
#'
#' @author Nick Socci <soccin@mskcc.org>
#'
#' @return
#'
#' doKSpectral returns two matrices
#'
#' kStats a dataframe containing number of stable K clusters (spectral K) for every kCC (K parameter), and the maximum delta Lamdba.
#'
#' eigenVals a data frame containing the eigen values for each kCC
#'
#' Lastly, it returns optK, a vector containing the recomended kCC and the value of `k` i.e. (spectral K). These are the kCC for the max number of k and max dLambdaMax
#' 
#' @references  Philip A. Knight (2008) The Sinkhorn–Knopp Algorithm: Convergence and Applications. SIAM Journal on Matrix Analysis and Applications 30(1), 261-275. doi: 10.1137/060659624
#' 
#' @import dplyr
#' 
#' @export
doKSpectral <- function(cnr, ...) {
    
    kStats <- list()
    eigenVals <- list()
    
    for(kCC in 2:length(cnr[["ccp"]])) {
    
        message("kCC = ", kCC, "...")
        S <- cnr[["ccp"]][[kCC]]$consensusMatrix
        kS <- kSpectral(S)
    
        kStats[[kCC]] <- as_tibble(kS)[1, ] %>%
            mutate(kCC = kCC) %>%
            select(-.data$topEigenValues)

        eigenVals[[kCC]] <- as_tibble(kS) %>%
            mutate(kCC = kCC) %>%
            select(.data$kCC, .data$k, .data$topEigenValues)
    
        if(kS$k < 4 & kCC > 15) {
            break
        }
    }
    
    cnr[["kStats"]] <- do.call(rbind, kStats)
    cnr[["eigenVals"]] <- do.call(rbind, eigenVals)
    
    cnr[["optK"]] <- c("kCC" = cnr[["kStats"]] %>%
                           top_n(1, .data$k) %>%
                           top_n(1, .data$dLambdaMax) %>%
                           pull(.data$kCC),
                       "sK" = cnr[["kStats"]] %>%
                           top_n(1, .data$k) %>%
                           top_n(1, .data$dLambdaMax) %>%
                           pull(.data$k))
    
    return(cnr)
    
}   ## doKSpectral

