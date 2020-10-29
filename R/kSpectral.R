#' kSpectral
#'
#' @param S similarity matrix
#'
#' @author Nick Socci <soccin@mskcc.org>
#'
#' @references Philip A. Knight (2008) The Sinkhorn–Knopp Algorithm: Convergence and Applications. SIAM Journal on Matrix Analysis and Applications 30(1), 261-275. doi: 10.1137/060659624
#' 
#' @export
kSpectral <- function(S) {

    tol <- eps(nrow(S))
    P <- symmetricSinkhornKnopp(S)
    ee <- eigen(P)
    deltaLambda <- -diff(head(ee$values, 25))
    kMax <- which.max(deltaLambda)
    dLambdaMax <- deltaLambda[kMax]

    numLambda1 <- sum(abs(ee$values-1) < tol)

    k <- kMax - numLambda1 + 1

    list(k = k,
        dLambdaMax = dLambdaMax,
        kMax = kMax,
        numLambda1 = numLambda1,
        topEigenValues = ee$values[1:min(3*kMax, nrow(S)/2)]
        )

}

#' eps : Spacing of floating point numbers
#'
#' @param N number of rows of S
#' #' 
#' @author Nick Socci <soccin@mskcc.org>
#' 
#' @references Philip A. Knight (2008) The Sinkhorn–Knopp Algorithm: Convergence and Applications. SIAM Journal on Matrix Analysis and Applications 30(1), 261-275. doi: 10.1137/060659624
#'
#' @export
eps <- function(N) {

    2^(-52+floor(log2(N)))

}


#' sinkhornKnopp
#'
#' @param A Asymmetric matrix
#'
#' @param tol tol default is NA, which computes `tol = eps(N)`
#'
#' @param maxiter maximum number of iterations, default is Inf
#'
#' @param debug logical, turn on/off debugging i.e. TRUE/FALSE; default is FALSE
#'
#' @references Philip A. Knight (2008) The Sinkhorn–Knopp Algorithm: Convergence and Applications. SIAM Journal on Matrix Analysis and Applications 30(1), 261-275. doi: 10.1137/060659624
#' 
#' @export
sinkhornKnopp <- function(A, tol = NA, maxiter = Inf, debug = FALSE) {

    ## N = size(A, 1);
    N <- nrow(A)
    if(is.na(tol)) {
        tol <- eps(N)
    }

    iter <- 1
    cl <- t(1/colSums(A))
    r <- 1/(A %*% t(cl))


    ## subsequent iterations include test
    while(iter < maxiter) {

        iter <- iter + 1
        cinv <- t(r) %*% A

        if(debug) print(c(iter, maxiter, max(abs(cinv * cl - 1)), tol))

        if(max(abs(cinv * cl - 1)) <= tol) {
            break
        }

        cl <- 1/cinv
        r <- 1/(A %*% t(cl))

    }

    A * (r %*% cl)

}

#' symmetricSinkhornKnopp
#' 
#' @param A symmetric matrix
#'
#' @param tol tol default is NA, which computes `tol = eps(N)`
#'
#' @param maxiter maximum number of iterations, default is Inf
#'
#' @param debug logical, turn on/off debugging i.e. TRUE/FALSE; default is FALSE
#'
#' @references Philip A. Knight (2008) The Sinkhorn–Knopp Algorithm: Convergence and Applications. SIAM Journal on Matrix Analysis and Applications 30(1), 261-275. doi: 10.1137/060659624
#'
#' @export
symmetricSinkhornKnopp <- function(A, tol = NA, maxiter = Inf, debug = FALSE) {

    N <- nrow(A)
    if(is.na(tol)) {
        tol <- eps(N)
    }

    ## Force matrix to be symmetric otherwise
    ## Return error it not
    B = (A + t(A))/2

    if( sum(abs(B-A)) > tol ) {
        stop("Passed a non-symmtrix matrix symmetricSinkhornKnopp")
    }

    skM = sinkhornKnopp(B,tol,maxiter,debug)

    ## Re-symmetrize output to deal with roundoff
    (skM + t(skM))/2

}

