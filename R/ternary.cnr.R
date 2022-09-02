#' build ternary matrix from integer copy number data
#'
#' This function builds a ternary matrix from the cnr$X or cnr$genes.  It was designed as a helper function to generate the input for infSCITE.
#'
#' Because infSCITE was developed for mutation data, this function attempts to recreate the scenario where 0 is no alteration, 1 is a low level / impact alteration, and 2 is high level / impact alteration.
#' 
#' By default, a deletion i.e. copy number of 0 will be treated as a high level deletion, and a loss i.e. one copy as a 1.  The diploid state is treated as 0 = No mutation
#'
#'
#' @param X a X matrix either the bins or genes
#'
#' @param gain an integer copy number value specifying the minimum number of copies to call a gain as 1
#'
#' @param amp an integer value to specify the minimum number of copies to call as a high level amplification e.g. double minute
#'
#' @examples
#' 
#' \dontrun{
#'
#' data(cnr)
#'
#' G <- ternary.cnr(cnr$genes[, c("CDK4", "MDM2")])
#'
#' }
#' 
#' @keywords internal
#' @noRd
ternary.cnr <- function(X, gain = 3, amp = 20) {
    G <- X
    G[G == 2] <- 0
    G[G == 0] <- 2
    G[G == 1] <- 1
    G[G >= gain & G <= amp] <- 1
    G[G > amp] <- 2
    G
} ## ternary.cnr
