#' Copy number call categories from ratio data
#'
#' @param cnr a cnr bundle
#'
#' @param cbioportal logical, weather to use the cBioPortal copy number
#' call categories of -2, -1, 0, 1, 2 to denote deletion, loss, neutral,
#' gain, amplificiation, respectively. Default is FALSE.  If TRUE,
#' hyper amplification calls will be joined with 2
#' 
#' @param hyper.amplification logical, weather to call hyper amplifications where
#' the number of copies exceeds the upper bound amplification threshold defaults
#' to TRUE, and uses 20 as the default. Parameter is ignored if cBioPortal
#' notation is TRUE
#' 
#' @param deletion.threshold integer to use as locus deletion, default = 0
#' 
#' @param loss.threshold number of copies to use as copy number loss, default = 1
#' 
#' @param base.ploidy integer to denote base ploidy to use as Neutral
#' 
#' @param gains.thresholds number of copies to use as copy number gains,
#'  default to a ower bound of 3 (+1 copy from 2N) and upper bound of 5 (+3 copies)
#'
#' @param amplification.thresholds number of copies to use as copy number
#' amplification defaults to a lower bound of 6 (+4 copies) and upper bound of
#' 20 (+18 copies)
#'
#' @return
#' With default parameters, the function returns a matrix with copy number
#' calls in the form of D, L, N, G, A and A2 from an genotype X matrix.  
#'
#' | CN    |  Call  | cBioPortal Notation |
#' |------:|: ---- :|: ----------------- :|
#' |     0 |   D    | -2                  |
#' |     1 |   L    | -1                  |
#' |     2 |   N    |  0                  |
#' |   3-5 |   G    |  1                  |
#' |  6-20 |   A    |  2                  |
#' |   >20 |   A2   |  2                  |
#'
#' @examples
#' data(cnr)
#'
#' Xc <- callCN(cnr)
#' Xc <- callCN(cnr, hyper.amplification = FALSE)
#'
#' Xp <- callCN(cnr, cbioportal = TRUE)
#'
#' @importFrom assertthat assert_that
#' 
#' @export
callCN <- function(cnr,
                   cbioportal = FALSE,
                   hyper.amplification = TRUE,
                   deletion.threshold = 0,
                   loss.threshold = 1,
                   base.ploidy = 2,
                   gains.thresholds = c(3, 5),
                   amplification.thresholds = c(6, 20)) {
    ## checks
    assertthat::assert_that(loss.threshold > deletion.threshold)
    assertthat::assert_that(base.ploidy > loss.threshold)
    assertthat::assert_that(gains.thresholds[1] > base.ploidy)
    assertthat::assert_that(gains.thresholds[1] <= gains.thresholds[2])
    assertthat::assert_that(amplification.thresholds[1] > gains.thresholds[2])
    assertthat::assert_that(amplification.thresholds[1] <=
                            amplification.thresholds[2])
    
    mm = matrix(NA,
                nrow= nrow(cnr$X),
                ncol = ncol(cnr$X),
                dimnames = dimnames(cnr$X))
    
    mm[cnr$X == deletion.threshold] = "D"

    mm[cnr$X == loss.threshold] = "L"

    mm[cnr$X == base.ploidy] = "N"

    mm[cnr$X >= gains.thresholds[1] & cnr$X <=  gains.thresholds[2]] = "G"

    if(!hyper.amplification) {
        mm[cnr$X >= amplification.thresholds[1]] <- "A"
    } else {
    
        mm[cnr$X >= amplification.thresholds[1] &
           cnr$X <= amplification.thresholds[2]] = "A"
        
        mm[cnr$X > amplification.thresholds[2]] = "A2"

    }

    if(cbioportal) {
        mm[mm == "D"] <- -2
        mm[mm == "L"] <- -1
        mm[mm == "N"] <-  0
        mm[mm == "G"] <-  1
        mm[mm == "A"] <-  2
        mm[mm == "A2"] <- 2
    }

    return(mm)
    
} ## end callCN


#' call copy numbers out of CNR data
#'
#' @param cnr a cnr bundle
#' 
#' @param cbioportal logical, weather to use the cBioPortal copy number
#' call categories of -2, -1, 0, 1, 2 to denote deletion, loss, neutral,
#' gain, amplificiation, respectively. Default is FALSE.  If TRUE,
#' hyper amplification calls will be joined with 2
#' 
#' @param deletion.threshold integer to use as locus deletion, default -0.5
#' 
#' @param loss.threshold copy number ratio to use as copy number loss,
#' default -0.1
#' 
#' @param gains.threshold copy number ratio to use as copy number gains,
#'  default 0.1
#'
#' @param amplification.threshold copy number ratio to use as copy number
#' amplification default 1
#'
#' @return
#' returns a matrix with copy number calls in the form of D, L, N, G, A1 and A2
#' from an X matrix
#'
#' | CN        |  Call  | cBioPortal Notation |
#' |----------:|: ---- :|: ----------------- :|
#' |     -0.50 |   D    |  -2                 |
#' |     -0.10 |   L    |  -1                 |
#' |-0.1 - 0.1 |   N    |  0                  |
#' | 0.1 - 1.0 |   G    |  1                  |
#' |      >1.0 |   A    |  2                  |
#'
#' @examples
#' data(cnr)
#'
#' Xc <- callDNA_CN(cnr)
#' 
#' Xp <- callDNA_CN(cnr, cbioportal = TRUE)
#'
#' @importFrom assertthat assert_that
#' @export
callDNA_CN <- function(cnr,
                       cbioportal = FALSE,
                       deletion.threshold = -0.5,
                       loss.threshold = -0.1,
                       gains.threshold = 0.1,
                       amplification.thresholds = 1.0) {

    assertthat::assert_that(loss.threshold > deletion.threshold)
    assertthat::assert_that(gains.threshold > loss.threshold)
    assertthat::assert_that(amplification.thresholds > gains.threshold)
    
    mm = matrix(NA, nrow= nrow(cnr$X), ncol = ncol(cnr$X), dimnames = dimnames(cnr$X))
    mm[cnr$X <= deletion.threshold] = "D"
    mm[cnr$X > deletion.threshold & cnr$X <= loss.threshold] = "L"
    mm[cnr$X > loss.threshold & cnr$X < gains.threshold ] = "N"
    mm[cnr$X >= gains.threshold & cnr$X < amplification.thresholds ] = "G"
    mm[cnr$X >= amplification.thresholds] = "A"

    if(cbioportal) {
            mm[mm == "D"] <- -2
            mm[mm == "L"] <- -1
            mm[mm == "N"] <-  0
            mm[mm == "G"] <-  1
            mm[mm == "A"] <-  2
    }
        
    return(mm)
    
}


