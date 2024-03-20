#' Copy number call categories from ratio data
#'
#' @param cnr a cnr bundle
#'
#' @param base.ploidy sample base ploidy e.g. 2, 4, 8.  Default is 2
#'
#' @param cbioportal logical, weather to use the cBioPortal copy number
#' call categories of -2, -1, 0, +1, +2 to denote deletion, loss, neutral or
#' base ploidy, gain, amplificiation, respectively. Default is FALSE.
#' We added one extra category `+4` representing high amplifications equivalent
#' to a locus having with >20 copies (log2 Ratio > 2.16). 
#' 
#'
#' @param ... Optional parameters passed down to
#'   \code{cbioportal_states_from_log2_ratio} or
#'   \code{cbioportal_states_from_integer_cn}.  These are  deletion.threshold,
#'   loss.threshold, gains.threshold, amplification.thresholds.
#'
#' amplification.thresholds can have one value e.g. 9.  In this case, all
#'   copy numbers >= to 9 will be considered amplifications, and will skip
#'   the high amplification call.  By default we pass two values, c(9, 20) for
#'   integer copy number, and 0.8 and 2.16 for log2 Ratio data common in
#'   bulk DNA.
#'
#' @return
#' Function returns a a list with two matrices with copy number
#' calls, one for bins, and one for genes.   Calls are interpreted as
#' deletion, loss, (ploidy)N, gain, amplification, high_amplicifation,
#' from an genotype X matrix.  Default thresholds are described below.
#'
#' Thresholds are chosen based on the cnr$bulk argument.  If bulk = FALSE,
#'   integer copy number thresholds are applied.  If bulk = FALSE,
#'   log2 ratio thresholds are applied.  
#' 
#' For integer copy number data thresholds are:
#'
#' |   CN | Call               | cBioPortal Notation |
#' |------+--------------------+---------------------|
#' |    0 | deletion           |                  -2 |
#' |    1 | loss               |                  -1 |
#' |    2 | (ploidy)N          |                   0 |
#' |  3-8 | gain               |                   1 |
#' | 9-20 | amplification      |                   2 |
#' |  >20 | high_amplicifation |                   4 |
#'
#' 
#' For bulk DNA log2 Ratio data thresholds are:
#' 
#' | log2(Ratio) | Call               | cBioPortal Notation |
#' |-------------+--------------------+---------------------|
#' |        -1.2 | deletion           |                  -2 |
#' |        -0.6 | loss               |                  -1 |
#' |  -0.6, +0.4 | (ploidy)N          |                   0 |
#' |     0.4-0.8 | gain               |                   1 |
#' |    0.8-2.16 | amplification      |                   2 |
#' |       >2.16 | high_amplicifation |                   4 |
#'
#' @references
#' Thresholds for DNA were chosen from Mina et al. (2020) PMID:32989323
#' 
#' 
#' @examples
#' data(cnr)
#'
#' Xc <- cbioportal_copy_number_states(cnr)
#' 
#' Xp <- cbioportal_copy_number_states(cnr, cbioportal = TRUE)
#'
#' @importFrom assertthat assert_that
#' 
#' @export
cbioportal_copy_number_states <-  function(cnr, base.ploidy =  2,
                                           cbioportal =  FALSE, ...) {
    
    if(cnr$bulk) {
        bx <- cbioportal_states_from_log2_ratio(cnr,
                                                base.ploidy =  base.ploidy,
                                                cbioportal = cbioportal,
                                                ...)
        gx <- cbioportal_states_from_log2_ratio(cnr,
                                                genes =  TRUE,
                                                base.ploidy =  base.ploidy,
                                                cbioportal = cbioportal,
                                                ...)
        
    } else {
        
        bx <- cbioportal_states_from_integer_cn(cnr,
                                                base.ploidy =  base.ploidy,
                                                cbioportal = cbioportal,
                                                ...)
        gx <- cbioportal_states_from_integer_cn(cnr,
                                                genes = TRUE,
                                                base.ploidy =  base.ploidy,
                                                cbioportal = cbioportal,
                                                ...)
    }
    
    out <- list("bins" =  bx,
                "genes" =  gx)
    
    return(out)
} ## end cbioportal_copy_number_states



#' Copy number call categories from ratio data
#'
#' @param cnr a cnr bundle
#'
#' @param genes If TRUE calls are infered on `genes` data, default is FALSE
#'
#' @param base.ploidy integer to denote base ploidy to use as Neutral
#' 
#' @param cbioportal logical, weather to use the cBioPortal copy number
#' call categories of -2, -1, 0, 1, 2 to denote deletion, loss, neutral,
#' gain, amplificiation, respectively. Default is FALSE.   We added an additional
#'   call `4` to denote high amplification for regions with > 20 copies.
#' 
#' @param deletion.threshold integer to use as locus deletion, default = 0
#' 
#' @param loss.threshold number of copies to use as copy number loss, default = 1
#' 
#' @param gains.thresholds number of copies to use as copy number gains,
#'  default to a ower bound of 3 (+1 copy from 2N) and upper bound of 8 i.e.
#'    +3 to +6 extra copies
#'
#' @param amplification.thresholds number of copies to use as copy number
#' amplification defaults to a lower bound of 9 (+7 copies). Optional upper bound
#' to call high amplifications (default is 20)
#' 
#'
#' @return
#' Function returns a a list with two matrices with copy number
#' calls, one for bins, and one for genes.   Calls are interpreted as
#' deletion, loss, (ploidy)N, gain, amplification, high_amplicifation,
#' from an genotype X matrix.  Default thresholds are described below.
#'
#' |   CN | Call               | cBioPortal Notation |
#' |------+--------------------+---------------------|
#' |    0 | deletion           |                  -2 |
#' |    1 | loss               |                  -1 |
#' |    2 | (ploidy)N          |                   0 |
#' |  3-5 | gain               |                   1 |
#' | 6-20 | amplification      |                   2 |
#' |  >20 | high_amplicifation |                   4 |
#'
#' 
#' @importFrom assertthat assert_that
#' 
#' @keywords internal
#' @noRd
cbioportal_states_from_integer_cn <- function(cnr,
                                          genes =  FALSE,
                                          cbioportal = FALSE,
                                          base.ploidy = 2,
                                          deletion.threshold = 0,
                                          loss.threshold = 1,
                                          gains.threshold = 3,
                                          amplification.thresholds = c(9,20)) {

    ## checks
    assertthat::assert_that(loss.threshold > deletion.threshold)
    assertthat::assert_that(base.ploidy > loss.threshold)
    assertthat::assert_that(gains.threshold[1] > base.ploidy)
    assertthat::assert_that(amplification.thresholds[1] > gains.threshold)
    if(length(amplification.thresholds) == 2) {
        assertthat::assert_that(amplification.thresholds[1] <=
                                amplification.thresholds[2])
    }
    
    ## choose dataset
    if(genes) {
        X <- t(cnr$genes)
    } else {
        X <- cnr$X
    }
    
    ## convert to cbioportal code
    mm <- apply(X, 2, function(i) {
        out <- ifelse(i == deletion.threshold, -2,
               ifelse(i == loss.threshold, -1,
               ifelse(i == base.ploidy, 0,
               ifelse(i >= gains.threshold[1] &
                      i <= amplification.thresholds[1], 1, 
               ifelse(i >  amplification.thresholds[1], 2, NA)))))
        })
    
    ## optional high_amplification
    if(length(amplification.thresholds) == 2) {
        mm[X > amplification.thresholds[2]] <- 4
    }
    
    ## convert to biological interpretation
    if(!cbioportal) {
        mm <- data.frame(apply(mm, 2, function(i) {
            out <-  droplevels(
                factor(i,
                       levels =  c(-2, -1, 0, 1, 2, 4),
                       labels = c("deletion", "loss",
                                  paste0(base.ploidy, "N"),
                                  "gain", "amplification",
                                  "high_amplification"))
            )
            return(out)
        }))
    }
    
    return(mm)
    
} ## end cbioportal_states_from_integer_cn


#' Copy number call categories from ratio data
#'
#' @param cnr a cnr bundle
#'
#' @param genes run calls on genes matrix
#'
#' @param cbioportal logical, weather to use the cBioPortal copy number
#'   call categories of -2, -1, 0, 1, 2 to denote deletion, loss, neutral,
#'   gain, amplificiation, respectively. Default is FALSE. We added an additional
#'   call `4` to denote high amplification for regions with log2R > 2.16, roughly
#'   equivalent to > 20 copies.
#' 
#' @param base.ploidy integer to denote base ploidy
#' 
#' @param deletion.threshold log2 ratio to call deletion, default < -1.2
#' 
#' @param loss.threshold log2 ratio to use as copy number loss, default < -0.6
#' 
#' @param gains.thresholds number of copies to use as copy number gains,
#'   default >= 0.4 
#'
#' @param amplification.thresholds number of copies to use as copy number
#'   amplification defaults to a lower bound of 0.8. Optional upper bound
#'   to call high amplifications, default is 2.16
#' 
#'
#' @return
#' Function returns a matrix with copy number calls.   Calls are interpreted as
#'   deletion, loss, (ploidy)N, gain, amplification, high_amplicifation,
#'   from an genotype X matrix.  Default thresholds are described below.
#' 
#' | log2(Ratio) | Call               | cBioPortal Notation |
#' |-------------+--------------------+---------------------|
#' |        -1.2 | deletion           |                  -2 |
#' |        -0.6 | loss               |                  -1 |
#' |  -0.6, +0.4 | (ploidy)N          |                   0 |
#' |     0.4-0.8 | gain               |                   1 |
#' |    0.8-2.16 | amplification      |                   2 |
#' |       >2.16 | high_amplicifation |                   4 |
#' 
#' @importFrom assertthat assert_that
#' 
#' @keywords internal
#' @noRd
cbioportal_states_from_log2_ratio <- function(cnr,
                                              genes =  FALSE,
                                              cbioportal = FALSE,
                                              base.ploidy = 2,
                                              deletion.threshold = -1.2,
                                              loss.threshold = -0.6,
                                              gains.threshold = 0.4,
                                              amplification.thresholds = c(0.8, 2.16)) {
    
    ## checks
    assertthat::assert_that(loss.threshold > deletion.threshold)
    assertthat::assert_that(amplification.thresholds[1] > gains.threshold)
    if(length(amplification.thresholds) == 2) {
        assertthat::assert_that(amplification.thresholds[1] <=
                                amplification.thresholds[2])
    }
    
    ## choose dataset
    if(genes) {
        X <- t(cnr$genes)
    } else {
        X <- cnr$X
    }
    
    ## convert to cbioportal code
    mm <- apply(X, 2, function(i) {
        out <- ifelse(i <= loss.threshold, -1,
               ifelse(i <= deletion.threshold, -2,
               ifelse(i > loss.threshold &
                      i < gains.threshold, 0, 
               ifelse(i >= gains.threshold[1] &
                      i < amplification.thresholds[1], 1,
               ifelse(i >= amplification.thresholds[1], 2, NA)))))
        })

     ## optional high_amplification
    if(length(amplification.thresholds) == 2) {
        mm[X > amplification.thresholds[2]] <- 4
    }

    ## convert to biological interpretation
    if(!cbioportal) {
        mm <- data.frame(apply(mm, 2, function(i) {
            out <-  droplevels(
                factor(i,
                       levels =  c(-2, -1, 0, 1, 2, 4),
                       labels = c("deletion", "loss",
                                  paste0(base.ploidy, "N"),
                                  "gain", "amplification",
                                  "high_amplification"))
            )
            return(out)
        }))
    }
    
    return(mm)
    
} ## end cbioportal_states_from_log2_ratio

