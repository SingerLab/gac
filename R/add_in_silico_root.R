
#' build an in-silico root cell and clone profile
#'
#' By default cell is a female diploid cell
#'
#' @param cnr a cnr
#'
#' @param cell.name name of root cell. default "diploid"
#' 
#' @param female logical, weather to build a female or male
#'  cell and clone, default is female. If female is NULL, no sex chromosomes
#' correction is performed, all bins will have copy number equal base.ploidy
#'
#' @param base.ploidy integer, base ploidy (2, 4, 6, 8), default 2, for diploid cell
#'
#' @return
#'
#' Append an in-silico diplod (or polyploid) cell to the cnr.  By default a female
#' human genome is created.  If FALSE, a male human genome is constructed.  NULL will
#' not peform a gender correction and the entire genome will be set to the base ploidy.
#'
#' In tetraploid males copy number is half of the base.ploidy e.g. for 4n : X = 2, Y = 2.
#' 
#' @examples
#'
#' data(cnr)
#'
#' cnr <- add_in_silico_root(cnr)
#'
#' head(cnr$Y)
#'
#' @importFrom assertthat assert_that
#' @export
add_in_silico_root <- function(cnr, cell.name = "diploid",
                               female = TRUE, base.ploidy = 2L) {
    ## 
    if(!is.integer(base.ploidy)) {
        base.ploidy <- as.integer(round(base.ploidy))
    }
    ## 
    assertthat::assert_that(nrow(cnr$X) == nrow(cnr$chromInfo))
    assertthat::assert_that(ncol(cnr$genes) == nrow(cnr$gene.index))
    ## 
    dX <- data.frame(rep(base.ploidy, times = nrow(cnr$chromInfo)))
    names(dX) <- cell.name
    dG <- data.frame(rep(base.ploidy, times = nrow(cnr$gene.index)))
    names(dG) <- cell.name
    ##
    if(!is.null(female)) {
        if(female) {
            dX[cnr$chromInfo$bin.chrom == "Y", cell.name] <- 0
            dG[cnr$gene.index$chrom == "Y", cell.name] <- 0
        } else {
            dX[cnr$chromInfo$bin.chrom %in% c("X", "Y"), cell.name] <- base.ploidy / 2
            dG[cnr$gene.index$chrom %in% c("X", "Y"), cell.name] <- base.ploidy / 2
        }
    }
    ##
    if(! cell.name %in% colnames(cnr$X)) {
        cnr$X <- cbind(cnr$X, dX)
    } else {
        warning(paste(cell.name, "exists in X"))
    }
    if(! cell.name %in% colnames(cnr$genes)) {
        cnr$genes[cell.name, ] <- as.integer(dG[,1])
    } else {
        warning(paste(cell.name, "exists in genes"))
    }
    ##
    if(length(cnr$DDRC.df)) {
        if(! cell.name %in% colnames(cnr$DDRC.df)) {
            cnr$DDRC.df <- cbind(cnr$DDRC.df, dX)
        } else {
            warning(paste(cell.name, "exists in DDRC.df"))
        }
    }
    ##    
    if(length(cnr$DDRC.g)) {
        if(! cell.name %in% colnames(cnr$DDRC.df)) {
            cnr$DDRC.g <- cbind(cnr$DDRC.g, dG)
        } else {
            warning(paste(cell.name, "exists in DDRC.g"))
        }
    }
    ##
    cnr$Y[cell.name, ] <- NA
    cnr$Y[cell.name, "cellID"] <- cell.name
    cnr$qc[cell.name, "cellID"] <- cell.name
    cnr$qc[cell.name, "qc.status"] <- "PASS"

    cnr$cells <- colnames(cnr$X)
    
    return(cnr)
}

