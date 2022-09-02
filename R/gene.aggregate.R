#' combine genes w/equal frequency to be only once
#'
#' This is a helper function to process the binary or ternary matrix.  Because
#' gene data is interpolated from the bin data, linked loci within a bin will
#' be duplicated.  Duplicate frequencies create infinite combinations in the
#' trees in infSCITE.
#'
#' This function checks if two rows have are identical and merges them into
#' one. It also creates unique rownames to know what was de-duplicated
#'
#' @param Z binary incidence matrix or G genotype matrix.  Internally, both
#' Z or G will be turned into a Z matrics for frequency calculations.
#'
#' @return
#' returns a de-duplicated incidence or ternary matrix with rownams naving
#' the multipe duplicated rows
#'
#' 
#' @examples \dontrun{
#' data(cnr)
#'
#' ## for binary data
#' Z <- binary.cnr(cnr$genes[, c("CDK4", "MDM2", "HMGA2")])
#'
#' aggrMat.Z <- gene.aggregate(Z = Z)
#'
#' ## for ternary data
#' G <- ternary.cnr(cnr$genes[, c("CDK4", "MDM2", "HMGA2")])
#'
#' aggrMat.G <- gene.aggregate(Z = G)
#' }
#' 
#' @keywords internal
#' @noRd
gene.aggregate <- function(Z) {

    fq <- rowSums(binary.cnr(Z))/ncol(binary.cnr(Z))
    dups <- fq[duplicated(Z) | duplicated(Z, fromLast = TRUE)]
    
    dgroups <- sapply(dups, function(i) names(fq)[fq == i])
    duse <- sapply(dgroups, `[`, 1)
    dnew <- as.character(sapply(dgroups, function(i) paste(dput(i), collapse = "_")))
    dd <- data.frame(duse = as.character(duse), dnew = as.character(dnew))
    dd

} 

