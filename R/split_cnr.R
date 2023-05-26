#' split a cnr by one variable
#' @param cnr a cnr object
#' @param split.by name of the column to split. Must be a column in
#'  the phenotype or qc tables
#'
#' @return
#' Function subsets based on `split.by` and returns a list of named cnr objects.
#'  Element names are taken from the unique values of split.by
#' 
#' @importFrom assertthat assert_that
#'
#' @examples
#'
#' data(cnr)
#'
#' cnrL <- split_cnr(cnr, split.by = "category1")
#'
#' lapply(cnrL, summary_cnr)
#' 
#' @export
split_cnr <- function(cnr, split.by) {

    assertthat::assert_that(split.by %in% union(names(cnr$Y), names(cnr$qc)))
    
    
    if(split.by %in% names(cnr$Y)) {

        ubio <- unique(cnr$Y[, split.by])

        cell.lists <- lapply(ubio, function(i) {
            cc <- cnr$Y[cnr$Y[, split.by] == i, "cellID"]
        })

    } else {

        if(split.by %in% names(cnr$qc)) {

            ubio <- unique(cnr$qc[, split.by])

            cell.lists <- lapply(ubio, function(i) {
                cc <- cnr$qc[cnr$qc[, split.by] == i, "cellID"]
            })

        }
    }

    out <- lapply(cell.lists, function(i) keepCells(cnr, keep = i))
    names(out) <- ubio

    return(out)
} # end split_cnr

