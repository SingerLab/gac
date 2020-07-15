#' optimizing clustering
#'
#' @param cnr a cnr bundle
#'
#' @param opt.range
#'
#' @export
optClust <- function(cnr, opt.range = seq(0, 0.6, by = 0.05)) {

    clL <- sapply(opt.range, function(h) {
        ctbl <- sort(table(cutree(cnr[["hcdb"]], h = h)),
                     decreasing = TRUE)
    })
    names(clL) <- opt.range

    omt <- matrix(t(sapply(clL, function(i) {
        cbind(sum(i == 1), sum(i != 1), sum(i != 1)/length(i))
    })),
    ncol = 3,
    dimnames = list(opt.range, c("One-cell", "Multi-cell", "%CMC")))
    
    omt
    
} # optClust
