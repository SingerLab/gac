#' plot cluster optimization
#'
#' @param optC output of optClust()
#'
#' @param type type of plot e.g point, line, both (p, l, b)
#' 
#' @param ... additional arguments passed to plot and points
#'
#' @examples
#'
#' data(cnr)
#'
#' cnr <- phylo_cnr(cnr)
#'
#' ( mopc <- optClust(cnr, opt.range = seq(0.01, 0.2, by = 0.005)) )
#'
#' plotCL(optC = mopc)
#' 
#' @importFrom graphics plot points legend
#' 
#' @export
plotCL <- function(optC, type = "b", ...) {
    graphics::plot(rownames(optC), optC[,1], pch = 1, ylim = c(0, max(optC)),
         type = type, ...)
    graphics::points(rownames(optC), optC[,2], pch = 2, type = type,  ...)
    graphics::legend("topright", legend = c("One-cell", "Multi-cell"), pch = 1:2)
}
