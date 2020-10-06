#' plot cluster optimization
#'
#' @param optC output of optClust()
#'
#' @param type type of plot e.g point, line, both (p, l, b)
#' 
#' @param ... additional arguments passed to plot and points
#'
#'
#' @import graphics
#' 
#' @examples
#'
#' data(cnr)
#'
#' cnr <- phyloCNR(cnr)
#'
#' ( mopc <- optClust(cnr, opt.range = seq(0.01, 0.2, by = 0.005)) )
#'
#' plotCL(optC = mopc)
#' 
#' @export
plotCL <- function(optC, type = "b", ...) {
    plot(rownames(optC), optC[,1], pch = 1, ylim = c(0, max(optC)),
         type = type, ...)
    graphics::points(rownames(optC), optC[,2], pch = 2, type = type,  ...)
    graphics::legend("topright", legend = c("One-cell", "Multi-cell"), pch = 1:2)
}
