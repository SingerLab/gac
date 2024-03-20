#' gac
#'
#' GAC is an open-source framework and toolkit for management and analysis of
#' genome wide copy number data.  It provides a set of standard matrices that
#' are symultaneously manipulated by various functions to maintain syncronized
#' bin, gene, annotation, and qc information from a single-cell DNA sequencing
#' experiment. Efforts will be made to add support for same-cell G+T and
#' log-ratio data for the analysis of bulk DNA and RNA data from the same 
#' samples.
#' 
#' GAC is actively being developed and is currently in ALPHA pre-release.
#' My hope is that others will find this helpful for their copy number needs,
#' and share new methods.
#'
#' The HeatmapCNR function in GAC implements Bray-Curtis disimilarity on
#' ComplexHeatmap as the default method to estimate cell-cell distances and
#' generate the heatmaps.
#' 
#' @author Rodrigo Gularte Merida \email{gularter@mskcc.org}
#' @keywords internal
#' @source \url{https://github.com/SingerLab/gac}
#' @name gac
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
