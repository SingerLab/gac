#' toSignac: the open-source simple integrated analysis of cells/copynumber
#'
#'
#' toSignac is an open-source framework and toolkit for managing genome wide copy
#' number data.  It provides a set of standard matrices that are symultaneously
#' manipulated by various functions to maintain syncronized bin, gene, annotation,
#' and qc information from a single-cell DNA sequencing experiment as same-cell
#' G+T.  Efforts will be made to add support for log-ratios for the analysis of
#' bulk DNA data.
#'
#' toSignac is actively being developed and is currently in ALPHA pre-release.
#' My hope is that others will find this helpful for their copy number needs,
#' and share new methods.
#'
#' The HeatmapCNR function in toSignac implements Bray-Curtis disimilarity on
#' ComplexHeatmap as the default method to estimate cell-cell distances and
#' generate the heatmaps.
#'
#'
#' @docType package
#' @author Rodrigo Gularte Merida \email{gularter@mskcc.org}
#' @keywords package
#' @source \url{https://github.com/SingerLab/toSignac}
#' @name toSignac
#'
#'
#' @import ComplexHeatmap
#' @import GenomicRanges
#' @import circlize
#' @import colourvalues
#' @import vegan
#' 
NULL
