% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gac.R
\docType{package}
\name{gac}
\alias{gac-package}
\alias{gac}
\title{gac}
\source{
\url{https://github.com/SingerLab/gac}
}
\description{
GAC is an open-source framework and toolkit for management and analysis of
genome wide copy number data.  It provides a set of standard matrices that
are symultaneously manipulated by various functions to maintain syncronized
bin, gene, annotation, and qc information from a single-cell DNA sequencing
experiment. Efforts will be made to add support for same-cell G+T and
log-ratio data for the analysis of bulk DNA and RNA data from the same 
samples.
}
\details{
GAC is actively being developed and is currently in ALPHA pre-release.
My hope is that others will find this helpful for their copy number needs,
and share new methods.

The HeatmapCNR function in GAC implements Bray-Curtis disimilarity on
ComplexHeatmap as the default method to estimate cell-cell distances and
generate the heatmaps.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/SingerLab/gac/}
  \item Report bugs at \url{https://github.com/SingerLab/gac/issues}
}

}
\author{
Rodrigo Gularte Merida \email{gularter@mskcc.org}
}
\keyword{internal}
