% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mark.genes.R
\name{mark.genes}
\alias{mark.genes}
\title{Mark genes}
\usage{
mark.genes(cnr, gene.list, identifier = "hgnc.symbol")
}
\arguments{
\item{cnr}{the cnr bundle}

\item{gene.list}{the list of genes you wish to mark}

\item{identifier}{name of the column with gene identifiers}
}
\value{
returns a vector of named gene.list (gg) w/their rownames/bin.id
}
\description{
Produces a named list with gene name and bin.id to use in ComplexHeatmap::anno_mark
}
\examples{

data(cnr)

aa <- mark.genes(cnr, gene.list = c("CDK4", "MDM2"))

\dontrun{
geneAnno <- ComplexHeatmap::rowAnnotation(genes = anno_mark(at = aa, labels = names(aa)))
}

}
