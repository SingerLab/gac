% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pull_gene_copy_numbers.R
\name{pull_gene_copy_numbers}
\alias{pull_gene_copy_numbers}
\title{pull a set of copy numbers from a cnr into a data.frame}
\usage{
pull_gene_copy_numbers(cnr, genes)
}
\arguments{
\item{cnr}{a cnr bundle}

\item{genes}{character, genes of interst}
}
\value{
A `data.frame` containing the copy numbers for selected genes.  The data set
 is restricted to the genes present in the `genes` matrix.  However, those not
 found will be listed out.
}
\description{
pull a set of copy numbers from a cnr into a data.frame
}
\examples{
data(cnr)

pull_gene_copy_numbers(cnr, c("MDM2", "CDK4"))

genes.of.interest <- list_gene_symbols(cnr, at = "12:58140000:58200000")

pull_gene_copy_numbers(cnr, genes.of.interest)

}
