% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_lookups.R
\name{pull_gene_details}
\alias{pull_gene_details}
\title{Pull gene details for a set of genes}
\usage{
pull_gene_details(
  cnr,
  genes = c("MDM2", "CDK4"),
  show.columns = NULL,
  identifier = "hgnc.symbol"
)
}
\arguments{
\item{cnr}{a cnr bundle}

\item{genes}{a list of genes}

\item{show.columns}{columns of gene.index to show}

\item{identifier}{gene identifier hgnc.symbol or ensembl_gene_id. default hgnc.symbol}
}
\value{
Returns the subset of the `gene.index` table for the genomic region.
}
\description{
This function subsets the gene index for a given set of genes
}
\examples{

data(cnr)

pull_gene_details(cnr)

pull_gene_details(cnr,
  genes = c("JUN", "MDM2", "CDK4"),
  show.columns = c("hgnc.symbol", "bin.id", "gene_biotype"))

pull_gene_details(cnr,
  genes = c("ENSG00000177606", "ENSG00000135446", "ENSG00000135679"),
  identifier = "ensembl_gene_id",
  show.columns = c("hgnc.symbol", "bin.id", "gene_biotype"))


}
