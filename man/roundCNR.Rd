% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/roundCNR.R
\name{roundCNR}
\alias{roundCNR}
\title{rounded quantal matrix}
\usage{
roundCNR(X, neut = 2.5, loss = 1.2, del = 0.2)
}
\arguments{
\item{X}{a numerical matrix composed of bins with copy numbers estimates integer or numerical}

\item{neut}{maximum value that will be floored to a neutral copy number for diploid species i.e. 2}

\item{loss}{lower bound for neutral, and upper bound threshold to be considered a loss of one copy}

\item{del}{maximum threshold to consider a total deletion of the bin i.e. 0 copies}
}
\value{
Rounds X based on quantal
}
\description{
Rounds an X matrix to the nearest integer with special
threasholds specific to single-cell quantal data from the Varbin algorithm
for deletions, losses, and neutral (2 copies)
}
\examples{

data(copynumbers)

head(copynumbers)

cni <- roundCNR(copynumbers)

head(cni)

}
