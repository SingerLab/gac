% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setClusters.R
\name{setKcc}
\alias{setKcc}
\title{Set cluster membership for K clusters}
\usage{
setKcc(cnr, kCC = NULL, prefix = "X")
}
\arguments{
\item{cnr}{a cnr object}

\item{kCC}{number of K clusters to use, 
select optimum or maximum k from kSpectral.  Default is NULL
which will pull the optK[kCC] from kStats}

\item{prefix}{prefix charcter to append to Consensus Clusters}
}
\value{
returns cluster membership based on consensus clustering for a
 specified kCC
}
\description{
Set cluster membership for K clusters
}
