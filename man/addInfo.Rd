% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addInfo.R
\name{addInfo}
\alias{addInfo}
\title{Adds info to the chromInfo}
\usage{
addInfo(cnr, df)
}
\arguments{
\item{cnr}{a cnr bundle}

\item{df}{a data.frame with the data to incorporate.  Particularly
useful for adding p-values, genetic effects, etc to the bins}
}
\value{
Returns a CNR object with added columns to the chromInfo. e.g. p-values for genome wide scans
}
\description{
The function adds information at the bin level.  The gene information is an
interpolation of the bin.  The number of independent tests you are performing
is the number of common, non-identical segments in the data.   However, to keep
the framework simple, it's best to append results from algorithms to the
chromInfo. E.g. if .X is 5000 bins, then your genome wide algorithim
is 5000 p-values.  Here is where you add it after generating it.
}
\examples{

data(cnr)

fakePval <- data.frame(pval = runif(5000))

cnr <- addInfo(cnr, df = fakePval)

head(cnr$chromInfo)

}
