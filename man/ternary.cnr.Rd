% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ternary.cnr.R
\name{ternary.cnr}
\alias{ternary.cnr}
\title{build ternary matrix from integer copy number data}
\usage{
ternary.cnr(X, gain = 3, amp = 20)
}
\arguments{
\item{X}{a X matrix either the bins or genes}

\item{gain}{an integer copy number value specifying the minimum number of copies to call a gain as 1}

\item{amp}{an integer value to specify the minimum number of copies to call as a high level amplification e.g. double minute}
}
\description{
This function builds a ternary matrix from the cnr$X or cnr$genes.  It was designed as a helper function to generate the input for infSCITE.
}
\details{
Because infSCITE was developed for mutation data, this function attempts to recreate the scenario where 0 is no alteration, 1 is a low level / impact alteration, and 2 is high level / impact alteration.

By default, a deletion i.e. copy number of 0 will be treated as a high level deletion, and a loss i.e. one copy as a 1.  The diploid state is treated as 0 = No mutation
}
\examples{

data(cnr)

G <- ternary.cnr(cnr$genes[, c("CDK4", "MDM2")])

head(G)

}
