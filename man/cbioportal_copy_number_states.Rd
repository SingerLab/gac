% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cbioportal_copy_number_states.R
\name{cbioportal_copy_number_states}
\alias{cbioportal_copy_number_states}
\title{Copy number call categories from ratio data}
\usage{
cbioportal_copy_number_states(cnr, base.ploidy = 2, cbioportal = FALSE, ...)
}
\arguments{
\item{cnr}{a cnr bundle}

\item{base.ploidy}{sample base ploidy e.g. 2, 4, 8.  Default is 2}

\item{cbioportal}{logical, weather to use the cBioPortal copy number
call categories of -2, -1, 0, +1, +2 to denote deletion, loss, neutral or
base ploidy, gain, amplificiation, respectively. Default is FALSE.
We added one extra category `+4` representing high amplifications equivalent
to a locus having with >20 copies (log2 Ratio > 2.16).}

\item{...}{Optional parameters passed down to
  \code{cbioportal_states_from_log2_ratio} or
  \code{cbioportal_states_from_integer_cn}.  These are  deletion.threshold,
  loss.threshold, gains.threshold, amplification.thresholds.

amplification.thresholds can have one value e.g. 9.  In this case, all
  copy numbers >= to 9 will be considered amplifications, and will skip
  the high amplification call.  By default we pass two values, c(9, 20) for
  integer copy number, and 0.8 and 2.16 for log2 Ratio data common in
  bulk DNA.}
}
\value{
Function returns a a list with two matrices with copy number
calls, one for bins, and one for genes.   Calls are interpreted as
deletion, loss, (ploidy)N, gain, amplification, high_amplicifation,
from an genotype X matrix.  Default thresholds are described below.

Thresholds are chosen based on the cnr$bulk argument.  If bulk = FALSE,
  integer copy number thresholds are applied.  If bulk = FALSE,
  log2 ratio thresholds are applied.  

For integer copy number data thresholds are:

|   CN | Call               | cBioPortal Notation |
|------+--------------------+---------------------|
|    0 | deletion           |                  -2 |
|    1 | loss               |                  -1 |
|    2 | (ploidy)N          |                   0 |
|  3-8 | gain               |                   1 |
| 9-20 | amplification      |                   2 |
|  >20 | high_amplicifation |                   4 |


For bulk DNA log2 Ratio data thresholds are:

| log2(Ratio) | Call               | cBioPortal Notation |
|-------------+--------------------+---------------------|
|        -1.2 | deletion           |                  -2 |
|        -0.6 | loss               |                  -1 |
|  -0.6, +0.4 | (ploidy)N          |                   0 |
|     0.4-0.8 | gain               |                   1 |
|    0.8-2.16 | amplification      |                   2 |
|       >2.16 | high_amplicifation |                   4 |
}
\description{
Copy number call categories from ratio data
}
\examples{
data(cnr)

Xc <- cbioportal_copy_number_states(cnr)

Xp <- cbioportal_copy_number_states(cnr, cbioportal = TRUE)

}
\references{
Thresholds for DNA were chosen from Mina et al. (2020) PMID:32989323
}
