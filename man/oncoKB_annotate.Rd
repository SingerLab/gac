% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oncoKB_annotate.R
\name{oncoKB_annotate}
\alias{oncoKB_annotate}
\title{annotate gene index with oncoKB}
\usage{
oncoKB_annotate(cnr, oncokb)
}
\arguments{
\item{cnr}{a cnr bundle}

\item{oncokb}{an oncokb table}
}
\value{
Function retunrs an annotated gene.index with OncoKB cancer.gene,oncogene, and
 tummor suppressor gene (tsg).  Currently does not bring in the Actionability
 OncoKB level annotations.
}
\description{
Cross-references a gene index with OncoKB.  Information about OncoKB is available
 at https://www.oncokb.org
}
\details{
When using OncoKB, please cite: [Chakravarty et al., JCO PO 2017.](https://ascopubs.org/doi/full/10.1200/PO.17.00011)
}
\examples{

\dontrun{
oncoKB <- read.delim("inst/extdata/oncokb.txt")

cnr <- oncoKB_annotate(cnr, oncokb = oncoKB)

}

}
