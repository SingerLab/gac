% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cluster_profiles.R
\name{get_cluster_profiles}
\alias{get_cluster_profiles}
\title{Build a consensus copy number profile for each cluster or clone}
\usage{
get_cluster_profiles(
  cnr,
  minimum_cells = 3,
  base.ploidy = 2,
  cluster.column = "final_cluster"
)
}
\arguments{
\item{cnr}{a cnr bundle containing a `final_cluster` column in Y.}

\item{minimum_cells}{minimum number of cells in a cluster, best if greater than 3
to estimate a median}

\item{base.ploidy}{base ploidy of the tumor, default 2 i.e. diploid}

\item{cluster.column}{column containing clusters}
}
\value{
Function returns the cnr with three additional tables. 
* uclust : number of cells in each final_cluster, only clusters
 greater than the specified minimum number of cells is shown
* DDRC.df: a matrix containing the representative profile for each
 `final_cluster` at the bin level.
* DDRC.g : interpolation of the DDRC.df at the gene level for each `final_cluster`.
}
\description{
Build a consensus copy number profile for each cluster or clone
}
\examples{

data(cnr)

cnr <- phylo_cnr(cnr)

cnr <- setBrayClusters(cnr, tree.height = 0.065)

cnr <- cluster_heterogeneity(cnr, cluster_column = "BrayC")

cnr <- get_cluster_profiles(cnr)

}
