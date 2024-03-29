% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phylo_cnr.R
\name{phylo_cnr}
\alias{phylo_cnr}
\title{Calculating cell-to-cell distances, heirarchical clustering, and generating a tree class `phylo`}
\usage{
phylo_cnr(
  cnr,
  root.cell = NULL,
  dist.method = "bray",
  hclust.method = "ward.D2",
  tree.method = "bal",
  ...
)
}
\arguments{
\item{cnr}{a cnr bundle}

\item{root.cell}{a cellID to root the three.}

\item{dist.method}{method for calculating cell-to-cell distance
(see \link[vegan]{vegdist})}

\item{hclust.method}{method for heirarchical clustering (see hclust)}

\item{tree.method}{minimum evolution phylogenetics method, can be
`bal`, `ols` or NULL. Default is `bal`.}

\item{...}{other parameters passed to \link[vegan]{vegdist}}
}
\value{
Creates a cell-to-cell distance matrix, runs heirarchical clustering.
 By defaul cell phylogenetics is infered by \link[ape]{fastme.bal},
 and alternatively by `fastme.ols`. When `tree.method = NULL, the
`hclust` object is converted  to an `ape` class `phylo` object to
 represent the cell phylogeny.

\itemize{
  \item cdb cell to cell Bray-Curtis dissimiarly 
  \item hcdb heirarchical clustering of distance matrix
  \item phylo ape class `phylo` object
}
}
\description{
Calculating cell-to-cell distances, heirarchical clustering, and generating a tree class `phylo`
}
\examples{

data(cnr)

## unrooted cell phylogenetic tree
cnr <- phylo_cnr(cnr)

cnr$phylo
plot(cnr$phylo)


## rooted cell phylogenetic tree
cnr <- phylo_cnr(cnr, root.cell = "cell0")

cnr$phylo
plot(cnr$phylo)

}
\references{
Vincent Lefort,  Richard Desper,  Olivier Gascue. 2015. "FastME 2.0:
  A Comprehensive, Accurate, and Fast Distance-Based Phylogeny Inference
  Program". Molecular Biology and Evolution, Volume 32, Issue 10, October
  2015, Pages 2798–2800. <https://doi.org/10.1093/molbev/msv150?>

Paradis E, Schliep K (2019). “ape 5.0: an environment for modern
  phylogenetics and evolutionary analyses in R.” Bioinformatics,
  *35*, 526-528.  <https://doi.org/10.1093/bioinformatics/bty633>.

Oksanen J, Simpson G, Blanchet F, Kindt R, Legendre P, Minchin P,
  O'Hara R, Solymos P, Stevens M, Szoecs E, Wagner H, Barbour M,
  Bedward M, Bolker B, Borcard D, Carvalho G, Chirico M, De Caceres M,
  Durand S, Evangelista H, FitzJohn R, Friendly M, Furneaux B, Hannigan
  G, Hill M, Lahti L, McGlinn D, Ouellette M, Ribeiro Cunha E, Smith T,
  Stier A, Ter Braak C, Weedon J (2022). _vegan: Community Ecology
  Package_. R package version 2.6-4,
  <https://CRAN.R-project.org/package=vegan>.
}
