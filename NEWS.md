# gac 0.0.9023

* minor fix to pass `roundCNR` arguments from `buildCNR`

# gac 0.0.9022

* added clone assignment ranked by FGA and/or frequency

* added ability to root the cell phylogenetic tree

* added function to run clone phylogenetic analysis

# gac 0.0.9021

* added VDJ clustering with Bray-Curtis Dissimilarity

* fixed vdj cell type in keepCells and excludeCells

# gac 0.0.9020

* added support for subbsetting based on bins, genes, full chromosomes and genomic coordinates

* debugged build notes
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

# gac 0.0.9019

* added plot_frequencies to plot amplification and deletion frequencies

* in setBrayClusters, if tree.height is NULL,  the minimum intersect between one-cell and multi-cell clusters is used by default as the tree.height

# gac 0.0.9018

* section on preparing your data added to `getting_started.Rmd`

* minor debugging on warnings, and unit tests

* moved badge to maturing

# gac 0.0.9017

* Added summaryCNR function

* debugged warnings and notes

# gac 0.0.9016

* Added VDJ genotyping and visualization

* Progress status in consensusClusterCNR

# gac 0.0.9013

* Added f(x) to pull information from the gene.index

# gac 0.0.9013

* Added estimation of alteration frequencies at gene and bin level

* Deduplicated example bin to gene tracks

* Added DATASET.R to show how example data is generated

# gac 0.0.9012

* Added chromosome tickmarks in HeatmapCNR

* Added support for consensus clustering

* Added support for K Spectral to select max number of stable clusters from consensus clustering

# gac 0.0.9000

* Added a `NEWS.md` file to track changes to the package.
