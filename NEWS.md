# gac 0.0.9037
* improved functions to call copy number states from bulk DNA or Integer Copy Number

* calls can be exported as text (deletion, loss, gain, amplification)  or with cbioportal notation, -2, -1, 0, +1, +2

* calls include an additional category representing > 20 copies of a gene.  Term is high_amplification, and in cbioportal is labeled as +4

* added functions to estimate alteration frequencies in bulk DNA or Integer Copy Number

* added option to plot bins or gene alteration frequences

* minor changes for compatibility

* bug fix in sync_cnr


# gac 0.0.9036
* allow arguments to be passed down tho HeatmapAnnotation when implementing create_chromosome_annotation
* simplified create_chromosome_annotation base functions

* `add_in_silico_root` now exits if "X" contains `cell.name` e.g. 'diploid'

* adds plot_sK to inspect kParameter and kStable

* new method for assigning clusters in Bray Curtis Dissimilarity

# gac 0.0.9035
* added function to add an in-silico root cell and clone

* bin.id checks between chrmoInfo and gene.index

* standardizing function names, renamed as export_cnr and summary_cnr

* pull_gene_details now pulls information from a gene list, rather than regions

# gac 0.0.9034
* added functions to order bins and genes, and are implemented by default in sync_cnr

* added default colors to vdjHeatmap

* modified names of functions related to gene and region information lookup

* getting_started vignette cleanup

* copynumbers example data now matches that of the example cnr

# gac 0.0.9033
* added features to allow FastME alogrithms for  phylogenetic cell and clone pseudobulk reconstruction.  Implemented via R/ape

# gac 0.0.9032
* new sync_cnr function to resyncronize cnr tables after merge, etc
  This is now applied by default after running addCells, addPheno, and buildCNR

* adding default colors to heatmap plots

* debuged DDRC.g creation as there was a bug introduced and went missing

* added function to append information to the gene.index

* removing color map objects 

# gac 0.0.9031
* debug pull_gene_details and get_gene_details

# gac 0.0.9030
* new `split_cnr` function to convert a cnr object into a list of cnr split by a categorical variable
 e.g. sample ID, treated vs control, etc.

# gac 0.0.9029
* added option to change gene.type column in genotype_vdj

* changed `seqnames` column to `chrom` in gene index and matches chromInfo

* debugged other functions after setting gene.type.column

* cnr now uses default bioMart `gene_biotype` column name

# gac 0.0.9028
* added code to estimate q-values in histology comparisons

# gac 0.0.9027
* introduction of gene lookup functions based on coordinates

* further debugging to pass CRAN

# gac 0.0.9026

* added functions to associate copy number to categorical phenotypes

* added functions to pull copy number froms specific genes

* added functions to visualize association resutls, e.g. manhattan plots, and effect plots

* added functions to convert from .seg Run Length Encoding (RLE) data to bin coordinates

* added functions to facilitate creating chromsome annotations for custom plots

* other helper functions

* removed @export tag in DepMap, OncoKB, and GISTIC2 specific functions

* tagged internal functions 

* use of vegan::diversity rather than entropy::entropy

# gac 0.0.9025

# gac 0.0.9024

* added phyloDDRC among other functions

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
