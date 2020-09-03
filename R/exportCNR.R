#' export CNR to flat text files
#'
#' @param cnr a cnr
#'
#' @param outdir path to output directory
#'
#' @param ... arguments passed to write.table and ape::write.tree
#'
#' @import ape
#' @import usethis
#'
#' @return
#'
#' Writes X, Y, genes, qc, gene.index, and chromInfo to `outdir`.  If
#' `phylo` exists, a `.newick` format tree is also exported. Furthermore,
#' if `DDRC` profiles are also present, files for these are also written
#' to the specified directory.
#' 
#' 
#' @examples
#'
#' data(cnr)
#' 
#' export(cnr, outdir = "cnr_out")
#' 
#' @export
exportCNR <- function(cnr, outdir, ...) {

    usethis::use_directory(outdir)
    
    write.table(cnr[["X"]], file = file.path(out.dir, "X.txt"), row.names = TRUE,
                col.names = TRUE, sep = "\t", quote = FALSE, ...)

    write.table(cnr[["Y"]], file = file.path(out.dir, "Y.txt"), row.names = TRUE,
                col.names = TRUE, sep = "\t", quote = FALSE, ...)

    write.table(cnr[["genes"]], file = file.path(out.dir, "geneCN.txt"),
                row.names = TRUE,col.names = TRUE,
                sep = "\t", quote = FALSE, ...)

    write.table(cnr[["qc"]], file = file.path(out.dir, "qc.txt"),
                row.names = TRUE, col.names = TRUE,
                sep = "\t", quote = FALSE, ...)
                
    write.table(cnr[["chromInfo"]], file = file.path(out.dir, "chromInfo.txt"), 
                row.names = TRUE, col.names = TRUE,
                sep = "\t", quote = FALSE, ...)

    write.table(cnr[["gene.index"]], file = file.path(out.dir, "gene.index.txt"),
                row.names = TRUE, col.names = TRUE,
                sep = "\t", quote = FALSE, ...)
    
    if(!is.null(cnr[["phylo"]])) {
        ape::write.tree(cnr[["phylo"]],
                        file = file.path(out.dir, "phylo.newick"), ...)
    }

    if(!is.null(cnr[["DDRC.df"]])) {
        
    write.table(cnr[["DDRC.df"]], file = file.path(out.dir, "DDRC.txt"),
                row.names = TRUE, col.names = TRUE,
                sep = "\t", quote = FALSE, ...)
    
    write.table(cnr[["DDRC.g"]], file = file.path(out.dir, "DDRC_genes.txt"),
                row.names = TRUE, col.names = TRUE,
                sep = "\t", quote = FALSE, ...)
    }
    
} # end exportCNR
    
