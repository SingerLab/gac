#' export CNR to flat text files
#'
#' @param cnr a cnr
#'
#' @param outdir path to output directory
#'
#' @param ... arguments passed to write.table and ape::write.tree
#'
#' @return
#'
#' Writes X, Y, genes, qc, gene.index, and chromInfo to `outdir`.  If
#' `phylo` exists, a `.newick` format tree is also exported. Furthermore,
#' if `DDRC` profiles are also present, files for these are also written
#' to the specified directory.
#' 
#' @examples
#'\dontrun{
#' data(cnr)
#' 
#' exportCNR(cnr, outdir = "cnr_out/")
#'}
#'
#' @importFrom utils write.table
#' @importFrom ape write.tree
#' @importFrom usethis use_directory
#' 
#' @export
exportCNR <- function(cnr, outdir = ".", ...) {

    usethis::use_directory(outdir)
    
    utils::write.table(cnr[["X"]], file = file.path(outdir, "X.txt"),
                       row.names = TRUE,
                       col.names = TRUE, sep = "\t", quote = FALSE, ...)

    utils::write.table(cnr[["Y"]], file = file.path(outdir, "Y.txt"),
                       row.names = TRUE,
                       col.names = TRUE, sep = "\t", quote = FALSE, ...)

    utils::write.table(cnr[["genes"]], file = file.path(outdir, "geneCN.txt"),
                row.names = TRUE,col.names = TRUE,
                sep = "\t", quote = FALSE, ...)

    utils::write.table(cnr[["qc"]], file = file.path(outdir, "qc.txt"),
                row.names = TRUE, col.names = TRUE,
                sep = "\t", quote = FALSE, ...)
                
    utils::write.table(cnr[["chromInfo"]], file = file.path(outdir, "chromInfo.txt"), 
                row.names = TRUE, col.names = TRUE,
                sep = "\t", quote = FALSE, ...)

    utils::write.table(cnr[["gene.index"]], file = file.path(outdir, "gene.index.txt"),
                row.names = TRUE, col.names = TRUE,
                sep = "\t", quote = FALSE, ...)
    
    if(!is.null(cnr[["phylo"]])) {
        ape::write.tree(cnr[["phylo"]],
                        file = file.path(outdir, "phylo.newick"),
                        ...)
    }

    if(!is.null(cnr[["DDRC.df"]])) {
        
        utils::write.table(cnr[["DDRC.df"]],
                           file = file.path(outdir, "DDRC.txt"),
                           row.names = TRUE, col.names = TRUE,
                           sep = "\t", quote = FALSE, ...)
        
        utils::write.table(cnr[["DDRC.g"]],
                           file = file.path(outdir, "DDRC_genes.txt"),
                           row.names = TRUE, col.names = TRUE,
                           sep = "\t", quote = FALSE, ...)
        
    }

    if(!is.null(cnr[["DDRC.nj"]])) {
        ape::write.tree(cnr[["DDRC.nj"]],
                        file = file.path(outdir, "DDRC.nj.newick"),
                        ...)
    }
    
} # end exportCNR
