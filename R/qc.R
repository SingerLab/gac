#' Quality control annotation
#'
#' @description
#' A dataset containing the the quality control of single-cell data. The Y matrix
#' and QC matrix are quite similar.  I like to think that Y are specific
#' experimental phenotypes that one would like to gain biological insights into.
#' QC is more adequate for technical variation, weather from manipulation and
#' batch processing of samples, and/or biological such as expected ploidy from
#' flow cytometry indexed FCS file.
#'
#' I think that good QC data is a reflection of good laboratory notes, and
#' practices.  They can provide additional support to make decisions such as
#' including or excluding certain samples. Hence why it's a required parameter in
#' the CNR.
#' 
#'
#' @format A data frame n.cell rows x qc.metrics columns
#'
#' \itemize{
#'   \item cellID, single-cell ID as the rownames ! this is imperative or it
#' will be re-written by the function
#' 
#'   \item ReadsKept, Number of reads kept after alignment
#' 
#'   \item MedianBinCount, Median number of reads per bin
#' 
#'   \item dna.ng.ul, DNA concentration of the cell post amplification
#' 
#'   \item sort.gate, 2N, 3N, or 4N gate sorted in FACS. Sets expected value of
#' ploidy. Not required but highly recomended
#' 
#'   \item qc.status, binary or multinomial PASS/FAIL/WARNING call for each cell;
#' it's based on your criteria, though some recomendations are provided in the
#' use vignette
#' }
#' @docType data
#' @keywords datasets
#' @usage data(qc)
#' @source \url{https://github.com/SingerLab/GAC}
"qc"
