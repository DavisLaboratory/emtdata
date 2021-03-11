#' TGFb stimulated cell lines from Foroutan et al. (2017)
#'
#' Foroutan et al. collected data from multiple different studies where various
#' cell lines (from various tissues of origin) were stimulated with transforming
#' growth factor beta (TGFb) to induce an epithelial to mesenchymal transition
#' (EMT). Since the data were from different studies (and platforms), the
#' removal of unwanted variation (RUV) approach was applied to discard batch
#' effects.
#'
#' Data from this publication were downloaded from figshare and processed into a
#' SummarizedExperiment object. Sample annotations were retrieved from the
#' original publication.
#'
#' @format A SummarizedExperiment object, containing gene expression data of
#'   human cell lines treated with TGFb. The [SummarizedExperiment::colData()]
#'   function can be used to access the sample annotations.
#' @references Foroutan, M., Cursons, J., Hediyeh-Zadeh, S., Thompson, E. W., &
#'   Davis, M. J. (2017). A transcriptional program for detecting TGFβ-induced
#'   EMT in Cancer. Molecular Cancer Research, 15(5), 619-631.
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' emt_datasets <- query(eh, "emtdata")
#'
"foroutan_se"
