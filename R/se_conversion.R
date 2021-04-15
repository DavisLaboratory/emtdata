#' Convert a SummarizedExperiment object to a DGEList object
#'
#' This function converts a SummarizedExperiment object to a DGEList object to
#' enhance differential expression analysis using the edgeR package.
#'
#' @param se a SummarizedExperiment object.
#' @param assay_name a character or a numeric, specifying the assay to retrieve.
#'
#' @return a DGEList object
#' @export
#'
#' @examples
#' library(ExperimentHub)
#'
#' eh = ExperimentHub()
#' query(eh, 'emtdata')
#'
#' cursons2018_se = eh[['EH5440']]
#' cursons2018_dge = asDGEList(cursons2018_se)
#'
asDGEList <- function(se, assay_name = 'counts') {
  dge = edgeR::DGEList(
    counts = SummarizedExperiment::assay(se, assay_name),
    genes = SummarizedExperiment::rowData(se),
    samples = SummarizedExperiment::colData(se)
  )
  return(dge)
}
