#' TGFb stimulated cell lines from Foroutan et al. (2017)
#'
#' Foroutan et al. collected data from multiple different studies where various
#' cell lines (from various tissues of origin) were stimulated with transforming
#' growth factor beta (TGFb) to induce an epithelial to mesenchymal transition
#' (EMT). Since the data were from different studies (and platforms), the
#' surrogate variable analysis (SVA) and ComBat approaches were applied to
#' correct for batch effects.
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
#' @name foroutan_se
#'
NULL

#' Breast Cancer Cell lines PMC42-ET and PMC42-LA from Cursons et al. (2015)
#'
#' Unpublished dataset from the laboratory of Rik Thompson at the Queensland University of Technology (QUT).
#' Dataset used in Cursons et al. 2015. where parallel ‘deep sequencing’ of RNA (RNA-Seq) were conducted
#' to examine the changes in expression profiles between PMC42-ET and PMC42-LA cells. Cells were treated
#' for 3 days in the presence or absence of 10 ng/ml EGF. ONLY the controls are processed.
#' There are 3 biological replicates per condition which were summed.
#' Sequenced on the Illumina HiSeq 2000, 100bp paired end.
#'
#' Data from this publication were downloaded from SRA and/or EMBL-EBI ENA and processed into a
#' SummarizedExperiment object. Sample annotations were modified from the original publication.
#'
#' @format A SummarizedExperiment object, containing gene expression data of
#'   different sub-lines of human breast cancer cell lines. The [SummarizedExperiment::colData()]
#'   function can be used to access the sample annotations.
#' @references Cursons, J., Leuchowius, KJ., Waltham, M., Tomaskovic-Crook, E., Foroutan, M.,
#' Bracken, CP., Redfern,A., Crampin, EJ., Street, I., Davis, MJ. & Thompson,EW. (2015)
#' Stimulus-dependent differences in signalling regulate epithelial-mesenchymal plasticity and
#' change the effects of drugs in breast cancer cell lines. Cell Commun Signal 13, 26 (2015).
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' emt_datasets <- query(eh, "emtdata")
#' @name thompson2015_se
#'
NULL
