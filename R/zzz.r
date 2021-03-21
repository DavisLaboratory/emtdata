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

#' EGF or hypoxia treatment of Breast Cancer Cell lines from Cursons et al. (2015)
#'
#' Dataset used in Cursons et al. 2015. where parallel ‘deep sequencing’ of RNA
#' (RNA-Seq) were conducted to examine the changes in expression profiles between
#' breast cancer cell lines PMC42-ET, PMC42-LA and MDA-MB-468 cells.
#'
#' For PMC42 cell lines, cells were treated for 3 or 7 days in the presence or
#' absence of 10 ng/ml EGF.
#' For MDA-MB-468 cells, cells were treated for 7 days in the presence or absence
#' of either 10 ng/ml EGF or kept under Hypoxia.
#'
#' There are 3 biological replicates per condition which were summed.
#' Samples were sequenced on the Illumina HiSeq 2000, 100bp paired end.
#'
#' Data from this publication were downloaded from EMBL-EBI ENA and processed into
#' a SummarizedExperiment object. Sample annotations were modified from the original
#' publication and SRA portal.
#'
#' @format A SummarizedExperiment object, containing gene expression data of
#'   different sub-lines of human breast cancer cell lines. The
#'   [SummarizedExperiment::colData()] function can be used to access the sample
#'   annotations.
#' @references Cursons, J., Leuchowius, KJ., Waltham, M., Tomaskovic-Crook, E.,
#' Foroutan, M., Bracken, CP., Redfern,A., Crampin, EJ., Street, I., Davis, MJ.
#' & Thompson,EW. (2015) Stimulus-dependent differences in signalling regulate
#' epithelial-mesenchymal plasticity and change the effects of drugs in breast
#' cancer cell lines. Cell Commun Signal 13, 26 (2015).
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' emt_datasets <- query(eh, "emtdata")
#' @name cursons2015_se
#'
NULL

#' Combinatorial miRNAs in breast cancer EMT from Cursons et al. (2018)
#'
#' Cursons et al. used the Human Mammary Epithelial Cells (HMLE) cell line dataset.
#' A mesenchymal HMLE (mesHMLE) phenotype was induced following treatment with
#' transforming growth factors (TGFb). The mesHMLE subline was then treated with
#' mir200c to reinduce an epithelial phenotype. All mRNA RNA-seq was collected
#' using the Illumina HiSeq 2500 with a paired end read length of 100bp.
#'
#' Data from this publication were downloaded from the European Nucleotide
#' Archive (ENA) and processed into a SummarizedExperiment object. Sample
#' annotations were modified from the original publication.
#'
#' @format A SummarizedExperiment object, containing gene expression data of
#'   the Human Mammary Epithelial Cells (HMLE) cell line. The
#'   [SummarizedExperiment::colData()] function can be used to access the sample
#'   annotations.
#' @references Cursons, J., Pillman, K. A., Scheer, K. G., Gregory, P. A.,
#' Foroutan, M., Hediyeh-Zadeh, S., ... & Davis, M. J. (2018). Combinatorial
#' targeting by microRNAs co-ordinates post-transcriptional control of EMT. Cell
#' systems, 7(1), 77-91.
#' @examples
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' emt_datasets <- query(eh, "emtdata")
#' @name Cursons2018_se
#'
NULL
