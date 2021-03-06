% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\name{cursons2015_se}
\alias{cursons2015_se}
\title{EGF or hypoxia treatment of Breast Cancer Cell lines from Cursons et al. (2015)}
\format{
A SummarizedExperiment object, containing gene expression data of
different sub-lines of human breast cancer cell lines. The
\code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment::colData()}} function can be used to access the sample
annotations.
}
\description{
Dataset used in Cursons et al. 2015. where parallel ‘deep sequencing’ of RNA
(RNA-Seq) were conducted to examine the changes in expression profiles between
breast cancer cell lines PMC42-ET, PMC42-LA and MDA-MB-468 cells.
}
\details{
For PMC42 cell lines, cells were treated for 3 or 7 days in the presence or
absence of 10 ng/ml EGF.
For MDA-MB-468 cells, cells were treated for 7 days in the presence or absence
of either 10 ng/ml EGF or kept under Hypoxia.

There are 3 biological replicates per condition which were summed.
Samples were sequenced on the Illumina HiSeq 2000, 100bp paired end.

Data from this publication were downloaded from EMBL-EBI ENA and processed into
a SummarizedExperiment object. Sample annotations were modified from the original
publication and SRA portal.
}
\examples{
library(ExperimentHub)
eh <- ExperimentHub()
emt_datasets <- query(eh, "emtdata")
}
\references{
Cursons, J., Leuchowius, KJ., Waltham, M., Tomaskovic-Crook, E.,
Foroutan, M., Bracken, CP., Redfern,A., Crampin, EJ., Street, I., Davis, MJ.
& Thompson,EW. (2015) Stimulus-dependent differences in signalling regulate
epithelial-mesenchymal plasticity and change the effects of drugs in breast
cancer cell lines. Cell Commun Signal 13, 26 (2015).
}
