meta = data.frame(
  Title = c('foroutan_se'),
  Description = c(
    'Gene expression data from Foroutan et al., MCR 2017. This gene expression data comes from multiple different studies (microarary and RNA-seq), with cell lines treated using TGFb to induce a mesenchymal shift. Data were combined using RUV to remove batch effects.'
  ),
  BiocVersion = c(3.13),
  Genome = NA,
  SourceType = c('TXT'),
  SourceUrl = c(
    'https://doi.org/10.4225/49/5a2a11fa43fe3'
  ),
  SourceVersion = c('1.0'),
  Species = c('Homo sapiens'),
  TaxonomyId = c(9606),
  Coordinate_1_based = TRUE,
  DataProvider = c('Walter and Eliza Hall Institute of Medical Research'),
  Maintainer = 'Malvika D. Kharbanda <kharbanda.m@wehi.edu.au>',
  RDataClass = 'GSEABase::SummarizedExperiment',
  DispatchClass = 'Rda',
  RDataPath = c(
    'emtdata/foroutan_se.rda'
  )
)

write.csv(meta, file = 'inst/extdata/metadata.csv', row.names = FALSE)

ExperimentHubData::makeExperimentHubMetadata('../emtdata')
