meta = data.frame(
  Title = c('foroutan2017_se', 'cursons2018_se', 'cursons2015_se'),
  Description = c(
    'Gene expression data from Foroutan et al., MCR 2017. This gene expression data comes from multiple different studies (microarary and RNA-seq), with cell lines treated using TGFb to induce a mesenchymal shift. Data were combined using SVA and ComBat to remove batch effects.',
    'Gene expression data from Cursons et al., Cell Syst 2018. This gene expression data comes from the human mammary epithelial (HMLE) cell line. A mesenchymal HMLE (mesHMLE) phenotype was induced following treatment with TGFb. The mesHMLE subline was then treated with mir200c to reinduce an epithelial phenotype.',
    'Gene expression data from Cursons et al., Cell Commun Signal 2015. This gene expression data comes from the PMC42-ET, PMC42-LA and MDA-MB-468 cell lines. Mesenchymal phenotype was induced in PMC42 cell lines with EGF treatment and in MDA-MB-468 with either EGF treatment or kept under Hypoxia.'
  ),
  BiocVersion = c(3.13),
  Genome = NA,
  SourceType = c('TXT'),
  SourceUrl = c(
    'https://doi.org/10.4225/49/5a2a11fa43fe3',
    'https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB25042',
    'https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA322427'
  ),
  SourceVersion = c('1.0'),
  Species = c('Homo sapiens'),
  TaxonomyId = c(9606),
  Coordinate_1_based = TRUE,
  DataProvider = c('Walter and Eliza Hall Institute of Medical Research', 'Queensland University of Technology', 'Center for Cancer Biology, Adelaide, Australia'),
  Maintainer = 'Malvika D. Kharbanda <kharbanda.m@wehi.edu.au>',
  RDataClass = 'GSEABase::SummarizedExperiment',
  DispatchClass = 'Rda',
  RDataPath = c(
    'emtdata/foroutan2017_se.rda',
    'emtdata/cursons2018_se.rda',
    'emtdata/cursons2015_se.rda'
  )
)

write.csv(meta, file = 'inst/extdata/metadata.csv', row.names = FALSE)

ExperimentHubData::makeExperimentHubMetadata('../emtdata')

