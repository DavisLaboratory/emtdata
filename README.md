# emtdata

The emtdata package is an ExperimentHub package for three data sets with an Epithelial to Mesenchymal Transition (EMT). This package provides pre-processed RNA-seq data where the epithelial to mesenchymal transition was induced on cell lines. These data come from three publications [Cursons et al. (2015)](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA322427), [Cursons etl al. (2018)](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB25042) and [Foroutan et al. (2017)](https://doi.org/10.4225/49/5a2a11fa43fe3). In each of these publications, EMT was induces across multiple cell lines following treatment by TGFb among other stimulants. This data will be useful in determining the regulatory programs modified in order to achieve an EMT. Data were processed by the Davis laboratory in the Bioinformatics division at WEHI.

This package can be installed using the code below:

---
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") BiocManager::install("emtdata")
---

## Documentation

The package comes with a vignette showing how the different functions in the package can be used to download and access the data. Pre-built vignettes can be accessed via [Bioconductor](https://www.bioconductor.org/packages/release/data/experiment/vignettes/emtdata/inst/doc/emtdataR.html).
