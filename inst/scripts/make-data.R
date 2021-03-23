library(SummarizedExperiment)
library(edgeR)
library(BiocFileCache)
library(Homo.sapiens)
library(stringr)
library(plyr)

#----Cursons et al. (2018) Cell Systems, HMLE system EMT----
str_dataset = 'Cursons2018'

#specify paths
fcpath = 'inst/extdata/Cursons2018_fc_counts.txt.gz'
annotpath = 'inst/extdata/Cursons2018_sample_annotations.csv'

#read in raw data
emat = read.delim(gzfile(fcpath), skip = 1, row.names = 1, check.names = FALSE)
sampleannot = read.csv(annotpath)
rownames(sampleannot) = sampleannot$Run

#split gene annotations from the emat
geneannot = emat[, 1:7]
emat = emat[, -(1:7)]
colnames(emat) = basename(dirname(colnames(emat)))
emat = emat[, rownames(sampleannot)]

#create DGEList object for pre-processing
dge = DGEList(counts = emat, genes = geneannot, samples = sampleannot)

#filter out lowly expressed genes
design = model.matrix(~ 0 + Subline + Treatment, data = dge$samples)
keep = filterByExpr(emat, design = design)

op = par(no.readonly = TRUE)
par(mfrow = c(1, 2))
hist(cpm(dge, log = TRUE), main = paste('Before Filter:', str_dataset), xlab = 'logCPM')
hist(cpm(dge[keep, ], log = TRUE), main = paste('After Filter:', str_dataset), xlab = 'logCPM')
par(op)

dge = dge[keep, ]

#compute normalisation factors and RPKM
dge = calcNormFactors(dge)
emat_rpkm = rpkm(dge, gene.length = 'Length', log = TRUE)

#create a SummarizedExperiment object
cursons2018_se = SummarizedExperiment(
  assays = list('counts' = dge$counts, 'logRPKM' = emat_rpkm),
  rowData = dge$genes,
  colData = dge$samples
)
colnames(cursons2018_se) = cursons2018_se$Sample.Name
save(cursons2018_se, file = 'cursons2018_se.rda')

#----Rik Thompson 2015 EMT----
str_dataset = 'Cursons2015'

#specify paths
fcpath = 'inst/extdata/Cursons2015_fc_counts.txt.gz'
annotpath = 'inst/extdata/Cursons2015_sample_annotations.csv'

#read in raw data
emat = read.table(gzfile(fcpath),row.names=1, check.names = FALSE)
sampleannot = read.csv(annotpath)
rownames(sampleannot) = sampleannot$Run

#split gene annotations from the emat
geneannot = emat[, 1:7]
emat = emat[, -(1:7)]
colnames(emat) = basename(dirname(colnames(emat)))
emat = emat[, rownames(sampleannot)]

#sum up technical replicates
emat = sumTechReps(emat, ID = sampleannot$Sample.Name)

#combine annotations for technical replicates
sampleannot = ddply(sampleannot, 'Sample.Name', function(df) {
  df$Run = paste(df$Run, collapse = ',')
  df$Experiment = paste(df$Experiment, collapse = ',')
  return(unique(df))
})
rownames(sampleannot) = sampleannot$Sample.Name
emat = emat[, rownames(sampleannot)]

#create DGEList object for pre-processing
dge_cursons = DGEList(counts = emat, genes = geneannot, samples = sampleannot)

#filter out lowly expressed genes
#no treatment differences in this experiment as the samples are all untreated controls
design = model.matrix(~ 0 + Cell.Line + Treatment, data = dge_cursons$samples)
keep = filterByExpr(emat, design = design)

op = par(no.readonly = TRUE)
par(mfrow = c(1, 2))
hist(cpm(dge_cursons, log = TRUE), main = paste('Before Filter:', str_dataset), xlab = 'logCPM')
hist(cpm(dge_cursons[keep, ], log = TRUE), main = paste('After Filter:', str_dataset), xlab = 'logCPM')
par(op)

dge_cursons = dge_cursons[keep, ]

#compute normalisation factors and RPKM
dge_cursons = calcNormFactors(dge_cursons)
emat_rpkm = rpkm(dge_cursons, gene.length = 'Length', log = TRUE)

#create a SummarizedExperiment object
cursons2015_se = SummarizedExperiment(
  assays = list('counts' = dge_cursons$counts, 'logRPKM' = emat_rpkm),
  rowData = dge_cursons$genes,
  colData = dge_cursons$samples
)
save(cursons2015_se, file = 'cursons2015_se.rda')

#----Foroutan et al. (2017) MCR, EMT compendium dataset----
figshare_link = 'https://ndownloader.figshare.com/files/9938620'
annotpath = 'inst/extdata/Foroutan_sample_annotations.csv'

#download data
bfc = BiocFileCache(tempfile(), ask = FALSE)
fpath = bfcrpath(bfc, figshare_link)
foroutan_data = read.delim(fpath, row.names = 1)

#read in sample annotations
sampleannot = read.csv(annotpath)
rownames(sampleannot) = sampleannot$SampleID

#create gene annotations
geneannot = data.frame(
  'EntrezID' = rownames(foroutan_data),
  'EnsemblID' = mapIds(Homo.sapiens, rownames(foroutan_data), 'ENSEMBL', 'ENTREZID'),
  'Chr' = mapIds(Homo.sapiens, rownames(foroutan_data), 'CDSCHROM', 'ENTREZID'),
  'Start' = mapIds(Homo.sapiens, rownames(foroutan_data), 'CDSSTART', 'ENTREZID'),
  'End' = mapIds(Homo.sapiens, rownames(foroutan_data), 'CDSEND', 'ENTREZID'),
  'Strand' = mapIds(Homo.sapiens, rownames(foroutan_data), 'CDSSTRAND', 'ENTREZID'),
  'gene_name' = mapIds(Homo.sapiens, rownames(foroutan_data), 'SYMBOL', 'ENTREZID')
)

#create a SummarizedExperiment object
foroutan2017_se = SummarizedExperiment(
  assays = list('logExpr' = as.matrix(foroutan_data)),
  rowData = geneannot,
  colData = sampleannot[colnames(foroutan_data), ]
)
save(foroutan2017_se, file = 'foroutan2017_se.rda')

