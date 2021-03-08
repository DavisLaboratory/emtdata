library(SummarizedExperiment)
library(edgeR)

#----Cursons et al. Cell Systems, HMLE system EMT----
#specify paths
fcpath = "inst/extdata/Cursons_fc_counts.txt.gz"
annotpath = "inst/extdata/Cursons_sample_annotations.csv"

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
hist(cpm(dge, log = TRUE))
hist(cpm(dge[keep, ], log = TRUE))
par(op)

dge = dge[keep, ]

#compute normalisation factors and RPKM
dge = calcNormFactors(dge)
emat_rpkm = rpkm(dge, gene.length = 'Length', log = TRUE)

#create a SummarizedExperiment object
cursons_se = SummarizedExperiment(
  assays = list('counts' = dge$counts, 'logRPKM' = emat_rpkm),
  rowData = dge$genes,
  colData = dge$samples
)
colnames(cursons_se) = cursons_se$Sample.Name
save(cursons_se, file = 'cursons_se.rda')

#----Rik Thompson EMT----

#----Foroutan et al. , EMT compendium dataset----

