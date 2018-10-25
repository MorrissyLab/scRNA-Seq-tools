

####################### HELPER FUNCTIONS #########################





####################### INSTALL HONEYBADGER ########################
install.packages("rjags")
library(rjags)
source("https://bioconductor.org/biocLite.R")
install.packages("digest")
library(digest)
biocLite(c("GenomicRanges", "AnnotationDbi", "GenomicFeatures", "TxDb.Hsapiens.UCSC.hg19.knownGene"))
install.packages("rlang")
library(rlang)
install.packages("devtools")
install.packages("Rcpp")
install.packages("HiddenMarkov")
require(devtools)
devtools::install_github('JEFworks/HoneyBADGER')
library(HoneyBADGER)

##################### INSTALLL CELLRANGER ##############

source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
library(cellrangerRkit)


################## LOAD DATA ###########
#the Gene / cell matrix (filtered) can be normalized to CPMs and log transformmed to serve as the gene expression matrix. 



# EDIT THIS PATH TO YOUR DATA FOLDER
p1217_path <- "/myeloma/scRNAseq/P1217_MM2"


mm_p1217 <- load_cellranger_matrix(p1217_path)
mt_p1217 <- exprs(mm_p1217)

################ NORMALIZE MATRIX ##############
use_genes <- get_nonzero_genes(mm_p1217)
mm_bcnorm <- normalize_barcode_sums_to_median(mm_p1217[use_genes,])
mm_log <- log_gene_bc_matrix(mm_bcnorm,base=10)
print(dim(mm_log))


## Convert sparse gene expression to dense 
gexp <- as.matrix(exprs(mm_log))
#rmt_gexp <- t(gexp)

################# GETTING STARTED and replace genes name ###############
#data(gexp) ## tumor cells

data(ref) 

require(biomaRt) ## for gene coordinates
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org")


go <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id"), filters = "ensembl_gene_id", values = rownames(gexp), mart =mart.obj)
rownames(gexp) <- mapvalues(rownames(gexp), go$ensembl_gene_id, to=go$hgnc_symbol)

############### CREATE HONEY BADGER #############
hb <- new('HoneyBADGER', name='P1217_MM2')
hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE, verbose=TRUE)


hb$plotGexpProfile() ## initial visualization




hb$setMvFit(verbose=TRUE) ## model variance
hb$setGexpDev(verbose=TRUE) ## model necessary expression deviation to identify CNVs
hb$calcGexpCnvBoundaries(init=TRUE, verbose=FALSE) ## HMM


########## IDENTIFY CNV ##################
bgf <- hb$bound.genes.final
genes <- hb$genes
regions.genes <- range(genes[unlist(bgf)])
print(regions.genes)
########### CLUSTER AND TEST