################# INSTALL PACKAGES ###############

source("http://bioconductor.org/biocLite.R")
biocLite("monocle")


source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")


################# Load Libraries ##################
library(monocle)
library(cellrangerRkit)


################# Load DATA ########################

#### EDIT THIS PATH TO YOUR DATA
acellranger_pipestance_path <- "E:/visiting_research/myeloma/scRNAseq/P1217_MM2"

mm2_data <- load_cellranger_matrix(acellranger_pipestance_path)
fd <- fData(mm2_data)
colnames(fd)[2] <- "gene_short_name"

HSMM <- newCellDataSet(exprs(mm2_data),
                      phenoData = new("AnnotatedDataFrame", data = pData(mm2_data)),
                      featureData = new("AnnotatedDataFrame", data = fd),
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

#################### Clean Data #########################

HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)


######################## Detect Genes ####################
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM)))

HSMM <- setOrderingFilter(HSMM, expressed_genes)
plot_ordering_genes(HSMM)



#################### Build Tree #################
HSMM <- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM,color_by = "time")

#######################################################