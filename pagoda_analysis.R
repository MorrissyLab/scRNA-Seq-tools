
#install.packages("devtools")
#require(devtools)
#devtools::install_version('flexmix', '2.3-13')
#devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")
#biocLite("biomaRt")
#install.packages('XML')
#biocLite('GO.db')
#biocLite('Rtsne')
#install.packages('tidyverse')


library(scde)

#### THE INPUT DATA SHOULD BE A DENSE MATRIX
fetal <- read.csv(file = 'E:/visiting_research/data/fetal_brain_all/fetal_samples.csv', row.names = 1)
fetal <- as.matrix(fetal)
dim(fetal)

cd <- clean.counts(fetal)
dim(cd)


library(biomaRt)
library(GO.db)
# Initialize the connection to the Ensembl BioMart Service
# Available datasets can be listed with 
# listDatasets(useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))
# Use mmusculus_gene_ensembl for mouse
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
# Constructs a dataframe with two columns: hgnc_symbol and go_id
# If rownames are Ensembl IDs, use ensembl_gene_id as filter value
go <- getBM(attributes = c("hgnc_symbol","ensembl_gene_id", "go_id"), filters = "ensembl_gene_id", values = rownames(cd), mart = ensembl)
# Use the GO.db library to add a column with the GO-term to the dataframe
go$term <- Term(go$go_id)




############## CONVERT GENE ENSEMBLE ID INTO NAMES ###############
genes_map <- as.list(go$hgnc_symbol)
names(genes_map)<-go$ensembl_gene_id
#NOT UNIQUE MAP 

for (i in 1:length(rownames(cd))){
  if(rownames(cd)[i] %in% names(genes_map)){
    rownames(cd)[i] <- genes_map[rownames(cd)[i]]
  }
}

library(tidyverse) 


no_names_genes <- setdiff(rownames(cd),go$hgnc_symbol)
no_names_genes <- c(no_names_genes,'')
clean_rownames <- rownames(cd)[!rownames(cd) %in% no_names_genes]
cd <- cd[clean_rownames,]

####################################################################


knn <- knn.error.models(cd, k = ncol(cd)/4, n.cores = 1, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)
varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = TRUE)
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))

# list top overdispersed genes
sort(varinfo$arv, decreasing = TRUE)[1:10]



# Create a named list of character vectors out of the df
s = split(go$hgnc_symbol, paste(go$go_id,go$term))
#s = split(go$hgnc_symbol, paste(go$go_id,go$term))
### OR TODO Convert rownames of varinfo  into go_id


# Saves the list as a R environment
#go.env <- clean.gos(go.env) # remove GOs with too few or too many genes
go.env <- list2env(s)

# convert GO lists from ids to gene names
#gos.interest <- unique(c(ls(org.Hs.egGO2ALLEGS)[1:100],"GO:0022008","GO:0048699", "GO:0000280", "GO:0007067")) 
#go.env <- lapply(mget(gos.interest, org.Hs.egGO2ALLEGS), function(x) as.character(na.omit(rids[x]))) 
library(parallel)

#we can calculate weighted first principal component magnitudes for each GO gene set in the provided environment.
### 2 hours 
pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 1)

#We can now evaluate the statistical significance of the observed overdispersion for each GO gene set.
df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
write.csv(df, file = "E:/visiting_research/data/fetal_brain_all/pagoda/go_set_overdispersion.csv")

head(df)

save.image(file='pagoda.RData')
#load('pagoda.RData')

#Evaluate overdispersion of 'de novo' gene sets
clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = 1, plot = TRUE)

# Now the set of top aspects can be recalculated taking these de novo gene clusters into account:
df <- pagoda.top.aspects(pwpca, clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
write.csv(df, file = "E:/visiting_research/data/fetal_brain_all/pagoda/de_novo_gene_clusters.csv")
head(df)



save.image(file='pagoda.RData')
load('pagoda.RData')


## Visualize significant aspects of heterogeneity
# get full info on the top aspects
tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
# determine overall cell clustering
hc <- pagoda.cluster.cells(tam, varinfo)

tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)

tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)



col.cols <- rbind(groups = cutree(hc, 3))
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = col.cols)



# compile a browsable app, showing top 10 clusters with the top color bar
app <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols = col.cols, cell.clustering = hc, title = "pathway clustering")
# show app in the browser (port 1468)
show.app(app, "fetal_brain", browse = TRUE, port = 1468) 


save.image(file='pagoda.RData')


###############
pagoda.show.pathways(c("ENSG00000001461"), varinfo, go.env, cell.clustering = hc, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE)





#### Adding 2D embedding
library(Rtsne)
# recalculate clustering distance .. we'll need to specify return.details=T
cell.clustering <- pagoda.cluster.cells(tam,varinfo,include.aspects=TRUE,verbose=TRUE,return.details=T)
# fix the seed to ensure reproducible results
set.seed(0); 
tSNE.pagoda <- Rtsne(cell.clustering$distance,is_distance=T,initial_dims=100,perplexity=10)
par(mfrow=c(1,1), mar = c(2.5,2.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0);
plot(tSNE.pagoda$Y,col=adjustcolor(col.cols,alpha=0.5),cex=1,pch=19,xlab="",ylab="")




pp <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols = col.cols, cell.clustering = hc, title = "pathway", embedding = tSNE.pagoda$Y)
# show app in the browser (port 1468)
show.app(app, "fetal_brain", browse = TRUE, port = 1469) 





##########################
# get cell cycle signature and view the top genes
cc.pattern <- pagoda.show.pathways(c("GO:0001525"), varinfo, go.env, show.cell.dendrogram = TRUE, cell.clustering = hc, showRowLabels = TRUE)
# subtract the pattern
varinfo.cc <- pagoda.subtract.aspect(varinfo, cc.pattern)



