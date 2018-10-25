
if (!require(viridis)) install.packages("viridis")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(readr)) install.packages("readr")
if (!require(phateR)) install.packages("phateR")

devtools::install_github("KrishnaswamyLab/magic/Rmagic")

library(ggplot2)
library(readr)
library(viridis)
library(Rmagic)
library(phateR)
library(biomaRt)
library(dplyr)


#### THE INPUT DATA SHOULD BE A DENSE TRANSPOSED MATRIX
bmmsc <- read_csv("E:/visiting_research/data/fetal_brain_all/fetal_samples_t.csv")
bmmsc <- bmmsc[,2:ncol(bmmsc)]
bmmsc[1:5,1:10]
# keep genes expressed in at least 10 cells
keep_cols <- colSums(bmmsc > 0) > 10
bmmsc <- bmmsc[,keep_cols]
# look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=rowSums(bmmsc)), bins=10) +
  geom_vline(xintercept = 33000, color='red')
bmmsc <- library.size.normalize(bmmsc)
bmmsc <- sqrt(bmmsc)


############ CONVERT GENES NAMES
human_mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
h_list<- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), filters="ensembl_gene_id",values=colnames(bmmsc),mart=human_mart, uniqueRows=TRUE)
genes_map <- as.list(h_list[,"hgnc_symbol"])
names(genes_map)<-h_list[,"ensembl_gene_id"]
genes_map['ENSG00000181449']

for (i in 1:length(colnames(bmmsc))){
  ### TODO KEEP NOT FOUND 
  if (colnames(bmmsc)[i] %in% names(genes_map)){
    colnames(bmmsc)[i] <- genes_map[colnames(bmmsc)[i]]
    
  }
  
}



############### MAGIC ON SEVERAL GENES
clust_genes1 <- c("SOX2","NES","OLIG2","OLIG1")
clust_genes2 <- c("CCL4","CCL3","AIF1")
clust_genes3 <- c("PLP1","MBP")
clust_genes4 <- c("RORA","LHX1","PCP4")

bmmsc_MAGIC <- magic(bmmsc,genes=clust_genes2)
### BEFORE MAGIC

p1<- ggplot(bmmsc) +
  geom_point(aes(x=CCL4, y=CCL3, color=AIF1)) +
  scale_color_viridis(option="B") +
  labs(color="AIF1") +
  ggtitle("Before MAGIC")

p2<- ggplot(bmmsc_MAGIC) +
  geom_point(aes(x=CCL4, y=CCL3, color=AIF1)) +
  scale_color_viridis(option="B") +
  labs(color="AIF1") +
  ggtitle("After MAGIC")

multiplot(p1, p2, cols=1)


bmmsc_MAGIC <- magic(bmmsc,genes=clust_genes1)
p1<- ggplot(bmmsc) +
  geom_point(aes(x=SOX2, y=OLIG2, color=NES)) +
  scale_color_viridis(option="B") +
  labs(color="NES") +
  ggtitle("Before MAGIC")
p2<- ggplot(bmmsc_MAGIC) +
  geom_point(aes(x=SOX2, y=OLIG2, color=NES)) +
  scale_color_viridis(option="B") +
  labs(color="NES") +
  ggtitle("After MAGIC")

multiplot(p1, p2, cols=1)

#########################
bmmsc_MAGIC_PCA <- magic(bmmsc, genes="pca_only", t=4)


ggplot(bmmsc_MAGIC_PCA) +
  geom_point(aes(x=PC1, y=PC2, color=bmmsc_MAGIC$result$NES)) +
  scale_color_viridis(option="B") +
  labs(color="NES")


ggplot(bmmsc_MAGIC_PCA) +
  geom_point(aes(x=PC1, y=PC2, color=bmmsc_MAGIC$result$NES)) +
  scale_color_viridis(option="B") +
  labs(color="SOX2")


ggplot(bmmsc_MAGIC_PCA) +
  geom_point(aes(x=PC1, y=PC2, color=bmmsc_MAGIC$result$NES)) +
  scale_color_viridis(option="B") +
  labs(color="OLIG2")

############## PHATE
bmmsc_PHATE <- phate(bmmsc, k=4, a=100, t=20)
ggplot(bmmsc_PHATE) +
  geom_point(aes(x=PHATE1, y=PHATE2, color=bmmsc_MAGIC$result$AIF1)) +
  scale_color_viridis(option="B") +
  labs(color="AIF1")


########### ALL GENES
bmmsc_MAGIC <- magic(bmmsc, genes="all_genes", t=4)
df_magic <- as.data.frame(bmmsc_MAGIC)
write.csv(df_magic, file = "E:/visiting_research/data/fetal_brain_all/magic_df.csv",row.names=F)

