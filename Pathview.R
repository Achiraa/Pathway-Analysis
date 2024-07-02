setwd("C:/geo_data/Deseq")

library(fgsea)
library(dplyr)
library(BiasedUrn)
library(org.Mm.eg.db)
library("AnnotationDbi")
library(DESeq2)
library(GenomicRanges)
library(GenomeInfoDb)
library(SummarizedExperiment)
library(MatrixGenerics)
library(matrixStats)


columns(org.Mm.eg.db)

gg=readRDS("C:/geo_data/Deseq/resu.01.rds")
gg

gg$symbol = mapIds(org.Mm.eg.db,
                    keys=row.names(gg), 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
gg$entrez = mapIds(org.Mm.eg.db,
                    keys=row.names(gg), 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
gg$name =   mapIds(org.Mm.eg.db,
                    keys=row.names(gg), 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")

head(gg, 20)

library(pathview)
library(gage)
library(gageData)

data(kegg.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm, 3)

foldchanges = gg$log2FoldChange
names(foldchanges) = gg$entrez
head(foldchanges)

# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=TRUE)
keggres

# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways

# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", new.signature=FALSE)
plot_pathway

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu"))













