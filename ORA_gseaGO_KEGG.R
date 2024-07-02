library(enrichplot)
library(ggplot2)
library("org.Mm.eg.db")
library(clusterProfiler)

setwd("D:/geo_data/stromalDeseq")

res=read.csv("forgsea.csv")
head(res)

# we want the log2 fold change 
original_gene_list1 = res$log2FoldChange
original_gene_list1

# name the vector
names(original_gene_list1) <- res$geneids

# sort the list in decreasing order (required for clusterProfiler)
gene_list1 = sort(original_gene_list1, decreasing = TRUE)
View(gene_list1)

gse= gseGO(geneList=gene_list1, 
      ont ="ALL", 
      keyType = "ENTREZID",
      minGSSize = 120, 
      maxGSSize = 800, 
      pvalueCutoff = 0.05, 
      verbose = TRUE, 
      OrgDb = org.Mm.eg.db, 
      pAdjustMethod = "none")
head(gse)

write.csv(gse,file = "GSEA excel stroma.csv")

dotplot(gse, showCategory=20) + ggtitle("Dotplot for GSEA:Gene Ontology")

res11=read.csv("forgseakegg.csv")
extr=res11$log2FoldChange
View(extr)

# name the vector
names(extr) <- res11$ENTREZID_id
head(extr)
View(extr)

kegg_gene_lists = sort(extr, decreasing = TRUE)

kegg_organism = "mmu"
gse1=gseKEGG(geneList = kegg_gene_lists,
             organism = kegg_organism,
             nPerm = 10000,
             minGSSize = 120,
             maxGSSize = 800,
             pvalueCutoff = 0.1,
             pAdjustMethod = "none",
             keyType = "ncbi-geneid",
             verbose = FALSE)
gse1


dotplot(gse1, showCategory=20) + ggtitle("Dotplot for GSEA:KEGG Pathway")
