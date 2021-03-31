#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('_____')


#Load the edgeR library
library(topGO)
library(edgeR)
library(GO.db)
library(reshape2)
library(ggplot2)

#Import gene count data
countsTable <- read.csv(file= "geneCounts_cleaned_PA42_v4.1.csv", 
                        row.names="gene")[ ,1:6]
#Add grouping factor
group <- factor(c(rep("ctrl",3),rep("treat",3)))
#Create DGE list object
list <- DGEList(counts=countsTable,group=group)

#There is no purpose in analyzing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
list <- list[keep, , keep.lib.sizes=FALSE]
#Calculate normalized factors
list <- calcNormFactors(list)

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)

#Perform an exact test for treat vs ctrl
tested <- exactTest(list, pair=c("ctrl", "treat"))

#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table

#GO enrichment
#Read in custom GO annotations
GOmaps <- readMappings(file='gene2GO_PA42_v4.1_transcripts.map',  sep='\t',  IDsep=',')

#Create named list of all genes (gene universe) and p-values. This time, the gene universe is set
#to be only the genes with p-values <= 0.05 that are also contained in the list of gene2GO annotated genes.
DGE_results_table <- tested$table
DGE_results_table_significant <- DGE_results_table[DGE_results_table$PValue <= 0.05, ]
list_genes <- as.numeric(DGE_results_table_significant$PValue)
list_genes <- setNames(list_genes, rownames(DGE_results_table_significant))
list_genes_filtered <- list_genes[names(list_genes) %in% names(GOmaps)]

#Read in list of interesting UV related genes for Daphnia that we want to do the analysis on. 
UV_gene_list <- read.csv('DDRGOTF_Dmel_PA42_v4.1_combined_geneIDs_uniq.csv')
UV_gene_list <- as.list(UV_gene_list)$gene

#Create function to return list of interesting DE genes (0 == gene in geneUniverse was not 
#contained in UV gene list, 1 == gene in geneUniverse was also contained in list of UV genes of interest)
get_interesting_DE_genes <- function(geneUniverse){
  interesting_DE_genes <- rep(0, length(geneUniverse))
  interesting_DE_genes[which(names(geneUniverse) %in% UV_gene_list)] = 1
  interesting_DE_genes <- setNames(interesting_DE_genes, names(geneUniverse))
  return(interesting_DE_genes)
}

#Create topGOdata objects for enrichment analysis (1 for each ontology)
BP_GO_data <- new('topGOdata', ontology = 'BP', allGenes = list_genes_filtered, 
                  geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                  gene2GO = GOmaps)
MF_GO_data <- new('topGOdata', ontology = 'MF', allGenes = list_genes_filtered, 
                  geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                  gene2GO = GOmaps)
CC_GO_data <- new('topGOdata', ontology = 'CC', allGenes = list_genes_filtered, 
                  geneSel = get_interesting_DE_genes, nodeSize = 10, annot = annFUN.gene2GO, 
                  gene2GO = GOmaps)

#Summary functions
#numGenes(BP_GO_data)
#sigGenes(BP_GO_data)

#PerformGO enrichment on topGOdata object
BP_GO_results <- runTest(BP_GO_data, statistic = 'Fisher')
MF_GO_results <- runTest(MF_GO_data, statistic = 'Fisher')
CC_GO_results <- runTest(CC_GO_data, statistic = 'Fisher')




