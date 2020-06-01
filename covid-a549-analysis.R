# Load packages
library(tximport)
library(DESeq2)
library(plotly)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

## NOT data.frame BUT tibble
sample_table = read_csv("SraRunTable.txt") %>% 
  select(`Sample Name`, source_name, treatment,
         Cell_Line, Cell_type, time_point) %>%
  slice(seq(1, 48, by=4))

sample_table = read_csv("SraRunTable.txt")
sample_table = select(sample_table, 
                      `Sample Name`, source_name, treatment,
                      Cell_Line, Cell_type, time_point)

sample_table = unique(sample_table)

## sample_table[sample_table$Cell_Line=='A549',]
sample_table_a549 = filter(sample_table, Cell_Line=='A549')

sample_files_a549 = paste0(
  pull(sample_table_a549, `Sample Name`), 
  '/quant.sf')

names(sample_files_a549) = pull(sample_table_a549,
                                `Sample Name`)

gene_map = read_csv("gene_map.csv", 
                    col_names = c('enstid', 'ensgid'))

count_data_a549 = tximport(files = sample_files_a549,
                           type = "salmon",
                           tx2gene = gene_map,
                           ignoreTxVersion = TRUE)


## 
sample_table_a549$conditions = factor(rep(c('mock', 'infected'), each=3), 
                                      levels=c('mock', 'infected'))
dds_a549 = DESeqDataSetFromTximport(txi=count_data_a549, 
                                    colData=sample_table_a549,
                                    design=~conditions)
#DESeq does:
# estimateSizeFactors
# estimateDispersions
# nbinomWaldTest
dds_a549 = DESeq(dds_a549)

vst_a549 = varianceStabilizingTransformation(dds_a549)
plotPCA(vst_a549, intgroup='conditions')

vst_mat = assay(vst_a549)
pca = prcomp(t(vst_mat))

df = as.data.frame(pca$x)
df$condition = sample_table_a549$conditions

pve = round(pca$sdev^2/sum(pca$sdev^2) * 100, 2)

ggplot(df, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour=condition)) +
  xlab(label = paste0("PC1 (", pve[1], "%)")) +
  ylab(label = paste0("PC2 (", pve[2], "%)"))


## Retrieve and filter results

a549_results = results(dds_a549, contrast = c('conditions',
                                              'infected',
                                              'mock'))
a549_df = data.frame(a549_results)
a549_df = rownames_to_column(a549_df, var='ensgene')
#a549_filter1 = a549_results[complete.cases(a549_results),]
a549_filter1 = filter(a549_df, complete.cases(a549_results))
a549_filter2 = filter(a549_filter1, padj < 0.05)
a549_filter3 = filter(a549_filter2, abs(log2FoldChange) > 1)

## plotting results

plotMA(a549_results)

a549_filter1$test = a549_filter1$padj < 0.05 & abs(a549_filter1$log2FoldChange) > 1

g = ggplot(a549_filter1, aes(x=log2FoldChange, 
                           y=-log10(padj), 
                           name=ensgene)) +
  geom_point(aes(colour=test), size=1, alpha=0.3) +
  scale_colour_manual(values=c('black', 'red')) +
  geom_vline(xintercept=1, colour='green', linetype=3) +
  geom_vline(xintercept=-1, colour='green', linetype=3) +
  geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
  theme_bw() +
  theme(legend.position = 'none')

ggplotly(g)


ensembl99 = useEnsembl(biomart="ensembl", version=99)

ensembl99 = useDataset("hsapiens_gene_ensembl", 
                       mart=ensembl99)

a549_anno = getBM(attributes=c('ensembl_gene_id',
                                'chromosome_name',
                                'start_position',
                                'end_position',
                                'strand',
                                'gene_biotype',
                                'external_gene_name',
                                'description'),
                   filters = c('ensembl_gene_id'),
                   values = a549_filter1$ensgene,
                   mart = ensembl99)

a549_anno_df = left_join(a549_filter1, a549_anno,
                         by=c('ensgene'='ensembl_gene_id'))
View(a549_anno_df)

g = ggplot(a549_anno_df, aes(x=log2FoldChange, 
                             y=-log10(padj), 
                             name=external_gene_name)) +
  geom_point(aes(colour=test), size=1, alpha=0.3) +
  scale_colour_manual(values=c('black', 'red')) +
  geom_vline(xintercept=1, colour='green', linetype=3) +
  geom_vline(xintercept=-1, colour='green', linetype=3) +
  geom_hline(yintercept=-log10(0.05), colour='blue', linetype=3) +
  theme_bw() +
  theme(legend.position = 'none')

ggplotly(g)

a549_anno2 = filter(a549_anno_df, padj < 0.05)
a549_anno3 = filter(a549_anno2, abs(log2FoldChange) > 1)

a549_degs = a549_anno3$ensgene
a549_hm = vst_mat[a549_degs,]
rownames(a549_hm) = a549_anno3$external_gene_name

heatmap(a549_hm)

pheatmap(a549_hm, fontsize_row=4, scale='row')

## GO enrichment

entrez_a549 = getBM(attributes=c('entrezgene_id'),
                 filters = c('ensembl_gene_id'),
                 values = a549_anno3$ensgene,
                 mart = ensembl99)
entrez_a549 = entrez_a549$entrezgene_id
entrez_a549 = as.character(entrez_a549)
a549_universe = getBM(attributes=c('entrezgene_id'),
                filters = c('ensembl_gene_id'),
                values = a549_anno_df$ensgene,
                mart = ensembl99)
a549_universe = a549_universe$entrezgene_id
a549_universe = as.character(a549_universe)

ego_a549 = enrichGO(gene = entrez_a549,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               universe = a549_universe,
               readable = TRUE)
barplot(ego_a549, showCategory=20)
dotplot(ego_a549, showCategory=20)

fold_changes = a549_anno3$log2FoldChange
names(fold_changes) = a549_anno3$external_gene_name

cnetplot(ego_a549, 
         showCategory=10, 
         foldChange = fold_changes)
goplot(ego_a549)








