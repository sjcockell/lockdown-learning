# Load packages
library(tximport)
library(DESeq2)
library(plotly)
library(tidyverse)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)

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

'GSM4432378/quant.sf'

paste0(c('Simon', 'S', 'SJ'), 'Cockell')

sample_files = paste0(
  pull(sample_table, `Sample Name`), 
  '/quant.sf')

names(sample_files) = pull(sample_table, `Sample Name`)

gene_map = read_csv("gene_map.csv", 
                    col_names = c('enstid', 'ensgid'))

count_data = tximport(files = sample_files,
         type = "salmon",
         tx2gene = gene_map,
         ignoreTxVersion = TRUE)

count_data[['counts']]
count_data$counts[1:6,]

sample_table = as.data.frame(sample_table)
colnames(sample_table)[1] = "Sample"

conditions = c('mock_nhbe', 'infected_nhbe', 
               'mock_a549', 'infected_a549')
conditions = rep(conditions, each=3)
conditions = factor(conditions)
sample_table$conditions = conditions

# y ~ x

deseq_dataset = DESeqDataSetFromTximport(txi=count_data, 
                                         colData=sample_table,
                                         design=~conditions)
counts(deseq_dataset)[1:6, 1:3]
count_data$counts[1:6, 1:3]

deseq_dataset = estimateSizeFactors(deseq_dataset)
normalizationFactors(deseq_dataset)
counts(deseq_dataset, normalized=TRUE)[1:6, 1:3]

boxplot(counts(deseq_dataset, normalized=TRUE))
vst = varianceStabilizingTransformation(deseq_dataset)
boxplot(assay(vst))

# observe the cell line effect
plotPCA(vst, intgroup='conditions') +
  theme_bw()

# split experiment in two
dds1 = deseq_dataset[,1:6]
dds2 = deseq_dataset[,7:12]

dds1 = estimateSizeFactors(dds1)
dds2 = estimateSizeFactors(dds2)

# PCA for 2 groups
vst1 = varianceStabilizingTransformation(dds1)
plotPCA(vst1, intgroup='conditions')

vst2 = varianceStabilizingTransformation(dds2)
plotPCA(vst2, intgroup='conditions')

d = assay(vst1)
d = t(d)
d = dist(d)

h = hclust(d)
plot(h)

k = kmeans(t(assay(vst1)), centers=2)
k$cluster

# 3 steps to DESeq2 analysis
# 1) estimate size factors (normalisation)
# 2) estimate dispersions
# 3) apply statistics (Wald Test)


dds1 = estimateDispersions(dds1)

## The above doesn't work, so split experiment
## BEFORE importing data

# sample_files = paste0(
#   pull(sample_table, `Sample Name`), 
#   '/quant.sf')
# names(sample_files) = pull(sample_table, `Sample Name`)
sample_files_nhbe = sample_files[1:6]
count_data_nhbe = tximport(files = sample_files_nhbe,
                      type = "salmon",
                      tx2gene = gene_map,
                      ignoreTxVersion = TRUE)
sample_table_nhbe = sample_table[1:6,]
sample_table_nhbe$conditions = factor(rep(c('mock', 'infected'), each=3), 
                                 levels=c('mock', 'infected'))
dds_nhbe = DESeqDataSetFromTximport(txi=count_data_nhbe, 
                                         colData=sample_table_nhbe,
                                         design=~conditions)
dds_nhbe = estimateSizeFactors(dds_nhbe)
dds_nhbe = estimateDispersions(dds_nhbe)
plotDispEsts(dds_nhbe)

# step 3
dds_nhbe = nbinomWaldTest(dds_nhbe)

## DESeq2 shortcut
# dds_nhbe = DESeq(dds_nhbe)
result_table = results(dds_nhbe)
summary(result_table)
View(as.data.frame(result_table))

# result_table is a DataFrame not a data.frame!
result_df = as.data.frame(result_table)
View(result_df)

plotCounts(dds_nhbe, gene='ENSG00000265794',
           intgroup='conditions')
sum(complete.cases(result_df))

filter_df1 = result_df[complete.cases(result_df),]
View(filter_df1)

# Filter results 
# padj < 0.05
# log2FoldChange > 1 < -1

filter_df1$padj < 0.05

filter_df2 = filter_df1[filter_df1$padj < 0.05,]

abs(filter_df2$log2FoldChange) > 1

filter_df3 = filter_df2[abs(filter_df2$log2FoldChange) > 1,]

View(filter_df3)

plotMA(result_table)

# volcano plot

filter_df1$test = filter_df1$padj < 0.05 & abs(filter_df1$log2FoldChange) > 1

filter_df1 = rownames_to_column(filter_df1, var='ensgene')

g = ggplot(filter_df1, aes(x=log2FoldChange, 
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

listMarts()
ensembl99 = useEnsembl(biomart="ensembl", version=99)
View(listDatasets(ensembl99))
ensembl99 = useDataset("hsapiens_gene_ensembl", 
                       mart=ensembl99)
View(listAttributes(ensembl99))
View(listFilters(ensembl99))
getBM(attributes=c('ensembl_gene_id', 'ensembl_gene_id_version',
                   'ensembl_transcript_id', 'ensembl_transcript_id_version',
                   'external_gene_name'), 
      filters = c('ensembl_gene_id'), 
      values = filter_df1$ensgene[1:6],
      mart = ensembl99)

annotation = getBM(attributes=c('ensembl_gene_id',
                                'chromosome_name',
                                'start_position',
                                'end_position',
                                'strand',
                                'gene_biotype',
                                'external_gene_name',
                                'description'),
                   filters = c('ensembl_gene_id'),
                   values = filter_df1$ensgene,
                   mart = ensembl99)
View(annotation)
View(filter_df1)
annotated_df = left_join(filter_df1, annotation,
                         by=c('ensgene'='ensembl_gene_id'))
View(annotated_df)

g = ggplot(annotated_df, aes(x=log2FoldChange, 
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

anno_df2 = annotated_df[annotated_df$padj < 0.05,]

anno_df3 = anno_df2[abs(anno_df2$log2FoldChange) > 1,]

degs = anno_df3$ensgene

vst_nhbe = varianceStabilizingTransformation(dds_nhbe)
vst_nhbe_mat = assay(vst_nhbe)

data_for_hm = vst_nhbe_mat[degs,]
rownames(data_for_hm) = anno_df3$external_gene_name

heatmap(data_for_hm)

pheatmap(data_for_hm, fontsize_row=4, scale='row')

greys = colorRampPalette(brewer.pal(9, "Greys"))(100)

pheatmap(data_for_hm, fontsize_row=4, scale='row',
         color=greys)

pairs = colorRampPalette(brewer.pal(12, "Paired"))(100)

pheatmap(data_for_hm, fontsize_row=4, scale='row',
         color=pairs)

last_scheme = colorRampPalette(brewer.pal(7, "Blues"))(100)

pheatmap(data_for_hm, fontsize_row=4, scale='row',
         color=last_scheme, cutree_cols = 2,
         cutree_rows = 2)






