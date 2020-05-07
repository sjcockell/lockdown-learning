# Load packages
library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(DESeq2)
library(ggplot2)

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
View(result_table)
