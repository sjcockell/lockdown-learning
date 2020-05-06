# Load packages
library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(DESeq2)

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







