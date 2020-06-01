# attach packages
library(tidyverse)
library(VennDiagram)
library(pheatmap)

source('covid-rna-analysis.R', echo = TRUE)
source('covid-a549-analysis.R', echo = TRUE)

#nhbe 1st filter results
#annotated_df
#a549 first filter results
#a549_anno_df

nhbe_ensgene = annotated_df$ensgene
a549_ensgene = a549_anno_df$ensgene

common_ensgene = intersect(nhbe_ensgene, 
                           a549_ensgene)
# union
# setdiff

nhbe_common = annotated_df[annotated_df$ensgene %in% common_ensgene, c('ensgene', 'log2FoldChange')]
a549_common = a549_anno_df[a549_anno_df$ensgene %in% common_ensgene, c('ensgene', 'log2FoldChange')]

common_fold_changes = left_join(nhbe_common,
                                a549_common, 
                                by=c('ensgene'='ensgene'))

ggplot(data=common_fold_changes, 
       aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_point(alpha=0.3) +
  geom_abline(slope=1, intercept=0, color='red') + 
  xlim(min=-5, max=5) +
  ylim(min=-5, max=5)

cor.test(common_fold_changes$log2FoldChange.x,
         common_fold_changes$log2FoldChange.y,
         method = 'spearman')

cor.test(common_fold_changes$log2FoldChange.x,
         common_fold_changes$log2FoldChange.y,
         method = 'pearson')

nhbe_diff = anno_df2$ensgene
a549_diff = a549_anno2$ensgene

common_diff = intersect(nhbe_diff, a549_diff)
nhbe_only = setdiff(nhbe_diff, a549_diff)
a549_only = setdiff(a549_diff, nhbe_diff)

venn.diagram(x=list(nhbe_diff, a549_diff))
plot.new()
draw.pairwise.venn(area1 = length(nhbe_diff),
                   area2 = length(a549_diff),
                   cross.area = length(common_diff),
                   scaled=TRUE, fill=c('red','blue'),
                   alpha=0.5)

union_diff = union(nhbe_diff, a549_diff)
union_fc_hm = filter(common_fold_changes, ensgene %in% union_diff)

union_fc_hm_mat = as.matrix(union_fc_hm[,2:3])

default_colours = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                    "RdYlBu")))(100)
my_colours = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                    "RdBu")))(100)

my_breaks = c(seq(-2, -0.01, length.out=50),
              0,
              seq(0.01, 2, length.out=50))
                 
pheatmap(union_fc_hm_mat, breaks = my_breaks, color = my_colours)
