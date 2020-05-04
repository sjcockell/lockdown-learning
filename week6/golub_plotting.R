## Attach required packages ----
library(multtest)
library(ggplot2)

data(golub)

head(golub)
head(golub.gnames)
golub.cl

nrow(golub)
ncol(golub)
dim(golub)
golub.gnames[1042,]
golub[1042,2]
golub[1042,]

golub.cl == 0

golub[1042, golub.cl == 0]

mean(golub[1042, golub.cl == 0])
mean(golub[1042, golub.cl == 1])

# plot CCND3 expression
ccnd3_exp = golub[1042,]
golub_factor = factor(golub.cl, levels=0:1, 
                      labels = c("ALL", "AML"))

df = data.frame('Expression'=ccnd3_exp, 
                'Class'=golub_factor)

g = ggplot(data=df, aes(x=Class, y=Expression))

g1 = g + geom_boxplot()
g2 = g1 + geom_point()

g3 = g + geom_point() + geom_boxplot()

# pick your own colour scheme
# scale_colour_manual(values=c('red', 'blue'))
ggplot(data=df, aes(x=Class, y=Expression)) +
  geom_violin(aes(fill=Class)) +
  geom_jitter(width=0.2) +
  ggtitle('A nice plot') +
  theme_bw()

## Day 27
## Plot 2 genes
golub.gnames[c(829, 1042), 2]

## NOT TIDY - df$CST3_exp = golub[829,]
gene_expression = golub[c(829, 1042),]
## as.vector(gene_expression)
gene_expression = as.vector(t(gene_expression))
# rep(1:5, times=2)
# rep(1:5, each=2)
cancer_type = rep(golub_factor, times=2)
gene_label = rep(c('CST3', 'CCND3'), each=38)

df2 = data.frame(gene_expression, 
                 gene_label, 
                 cancer_type)
ggplot(data=df2, aes(x=cancer_type, 
                     y=gene_expression,
                     colour=cancer_type)) +
  geom_violin(aes(fill=cancer_type), alpha=0.1) +
  geom_jitter(width=0.1) +
  facet_grid(cols = vars(gene_label)) +
  ggtitle('My ultra-professional plot') +
  theme_bw()

