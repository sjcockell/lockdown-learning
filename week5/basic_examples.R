# Fundamentals of R

## Calculations ----
2 + 1
2 + 5 * 7  # 2 + (5 * 7) (order of operations honoured)

## variable assignment ----
x = 10
# x <- 10 will also work

## functions ----
log(10)
log(10, base=2)  #log2
log(10, base=10) #log10
?log
x = log2(16)
y = x^2
x + y
x = log2(x)

## data types ----
x = 3                 # double
y = "Bioinformatics"  # character
z = FALSE             # logical

is.double(x)
is.double(y)
is.character(y)

# special data types
2/0                   # Inf
Inf - Inf             # NaN
NA                    # NA

# this is a comment
# Vectors ----

# c() to construct
my_first_vector = c(2, 3, 5, 6, 2)
is.double(my_first_vector)
is.vector(my_first_vector)
my_second_vector = c("X", "Y", "M")
is.character(my_second_vector)

my_third_vector = c(2, 3, "a", "b")
is.double(my_third_vector)
is.character(my_third_vector)

# logical < double < character

# Making sequences of doubles
seq(1, 6)
?seq
seq(1, 6, by=2)
seq(6, 1, by=-1)
1:6
?':'
6:1

myvec = c(1, 2, 3, 5, 8)
myvec2 = c(6, 9, 12, 15, 18)

# vector arithmetic
myvec + myvec2
myvec * myvec2

# vector recycling
myvec = c(8, 11)
myvec2 = c(2, 3, 5, 7)

myvec + myvec2

## Day 23 ----

x = c(6, 4, 3, 8, 1, 3, 10, 5)

## Vector functions ----
length(x)
unique(x)

which(x == 3)

rev(x)

sort(x)

sum(x)
mean(x)
median(x)
quantile(x)
summary(x)

# nesting functions
length(unique(x))

paste0('Simon', 'Cockell')
paste('Simon', 'Cockell')

paste0('chr', seq(1, 5))

## Subsetting vectors ----
x[5]
x[c(1,5)]

vec = c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE)

x[vec]

## Relational operators ----

10 == 10
10 == 9
'cat' == 'cat'
'cat' == 'dog'
10 != 10
10 != 9
'cat' != 'cat'
'cat' != 'dog'
10 > 10
10 > 9
'cat' > 'cat'
'cat' > 'dog'
10 < 10
10 < 9
10 >= 10
10 <= 10

## Subsetting using relational operators ----
x < 5
x[x < 5]
x[x > 'cat']
length(x[x==3])

x %in% c(2, 3, 5)
x[x %in% c(2, 3, 5)]
x[x %in% seq(2, 10, by=2)]
x[x %in% seq(0, max(x), by=2)]


1 < 'a'
'a' < 'A'

9 %% 4

x[x %% 2 == 0]

## Day 24 ---- 

## Factors ----
# factor() to construct
chr_vec = c("chr3", "chrX", "chr2",
            "chr1", "chr20", "chr10")
sort(chr_vec)

ordered_vec = paste0("chr", c(1:22, 'X', 'Y', 'M'))
# construct factor with specified sorting order
chr_fac = factor(chr_vec, 
                 levels=ordered_vec, 
                 ordered=TRUE)
sort(chr_fac)

## Matrices ----
# matrix() to construct
my_first_matrix = matrix(1:12, nrow = 4, ncol = 3)

matrix(1:12, nrow = 4, ncol = 3, byrow=TRUE)
matrix(c(1:5, 'a'), nrow = 4, ncol = 3)

## matrix functions
length(x)
length(my_first_matrix)
dim(my_first_matrix)
t(my_first_matrix)

as.vector(my_first_matrix)
t_matrix = t(my_first_matrix)

as.vector(t_matrix)

my_second_vector = matrix(1:12, nrow = 4, ncol = 3, byrow=TRUE)
as.vector(my_second_vector)

# subsetting matrices
my_first_matrix[3:4, 2:3]

my_first_matrix[c(TRUE, FALSE, TRUE, FALSE),]

## Lists ----
# list() to construct
gene = list(name="GAPDH", 
            protein_coding=TRUE,
            chromosome=13)

x[1]
# subsetting lists
gene[[1]]
gene[["name"]]
gene$name

## Data Frames ----
df = data.frame(gene_symbol = c("TP53", "GAPDH", "FOSB1"),
                fold_change = c(2.1, 0.45, 1.3),
                p_value = c(0.02, 0.04, 0.23),
                significant = c(TRUE, TRUE, FALSE),
                stringsAsFactors=FALSE)

# subsetting data.frames
df[["fold_change"]]
df[[2]]
df$gene_symbol
