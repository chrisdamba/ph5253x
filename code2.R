# library(devtools)
# install_github("genomicsclass/GSE5859Subset")

library(GSE5859Subset)
data(GSE5859Subset)

dim(geneExpression)
dim(geneAnnotation)
dim(sampleInfo)

head(geneExpression)
head(geneAnnotation)
head(sampleInfo)

# Check if sample filename (used here as sample ID) match columns in geneExpression
identical(colnames(geneExpression), sampleInfo$filename)

# Check if probeset ID match rows in geneExpression
identical(rownames(geneExpression), geneAnnotation$PROBEID)

# Extract group index from sampleInfo
g = sampleInfo$group

# Extract one row from geneExpression (gene expression data for one gene across different samples)
e = geneExpression[25,]

# Perform a t-test to check if the expression of this gene differs between the two groups of samples

# First, check that data from each group are normally distributed
par(mfrow=c(1,2))

qqnorm(e[g == 1])
qqline(e[g == 1])

qqnorm(e[g == 0])
qqline(e[g == 0])

par(mfrow=c(1,1))

# Second, perform the t-test
alpha = 0.05

t.test(e[g == 1], e[g == 0])

# Perform the same t-test as before, but this time for each gene/row (dimension 1 for apply) in geneExpression
pvals = apply(geneExpression, 1, function(e){
  t.test(e[g == 1], e[g == 0], var.equal = TRUE)$p.value
})

# number of differentially expressed genes (without correction for multiple testing)
sum(pvals < alpha)

# Generate a matrix with the same dims as geneExpression and fill it with random normal variates
m = nrow(geneExpression)
n = ncol(geneExpression)
random.geneExpression = matrix(rnorm(n*m),m,n)

# Perform the same t-test as before for each gene/row in this matrix containing random gene expression data (null)
null.pvals = apply(random.geneExpression, 1, function(e){
  t.test(e[g == 1], e[g == 0], var.equal = TRUE)$p.value
})

# number of differentially expressed genes given the null (without correction for multiple testing)
sum(null.pvals < alpha)

# the number of positive tests given the null is the same as the number of tests times alpha
nrow(random.geneExpression) * alpha


# the genefilter package contains a much more optimized function to perform multiple t-tests
library(genefilter)

ttests = rowttests(geneExpression, factor(g))
head(ttests)

plot(ttests$p.value, pvals)
