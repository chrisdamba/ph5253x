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

# Extract group index from sampleInfo
g = sampleInfo$group

# the genefilter package contains an optimized function to perform multiple t-tests
library(genefilter)
ttests = rowttests(geneExpression, factor(g))
pvals = ttests$p.value

# Generate a matrix with the same dims as geneExpression and fill it with random normal variates
m = nrow(geneExpression)
n = ncol(geneExpression)
random.geneExpression = matrix(rnorm(n*m),m,n)

# Perform the same t-test as before for each gene/row in this matrix containing random gene expression data (null)
null.pvals = apply(random.geneExpression, 1, function(e){
  t.test(e[g == 1], e[g == 0], var.equal = TRUE)$p.value
})

alpha = 0.05

sum(pvals < alpha)
sum(null.pvals < alpha)

# volcano plot
plot(ttests$dm, -log10(ttests$p.value), pch=16, col=rgb(0,0,0,0.2)) 

# histogram of p-values
hist(null.pvals)
hist(pvals)

# boxplot of gene expression across samples
boxplot(geneExpression)

# smooth histogram of gene expression across samples
library(rafalib)
shist(geneExpression)

# M-A plot
x = geneExpression[,1]
y = geneExpression[,13]

plot(x,y)
plot((x+y)/2, y-x)


# library(lattice)
# mat = geneExpression
# A = rowMeans(log2(mat))
# M = log2(unlist(mat)) - A
# sample = rep(colnames(mat), each=nrow(mat))
# df = data.frame(M, A, sample, row.names=NULL, check.names=FALSE)
# plt = xyplot(M ~ A | sample, df, panel=panel.smoothScatter)
# plt
