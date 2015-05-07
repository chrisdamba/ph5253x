# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)

data(tissuesGeneExpression)

# gene expression data
dim(e)
class(e)

# sample info
head(tissue)
table(tissue)

# extract 3 different samples
tissue[c(1,2,87)]

x = e[,1]
y = e[,2]
z = e[,87]

sqrt(sum((x-y)^2))
sqrt(sum((x-z)^2))

sqrt(crossprod(x-y))
sqrt(crossprod(x-z))

# distance between all samples
# the dist() function calculate the distance between all rows
# therefore we transpose the matrix e to calculate the distance between all columns/samples
d = dist(t(e)) # by sample/column
class(d)

as.matrix(d)[1,2] # gives the distance between the first and second column/sample
as.matrix(d)[1,87]

image(as.matrix(d))
