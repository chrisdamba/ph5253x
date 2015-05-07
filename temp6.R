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

as.matrix(d)[3,45]

(dist(e[c('210486_at', '200805_at'),]))

nrow(e)^2

d = dist(t(e))

?dist

length(d)
ncol(e)^2
