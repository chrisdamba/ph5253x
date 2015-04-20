# source("http://bioconductor.org/biocLite.R")
# biocLite("Biobase")

library(Biobase)

# library(devtools)
# install_github("genomicsclass/GSE5859")

library(GSE5859)
data(GSE5859) # => ExpressionSet e

class(e)

dim(exprs(e))
dim(pData(e))

head(exprs(e)[,1:6])
head(pData(e))

