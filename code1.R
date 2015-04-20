library(dplyr)
library(magrittr)

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

