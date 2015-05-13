# library(devtools)
# install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

library(matrixStats)

n = 25

rms = rowMads(geneExpression)
idx = order(rms, decreasing=TRUE)[1:n]

heatmap(geneExpression[idx,])

library(RColorBrewer)
hmcol = colorRampPalette(rev(brewer.pal(11,"RdBu")))(n)

heatmap(geneExpression[idx,], col=hmcol)

###############################################################################

library(gplots)
library(rafalib)
library(lubridate)

cols = brewer.pal(3,"Dark2")
cols = cols[sampleInfo$group + 1]

cbind(colnames(geneExpression), sampleInfo$group, cols)

heatmap.2(geneExpression[idx,],
          col=hmcol,
          trace="none",
          scale="row",
          labRow=geneAnnotation$CHR[idx],
          labCol=month(sampleInfo$date),
          ColSideColors=cols,
          key=FALSE)


###############################################################################

# Create a large data set of random data that is completely independent of sampleInfo$group like this:

set.seed(17)
m = nrow(geneExpression)
n = ncol(geneExpression)
x = matrix(rnorm(m*n),m,n)
g = factor(sampleInfo$group)


library(matrixStats)

n = 50

tts = rowttests(x, g)
tts.idx = order(tts$p.value)[1:n]

sds = rowSds(x)
sds.idx = order(-sds)[1:n]

library(RColorBrewer)
hmcol = colorRampPalette(rev(brewer.pal(11,"RdBu")))(n)

heatmap.2(x[tts.idx,],
          col=hmcol,
          trace="none",
          scale="row",
          labCol=g,
          key=FALSE)

heatmap.2(x[sds.idx,],
          col=hmcol,
          trace="none",
          scale="row",
          labCol=g,
          key=FALSE)

