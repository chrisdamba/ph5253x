# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)

data(tissuesGeneExpression)

image(e[1:100,])

# for a clustering heatmap, it is customary to plot
# only the top differentially expressed genes
# or
# only the genes that vary a lot

rvs = rowVars(e)
idx = order(-rvs)[1:100]

heatmap(e[idx,])

library(RColorBrewer)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

heatmap(e[idx,], col=hmcol)

###############################################################################

library(gplots)
library(rafalib)

cols = palette(brewer.pal(7, "Dark2"))[as.fumeric(tissue)]

cbind(colnames(e), tissue, cols)

heatmap.2(e[idx,],
          labCol=tissue,
          trace='none',
          ColSideColors=cols,
          col=hmcol)
