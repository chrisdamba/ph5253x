# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)

data(tissuesGeneExpression)

image(e[1:100,])

# for a clustering heatmap, it is customary to plot
# only the top differentially expressed genes
# or
# only the genes that vary a lot

library(genefilter)

rv = rowVars(e)
idx = order(-rv)[1:100]

heatmap(e[idx,])

library(RColorBrewer)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

heatmap(e[idx,], col=hmcol)
