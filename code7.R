# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)

data(tissuesGeneExpression)

d = dist(t(e))

hc = hclust(d)

plot(hc, cex=0.5, label=tissue)

library(rafalib)

myplclust(hc, cex=0.5, label=tissue, lab.col=as.fumeric(tissue))

abline(h=120)

cl = cutree(hc, h=120) # => clusters; alternatively, number of clusters => h

table(tissue, cluster=cl)

###############################################################################

km = kmeans(t(e), centers=7)

table(tissue, cluster=km$cluster)

d = dist(t(e))
mds = cmdscale(d)
qplot(mds[, 1], mds[, 2], col = km$cluster, xlab = "Dimension 1", ylab = "Dimension 2") + theme_bw()

