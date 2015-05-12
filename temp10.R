# library(devtools)
# install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

set.seed(10)

km = kmeans(t(geneExpression), centers=5)

library(lubridate)

table(group=sampleInfo$group, cluster=km$cluster)
table(group=sampleInfo$date, cluster=km$cluster)
table(group=year(sampleInfo$date), cluster=km$cluster)

d = dist(t(geneExpression))
mds = cmdscale(d)
plot(mds[, 1], mds[, 2], col = km$cluster, xlab = "Dimension 1", ylab = "Dimension 2")

