# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)

data(tissuesGeneExpression)

colind = tissue %in% c('colon', 'kidney', 'liver')

mat = e[, colind]

ftissue = factor(tissue[colind])

dim(mat)


###############################################################################
library(ggplot2)

s = svd(mat - rowMeans(mat))

plot(s$d^2/sum(s$d^2)*100, ylab='% variance explained')

z = diag(s$d[1:2]) %*% t(s$v[, 1:2])

dim(z)

qplot(z[1, ], z[2, ], col = ftissue, xlab = "Dimension 1", ylab = "Dimension 2") + theme_bw()

# or

d = dist(t(mat))

mds = cmdscale(d)

qplot(mds[, 1], mds[, 2], col = ftissue, xlab = "Dimension 1", ylab = "Dimension 2") + theme_bw()
