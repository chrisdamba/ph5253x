library(ggplot2)

# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)

data(tissuesGeneExpression)

s = svd(e - rowMeans(e))

plot(s$d^2/sum(s$d^2)*100, ylab='% variance explained')

z = s$d * t(s$v)

dim(z)

qplot(z[1, ], z[2, ], col = factor(tissue), xlab = "Dimension 1", ylab = "Dimension 2") + theme_bw()

# or

d = dist(t(e))

mds = cmdscale(d)

qplot(mds[, 1], mds[, 2], col = factor(tissue), xlab = "Dimension 1", ylab = "Dimension 2") + theme_bw()


cor(z[1,], mds[, 1])
cor(z[2,], mds[, 2])


###############################################################################

library(GSE5859Subset)

data(GSE5859Subset)

s = svd(geneExpression-rowMeans(geneExpression))

z = s$d * t(s$v)

which.max(cor(sampleInfo$group, t(z)))

max(cor(sampleInfo$group, t(z)))

which.max(cor(sampleInfo$group, t(z))[-1]) + 1

#

month = format(sampleInfo$date, "%m")
month = factor(month)

which.max(cor(as.numeric(month), t(z)))

max(cor(as.numeric(month), t(z)))

table(sampleInfo$g,month)

###############################################################################

df = data.frame(u6 = s$u[,6], chr = geneAnnotation$CHR)
boxplot(u6, chr, data = df)
