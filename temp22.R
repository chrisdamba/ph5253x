# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)

data("tissuesGeneExpression")

y = e[, which(tissue=="endometrium")]

library(genefilter)

head(y)

mypar(1,2)
vars = rowVars(y)
qqnorm(vars)
qqline(vars)

sds = rowSds(y)
qqnorm(sds)
qqline(sds)

###############################################################################

library(limma)

df1 = ncol(y)-1

fit = fitFDist(y, df1)

scale = fit$scale
df2 = fit$df2

mypar(1,1)
qqplot(qf(ppoints(length(counts)), df1, df2), y*scale)
qqline(y*scale, distribution = function(p) qf(p, df1, df2))
