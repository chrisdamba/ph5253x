# install.packages("devtools")
# library(devtools)
# install_github("ririzarr/rafalib")
library(rafalib)
library(Biobase)
# install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

y = geneExpression - rowMeans(geneExpression)

library(lubridate)

sex = sampleInfo$group
month = month(sampleInfo$date)
date = sampleInfo$date

# samples ordered by group/sex
image(cor(y[, order(sex)]))

# samples ordered by month
image(cor(y[, order(month)]))

# samples ordered by date
image(cor(y[, order(date)]))

s = svd(y)

o = order(date)

cols = as.numeric(month)[o]

plot(s$v[o,1], col=cols, pch=16)
plot(s$v[o,2], col=cols, pch=16)

boxplot(s$v[o,1] ~ date[o])
boxplot(s$v[o,2] ~ date[o])

varexplained = s$d^2/sum(s$d^2)
plot(varexplained)
sum(varexplained > 0.1)

cors = cor(month,s$v) # cors between month and cols of s$v
plot(t(cors))
which.max(abs(cors))
max(abs(cors))

cors = cor(sex,s$v) # cors between sex and cols of s$v
plot(t(cors))
which.max(abs(cors))
max(abs(cors))

###############################################################################

X <- model.matrix(~sex+s$v[,1:2])

pvals = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,4]
})

library(qvalue)
qvals = qvalue(pvals)$qvalue
sum(qvals < 0.1)
mean(geneAnnotation$CHR[qvals < 0.1] %in% c("chrX", "chrY"))

