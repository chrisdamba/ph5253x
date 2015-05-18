library(GSE5859Subset)
data(GSE5859Subset)

# Note that sampleInfo$group here represents males and females. Thus we expect
# differences to be on chrY and, for genes that escape inactivation, chrX. Note
# that we do not expect many autosomal genes to be different between males and
# females. This gives us an opportunity to evaluate false and true positives
# with experimental data.

sex = as.factor(sampleInfo$group)
library(lubridate)
month = as.factor(month(sampleInfo$date))
table(sex, month)

library(genefilter)
library(qvalue)
pvals = rowttests(geneExpression, sex)$p.value
qvals = qvalue(pvals)$qvalues
sum(qvals < 0.1)

###############################################################################

sum(geneAnnotation$CHR %in% c("chrX", "chrY") & qvals < 0.1)/sum(qvals < 0.1)
# same as
mean(geneAnnotation$CHR[qvals < 0.1] %in% c("chrX", "chrY"))

###############################################################################

ind = which(!(geneAnnotation$CHR %in% c("chrX", "chrY")) & qvals < 0.1)
pvals = rowttests(geneExpression[ind, ], month)$p.value
mean(pvals < 0.05)

###############################################################################

# The above result shows that the great majority of the autosomal genes show
# differences due to processing date (month). This provides further evidence that
# confounding is resulting in false positives. So we are going to try to model
# the month effect to better estimate the sex effect. We are going to use a
# linear model:

X = model.matrix(~ sex + month)

pvals = sapply(1:nrow(geneExpression), function(i){
  y = geneExpression[i,]
  fit = lm(y ~ X-1)
  return(summary(fit)$coef[2,4])
})

qvals = qvalue(pvals)$qvalue
sum(qvals < 0.1)

sum(geneAnnotation$CHR %in% c("chrX", "chrY") & qvals < 0.1)/sum(qvals < 0.1)
# same as
mean(geneAnnotation$CHR[qvals < 0.1] %in% c("chrX", "chrY"))

pvals = sapply(1:nrow(geneExpression), function(i){
  y = geneExpression[i,]
  fit = lm(y ~ X-1)
  return(summary(fit)$coef[3,4])
})

qvals = qvalue(pvals)$qvalue
sum(qvals < 0.1)
