# install.packages("devtools")
# library(devtools)
# install_github("ririzarr/rafalib")
library(rafalib)
library(Biobase)
# install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

library(sva)

s = svd(geneExpression-rowMeans(geneExpression))
cor(sampleInfo$group,s$v[,1])

sex = sampleInfo$group
mod = model.matrix(~sex)
svafit = sva(geneExpression, mod)
head(svafit$sv)

for(i in 1:ncol(svafit$sv)){
  print( cor(s$v[,i],svafit$sv[,i]) )
}

###############################################################################

X <- model.matrix(~sex+svafit$sv)

pvals = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,4]
})

library(qvalue)
qvals = qvalue(pvals)$qvalue
sum(qvals < 0.1)
mean(geneAnnotation$CHR[qvals < 0.1] %in% c("chrX", "chrY"))

###############################################################################

res = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,c(1,4)]
})

qvals = qvalue(res[2,])$qvalue
pcutoff = max( res[2,qvals < .1] )
library(rafalib)
mypar2(1,1)

plot(res[1,],-log10(res[2,]),xlab="M",ylab="log10 p-value")

ind = which(geneAnnotation$CHR=="chrY")
points(res[1,ind],-log10(res[2,ind]),col=1,pch=16)

ind = which(geneAnnotation$CHR=="chrX")
points(res[1,ind],-log10(res[2,ind]),col=2,pch=16)

abline(h=-log10(pcutoff))
legend("bottomleft",c("chrX","chrY"),col=c(2,1),pch=16)
