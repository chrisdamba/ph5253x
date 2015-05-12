# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)

data(tissuesGeneExpression)

s = svd(e)

# signflips = sample(c(-1,1), ncol(e), replace=TRUE)
# signflips
# newu= sweep(s$u, 2, signflips, FUN="*")
# newv= sweep(s$v, 2, signflips, FUN="*" )
# all.equal( s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))

m = rowMeans(e)

cor(s$u[,1], m)

newmeans = rnorm(nrow(e))
newe = e + newmeans
sqrt(crossprod(e[,3]-e[,45])) 
sqrt(crossprod(newe[,3]-newe[,45])) 

y = e - rowMeans(e)
s = svd(y)

resid = y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))

all.equal(diag(s$d)%*%t(s$v), s$d * t(s$v))
  

z = s$d * t(s$v)

dim(z)

sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(y[,3]-y[,45]))
sqrt(crossprod(z[,3]-z[,45]))

realdistance = sqrt(crossprod(e[,3]-e[,45]))
approxdistance = sqrt(crossprod(z[1:2,3]-z[1:2,45]))
abs(realdistance - approxdistance)

ks = 1:189
realdistance = sqrt(crossprod(e[,3]-e[,45]))
approxdistances = sapply(ks,function(k){
  sqrt(crossprod(z[1:k,3,drop=FALSE]-z[1:k,45,drop=FALSE] )) 
})
percentdiff = 100*abs(approxdistances - realdistance)/realdistance
plot(ks,percentdiff) ##take a look
min(ks[which(percentdiff < 10)])


realdistances = sqrt(apply(e[,-3]-e[,3],2,crossprod))
approxdistances = sqrt(apply(z[1:2,-3]-z[1:2,3],2,crossprod))
cor(realdistances, approxdistances, method = 'spearman')

