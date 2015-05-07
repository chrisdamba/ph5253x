# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)

data(tissuesGeneExpression)

s = svd(e)
signflips = sample(c(-1,1), ncol(e), replace=TRUE)
signflips

newu= sweep(s$u,2,signflips,FUN="*")
newv= sweep(s$v,2,signflips,FUN="*" )
all.equal( s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))

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

x=matrix(rep(c(1,2),each=5),5,2)
x

x*c(1:5)

sweep(x,1,1:5,"*")

all.equal(diag(s$d)%*%t(s$v), s$d * t(s$v))
  
vd = t(s$d * t(s$v))
