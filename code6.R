# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)

data(tissuesGeneExpression)

set.seed(1)

# take a random sample of the original data set to make computations faster for this demo!
ind = sample(nrow(e), 500)
Y = e[ind,]

# standardize each row to make explanations simpler (not necessary)!
Y = t(apply(Y,1,scale))

dim(Y) # -> 500 (features/genes) x 189 (samples)

head(Y)

s = svd(Y)

U = s$u
V = s$v
D = diag(s$d) # s$d is a diagonal 'matrix', a vector containing the elements on the diagonal

plot(s$d)

Yhat = U %*% D %*% t(V)
resid = Y - Yhat
boxplot(resid, ylim=c(-2,2))

k = ncol(Y) - 4

Yhat = U[, 1:k] %*% D[1:k, 1:k] %*% t(V[, 1:k])
resid = Y - Yhat
boxplot(resid, ylim=c(-2,2))

plot(s$d^2/sum(s$d^2)*100, ylab='% variance explained')

k = ncol(Y) - 94

Yhat = U[, 1:k] %*% D[1:k, 1:k] %*% t(V[, 1:k])
resid = Y - Yhat
boxplot(resid, ylim=c(-2,2))

# % variance explained by the columns/samples that were removed
var(as.vector(resid)/var(as.vector(Y)))
1 - sum(s$d[1:k]^2)/sum(s$d^2)

dim(Yhat)
