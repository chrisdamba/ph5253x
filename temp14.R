# library(devtools)
# install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

# remove genes on sex chromosomes
table(sampleInfo$group)
y = factor(sampleInfo$group)
X = t(geneExpression)
out = which(geneAnnotation$CHR %in% c("chrX","chrY"))
X = X[, -out]
length(y)
dim(X)

# create five folds for k-fold cross-validation
library(caret)
set.seed(1)
idx = createFolds(y, k=10)
# make sure every fold has 0s and 1s
sapply(idx, function(ind) table(y[ind])) 

idx[[3]][2]

###############################################################################
# KNN prediction using the second fold as test/validation data
library(class)

k=5
m=8

ind = idx[[2]]

# To reduce the number of predictors/dimensions,
# we perform a t-test **on the training data**
# and select the m genes with the smallest p-values

library(genefilter)
p.vals = colttests(X[-ind, ], factor(y[-ind]))$p.val
ind2 = order(p.vals)[1:m]

pred = knn(train = X[-ind, ind2],
            test = X[ind, ind2],
            cl = y[-ind],
            k = k)
table(true=y[ind], pred)
sum(y[ind] != pred)

###############################################################################

res = sapply(idx, function(ind) {
  # loop over each of the 10 cross-validation folds
  # predict the held-out samples using k nearest neighbors
  p.vals = colttests(X[-ind, ], factor(y[-ind]))$p.val
  ind2 = order(p.vals)[1:m]
  
  pred = knn(train = X[-ind, ind2],
             test = X[ind, ind2],
             cl = y[-ind],
             k = k)
  # number of prediction errors
  sum(y[ind] != pred)
})

# average error rare over the 10 folds
sum(res)/length(y)

###############################################################################

# k-fold cross-validation to optimize the number of nearest neighbors to use for KNN prediction
ms=2^c(1:11)
ks=seq(1,9,2)
params = expand.grid(k=ks,m=ms)

res = apply(params, 1, function(param) {
  k = param[1]
  m = param[2]
  
  res.km = sapply(idx, function(ind) {
    # loop over each of the 10 cross-validation folds
    # predict the held-out samples using k nearest neighbors
    p.vals = colttests(X[-ind, ], factor(y[-ind]))$p.val
    ind2 = order(p.vals)[1:m]
    
    pred = knn(train = X[-ind, ind2],
               test = X[ind, ind2],
               cl = y[-ind],
               k = k)
    # number of prediction errors
    sum(y[ind] != pred)
  })
  
  # average error rate over the 10 folds
  sum(res.km)/length(y)
})

params[which.min(res), ]

# make a plot and confirm its just one minimum
res = matrix(res, 5, 11)
matplot(ms, t(res), type="l", log="x", ylab="error rate")
legend("topright", as.character(ks), lty=seq_along(ks), col=seq_along(ks), title = "ks")


###############################################################################

# k-fold cross-validation to optimize the number of nearest neighbors to use for KNN prediction
ms=2^c(1:11)
ks=seq(1,9,2)
params = expand.grid(k=ks,m=ms)

res = apply(params, 1, function(param) {
  k = param[1]
  m = param[2]
  
  # FILTERING BEFORE CROSS-VALIDATION!!!
  # Note how this biases the entire result and gives us much lower estimated
  # error rates. The filtering must be applied without the test set data.
  p.vals = colttests(X, factor(y))$p.val
  ind2 = order(p.vals)[1:m]
  
  res.km = sapply(idx, function(ind) {
    # loop over each of the 10 cross-validation folds
    # predict the held-out samples using k nearest neighbors    
    pred = knn(train = X[-ind, ind2],
               test = X[ind, ind2],
               cl = y[-ind],
               k = k)
    # number of prediction errors
    sum(y[ind] != pred)
  })
  
  # average error rate over the 10 folds
  sum(res.km)/length(y)
})

min(res)

# make a plot and confirm its just one minimum
res = matrix(res, 5, 11)
matplot(ms, t(res), type="l", log="x", ylab="error rate")
legend("topright", as.character(ks), lty=seq_along(ks), col=seq_along(ks), title = "ks")

###############################################################################

y = factor(as.numeric(format( sampleInfo$date, "%m")=="06"))

# Note that we achieve much lower error rate when predicting date than when
# predicting the group. Because group is confounded with date, it is very
# possible that these predictors have no information about group and that our
# lower 0.5 error rates are due to the confounding with date. We will learn more
# about this when discussing batch effects.

# k-fold cross-validation to optimize the number of nearest neighbors to use for KNN prediction
ms=2^c(1:11)
ks=seq(1,9,2)
params = expand.grid(k=ks,m=ms)

res = apply(params, 1, function(param) {
  k = param[1]
  m = param[2]
  
  res.km = sapply(idx, function(ind) {
    # loop over each of the 10 cross-validation folds
    # predict the held-out samples using k nearest neighbors
    p.vals = colttests(X[-ind, ], factor(y[-ind]))$p.val
    ind2 = order(p.vals)[1:m]
    
    pred = knn(train = X[-ind, ind2],
               test = X[ind, ind2],
               cl = y[-ind],
               k = k)
    # number of prediction errors
    sum(y[ind] != pred)
  })
  
  # average error rate over the 10 folds
  sum(res.km)/length(y)
})

min(res)

# make a plot and confirm its just one minimum
res = matrix(res, 5, 11)
matplot(ms, t(res), type="l", log="x", ylab="error rate")
legend("topright", as.character(ks), lty=seq_along(ks), col=seq_along(ks), title = "ks")
