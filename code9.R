# library(devtools)
# install_github("genomicsclass/tissuesGeneExpression")

library(tissuesGeneExpression)
data(tissuesGeneExpression)

# remove placenta samples
table(tissue)
ind = which(tissue != "placenta")
y = tissue[ind]
X = t(e[,ind])
length(y)
dim(X)

# create five folds for k-fold cross-validation
library(caret)
set.seed(1)
idx = createFolds(y, k=5)
sapply(idx, function(ind) table(y[ind]))

# dimension reduction
library(rafalib)
Xsmall = cmdscale(dist(X), k=2)
plot(Xsmall,col=as.fumeric(y))
legend("topleft",levels(factor(y)),fill=seq_along(levels(factor(y))))

# KNN prediction using the first fold as test/validation data
library(class)

k = 5

ind = idx[[1]]

pred = knn(train = Xsmall[-ind, ],
            test = Xsmall[ind, ],
            cl = y[-ind],
            k = k)
table(true=y[ind], pred)
mean(y[ind] != pred)


# k-fold cross-validation to optimize the number of nearest neighbors to use for KNN prediction
set.seed(1)
ks = 1:12
res = sapply(ks, function(k) {  
  res.k = sapply(idx, function(ind) {
    # loop over each of the 5 cross-validation folds
    # predict the held-out samples using k nearest neighbors
    pred = knn(train = Xsmall[-ind, ],
                test = Xsmall[ind, ],
                cl = y[-ind],
                k = k)
    
    # prediction error rate
    mean(y[ind] != pred)
  })
  
  # average error rare over the 5 folds
  mean(res.k)
})

plot(ks, res)
