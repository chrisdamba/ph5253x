n = 10000
set.seed(1)
men = rnorm(n,176,7) # height in centimeters
women = rnorm(n,162,7) # height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))

# mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]

set.seed(5)
N = 250
ind = sample(length(y),N)
Y = y[ind]
X = x[ind]

fit = loess(Y ~ X)
predict(fit, 168)

set.seed(5)
prs = replicate(1000, (function(){
  N = 250
  ind = sample(length(y),N)
  Y = y[ind]
  X = x[ind]
  fit = loess(Y ~ X)
  predict(fit, 168)
})())

library(rafalib)
popsd(prs)
