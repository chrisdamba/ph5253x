# When data are 0s and 1s probabilities and expectations are the same thing!!!
n = 1000
y = rbinom(n,1,0.25)
# proportion of ones Pr(Y)
sum(y==1)/length(y)
# expectaion of Y
mean(y)


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

mean(y[x==176])

plot(y[seq(160,178)], x[seq(160,178)])

idx = sapply(unique(x), function(h){
  mean(y[x==h]) > 0.5
})

max(unique(x)[idx])
