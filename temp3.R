# p-values are random variables.
# 
# Note that just like the sample average is a random variable because it is
# based on a random sample, p-values are based on random variables (sample mean,
# sample standard deviation) so they are also a random variable.

# the probability of two random events that are statistically independent occurring is
# P(A and B) = P(A) * P(B). Note that this is a consequence of the more general formula
# P(A and B) = P(A) P(B | A )

# Sidak error controlling procedure vs Bonferroni error controlloing procedure

alphas <- seq(0,0.25,0.01)

par(mfrow=c(2,2))

for(m in c(2,10,100,1000)){
  plot(alphas, 1-(1-alphas)^(1/m), type="l", main=paste("m =", m), ylab="k", col="red")
  lines(alphas, alphas/m, col="blue")
  legend("topleft", legend=c("sidak", "bonferroni"), fill = c("red", "blue"))
}

par(mfrow=c(1,1))

###############################################################################

set.seed(1)

alpha = 0.05
m = 8793
k.bonferroni = alpha/m
k.bonferroni

vs = replicate(10000, (function(){  
  pvals = runif(m,0,1)
  v = sum(pvals < k.bonferroni)
  v
})())

fwer.bonferroni = mean(vs >= 1)
fwer.bonferroni

k.sidak = 1-(1-alpha)^(1/m)
k.sidak

vs = replicate(10000, (function(){  
  pvals = runif(m,0,1)
  v = sum(pvals < sidak)
  v
})())

fwer.sidak = mean(vs >= 1)
fwer.sidak

