# library(devtools)
# install_github("genomicsclass/dagdata")

library(dagdata)
data(hcmv)

library(rafalib)
mypar2()
plot(locations,rep(1,length(locations)),ylab="",yaxt="n")

breaks = seq(0, 4000*round(max(locations)/4000) , 4000)
tmp = cut(locations, breaks)
counts= as.numeric(table(tmp))
counts

hist(counts)

logprobs = dpois(counts, 4, log=TRUE)
loglikelihood = sum(logprobs)
loglikelihood

pois.loglikelihood <- function(lambda, counts) {
  logprobs = dpois(counts, lambda, log=TRUE)
  loglikelihood = sum(logprobs)
  loglikelihood
}

lambdas = seq(0, 15, len=300)

pois.loglikelihoods = sapply(lambdas, function(lambda) pois.loglikelihood(lambda, counts))

plot(pois.loglikelihoods ~ lambdas)

mle = lambdas[which.max(pois.loglikelihoods)]

abline(v=mle)

print(mle)

mle = mean(counts)
mle

###############################################################################

breaks = seq(0, 4000*round(max(locations)/4000) , 4000)
tmp = cut(locations, breaks)
counts= as.numeric(table(tmp))
counts

binCenters = (breaks[-1] + breaks[-length(breaks)])/2
plot(binCenters, counts, type="l", xlab=)
binCenters[which.max(counts)]
max(counts)

lambda = mean(counts[-which.max(counts)])
pval = ppois(13, lambda, lower.tail = FALSE)

pval < 0.05/57

qqplot(qpois(ppoints(length(counts)), mle), counts)
qqline(counts, distribution = function(p) qpois(p, mle))
