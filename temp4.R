library(GSE5859Subset)
data(GSE5859Subset)

library(genefilter)

g = sampleInfo$group

ttests = rowttests(geneExpression, factor(g))
head(ttests)

pvals = ttests$p.value

m = length(pvals)
alpha = 0.05

# empirical distribution of p-values
hist(pvals, ylim=c(0,1500), xlab="pval", main=NULL,
     col=rgb(1,1,0,0.4))

par(new=TRUE) # to overlay the next plot

# null distribution of p-values
hist(runif(m, 0, 1), ylim=c(0,1500), xlab=NULL, ylab=NULL, main=NULL,
     col=rgb(1,1,0,0.6))


# How many genes have p-values smaller than alpha = 0.05?
sum(pvals < alpha)

# Apply the Bonferroni correction to the p-values to achieve a FWER of 0.05.
# How many genes are called significant under this procedure?

k = alpha/m
sum(pvals < k)

# same as

pvals.bonferroni = p.adjust(pvals, method = "bonferroni")
sum(pvals.bonferroni < 0.05)

# Apply the 'FDR' correction to the p-values to achieve a FDR of 0.05.
# How many genes are called significant under this procedure?
pvals.fdr = p.adjust(pvals, method = "fdr")
sum(pvals.fdr < 0.05)


# Use the qvalue function, in the Bioconductor qvalue package, to estimate q-values using the procedure described by Storey.
# Using this estimate how many genes have q-values below 0.05?


# source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue")
library(qvalue)

res = qvalue(pvals)
qvals = res$qvalues
sum(qvals < 0.05)
res$pi0

###############################################################################

# Create a Monte Carlo Simulation in which you simulate measurements from 8,793 genes
# for 24 samples: 12 cases and 12 controls

experiment = function() {
  n = 24
  m = 8793
  mat = matrix(rnorm(n*m),m,n)  

  # Now for the first 500 genes, we add a difference of 2 between cases and controls
  delta = 2
  positives = 500
  mat[1:positives,1:(n/2)] = mat[1:positives,1:(n/2)] + delta
  
  # So the null hypothesis is true for 8793-500 genes, i.e., m=8793, m0=8293 and m1=500
  m0 = 8793
  m1 = 500
  
  # Compute p-values using a t-test (using rowttests in the genefilter package)  
  ttests = rowttests(mat, factor(g))
  pvals = ttests$p.value
  
  # Create three lists of genes using:  
  # 1) Bonferroni correction to achieve an FWER of 0.05
  pvals.bonferroni = p.adjust(pvals, 'bonferroni')
  
  # 2) p-adjust estimates of FDR to achieve an FDR of 0.05 
  pvals.fdr = p.adjust(pvals, 'fdr')
  
  # 3) qvalue estimates of FDR to to achieve an FDR of 0.05
  res = qvalue(pvals)
  qvals = res$qvalues
  
  k = 0.05
    
  # For each of these three lists compute the number of false positives in the list and the number
  # of false negatives: genes not in the list that should have been because the null hypothesis is not true

  pvals.bonferroni.fps = sum(pvals.bonferroni[-(1:positives)] < k)
  pvals.bonferroni.fns = sum(pvals.bonferroni[1:positives] >= k)

  pvals.fdr.fps = sum(pvals.fdr[-(1:positives)] < k)
  pvals.fdr.fns = sum(pvals.fdr[1:positives] >= k)
  
  qvals.fps = sum(qvals[-(1:positives)] < k)
  qvals.fns = sum(qvals[1:positives] >= k)
  
  # return the false positive rate = false positive / m0 and the false negative rate = false negative / m1
  
  return(c(pvals.bonferroni.fps/m0, pvals.bonferroni.fns/m1,
                             pvals.fdr.fps/m0, pvals.fdr.fns/m1,
                             qvals.fps/m0, qvals.fns/m1))
}

set.seed(1)

experiments = replicate(1000, experiment())

rates = rowMeans(experiments)
rates
