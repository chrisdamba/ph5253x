# p-values are random variables.
# 
# Note that just like the sample average is a random variable because it is
# based on a random sample, p-values are based on random variables (sample mean,
# sample standard deviation) so they are also a random variable.

set.seed(1)
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download.file(url,destfile=filename, method="curl")
population = read.csv(filename)
pvals <- replicate(1000,{
  control = sample(population[,1],12)
  treatment = sample(population[,1],12)
  t.test(treatment,control)$p.val
})
head(pvals)
hist(pvals)

mean(pvals < 0.05)
mean(pvals < 0.01)

###############################################################################

set.seed(100)

signif.pvalues = replicate(1000, (function(){
  pvals = replicate(20, (function(){
    cases = rnorm(10,30,2)
    controls = rnorm(10,30,2)
    t.test(cases,controls)$p.val
  })())
  sum(pvals < 0.05)  
})())

mean(signif.pvalues)

sum(signif.pvalues > 0)/1000
