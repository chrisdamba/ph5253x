library(dplyr)
library(magrittr)
library(devtools)
install_github("genomicsclass/GSE5859Subset")

library(GSE5859Subset)
data(GSE5859Subset)

head(geneAnnotation)
head(geneExpression)
head(sampleInfo)

# How many samples where processed on 2005-06-27?
sampleInfo %>%
  filter(date == "2005-06-27") %>%
  summarise(N = n()) %>%
  print

# How many of the genes represented in this particular technology are on chromosome Y?
geneAnnotation %>%
  filter(CHR == "chrY") %>%
  summarise(N = n()) %>%
  print

# What is the log expression value of the for gene ARPC1A on the one subject
# that we measured on 2005-06-10?
filename = sampleInfo %>%
  filter(date == "2005-06-10") %$% filename

probeid = geneAnnotation %>%
  filter(SYMBOL == "ARPC1A") %$% PROBEID

as.data.frame(geneExpression) %>%
  add_rownames() %>%
  filter(rowname == probeid) %>%
  select(rowname, contains(filename)) %>%
  print

g = sampleInfo$group

e = geneExpression[25,]

par(mfrow=c(1,2))

qqnorm(e[g == 1])
qqline(e[g == 1])

qqnorm(e[g == 0])
qqline(e[g == 0])

par(mfrow=c(1,1))

alpha = 0.05

t.test(e[g == 1], e[g == 0])

pvals = apply(geneExpression, 1, function(e){
  t.test(e[g == 1], e[g == 0], var.equal = TRUE)$p.value
})

# number of differentially expressed genes (without correction for multiple testing)
sum(pvals < alpha)

m = nrow(geneExpression)
n = ncol(geneExpression)
random.geneExpression = matrix(rnorm(n*m),m,n)

null.pvals = apply(random.geneExpression, 1, function(e){
  t.test(e[g == 1], e[g == 0], var.equal = TRUE)$p.value
})

sum(null.pvals < alpha)

nrow(random.geneExpression) * alpha

library(genefilter)

ttests = rowttests(geneExpression, factor(g))
head(ttests)

plot(ttests$p.value, pvals)
