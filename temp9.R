set.seed(1)
m = 10000
n = 24
x = matrix(rnorm(m*n),m,n)
colnames(x)=1:n

d = dist(t(x))

hc = hclust(d)

plot(hc, cex=0.5)

abline(h=143)

cl = cutree(hc, h=143) # => k = number of clusters; alternatively, k => h = height of cut

length(unique(cl))

###############################################################################

set.seed(1)

cuts = replicate(n = 100, expr = (function(){
  m = 10000
  n = 24
  x = matrix(rnorm(m*n),m,n)
  colnames(x)=1:n  
  d = dist(t(x))  
  hc = hclust(d)
  cl = cutree(hc, h=143)
  return(length(unique(cl)))
})())

plot(table(cuts))
popsd(cuts)
