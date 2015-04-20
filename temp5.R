# source("http://www.bioconductor.org/biocLite.R")
# biocLite("SpikeInSubset")

library(SpikeInSubset)
data(mas133)

e = exprs(mas133)
x = e[,1]
y = e[,2]
plot(x, y, main=paste0("corr=", signif(cor(e[,1],e[,2]), 3)), cex=0.5)
k = 3000
b = 1000 # a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k), col="red", density=0, border="red")

# What proportion of the points are inside the red box?
mean(x<k & y<k)


plot(log2(x),log2(y), main=paste0("corr=",signif(cor(log2(x),log2(y)), 3)), cex=0.5)
k = log2(3000)
b = log2(0.5) # a buffer
polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")


# M-A plot

e = log2(exprs(mas133))
x = e[,1]
y = e[,2]
plot((x+y)/2,y-x, cex=0.5)


e = exprs(mas133)
x = e[,1]
y = e[,2]


sd(log2(x/y))

sum(abs(log2(x/y)) > 1)
