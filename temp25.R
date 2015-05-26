# source("http://www.bioconductor.org/biocLite.R")
# biocLite("SpikeInSubset")

library(Biobase)
library(SpikeInSubset)

data(rma95)
y = exprs(rma95)

pData(rma95)

g = as.factor(rep(0:1, each=3))

rownames(y) %in% colnames(pData(rma95))

spike = rownames(y) %in% colnames(pData(rma95))

library(genefilter)

rtts = rowttests(y, g)

dim(rtts)
dim(y)

signif.rtts.idx = rtts$p.value < 0.01

print(sum(signif.rtts.idx & !spike)/sum(signif.rtts.idx))

# volcano plot
mask = abs(rtts$dm) < .2 & rtts$p.value < .01
cols = ifelse(mask, "red", ifelse(spike, "dodgerblue", "black"))
plot(-rtts$dm, -log10(rtts$p.value), cex=.8, pch=16,
              xlim=c(-1,1), ylim=c(0,5),
              xlab="difference in means",
              col=cols)
abline(h=2,v=c(-.2,.2), lty=2)

library(matrixStats)

rsds = rowSds(y[, g==0])

tps.idx = signif.rtts.idx & spike
fps.idx = signif.rtts.idx & !spike
tns.idx = !signif.rtts.idx & !spike
fns.idx = !signif.rtts.idx & spike

tps.sds = rsds[tps.idx]
fps.sds = rsds[fps.idx]
tns.sds = rsds[tns.idx]
fns.sds = rsds[fns.idx]

boxplot(tps.sds, fps.sds, tns.sds, fns.sds, names = c('tps', 'fps', 'tns', 'fns'))

# same as

ind = paste0(as.numeric(spike), as.numeric(signif.rtts.idx))
ind = factor(ind, levels=c("11","01","00","10"), labels=c("TP","FP","TN","FN"))
boxplot(split(rsds, ind))

###############################################################################

library(limma)

fit = lmFit(y, design=model.matrix(~ g))
colnames(coef(fit))
fit = eBayes(fit)

sampleSD = fit$sigma
posteriorSD = sqrt(fit$s2.post)

lim = range(sampleSD, posteriorSD)
plot(sampleSD, posteriorSD, ylim=lim, xlim=lim)
abline(0, 1)
abline(v=sqrt(fit$s2.prior))

###############################################################################

fit.pvals = fit$p.value[, 2]

signif.fit.idx = fit.pvals < 0.01

print(sum(signif.fit.idx & !spike)/sum(signif.fit.idx))

# volcano plot
mask = abs(fit$coef[,2]) < .2 & fit$p.value[,2] < .01
cols = ifelse(mask, "red", ifelse(spike, "dodgerblue", "black"))
plot(fit$coef[, 2], -log10(fit$p.value[, 2]), cex=.8, pch=16,
     xlim=c(-1,1), ylim=c(0,5),
     xlab="difference in means",
     col=cols)
abline(h=2,v=c(-.2,.2), lty=2)
