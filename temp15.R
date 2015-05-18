# library(devtools)
# install_github("genomicsclass/GSE5859")
library(Biobase)
library(GSE5859)
data(GSE5859)

geneExpression = exprs(e)
sampleInfo = pData(e)

head(sampleInfo)

year = format(sampleInfo$date,"%y")
unique(length(year))

tbl = table(year, sampleInfo$ethnicity)
print(tbl)

x = rowSums(tbl != 0)
sum(x > 1)


month.year = format(sampleInfo$date,"%m%y")
unique(length(month.year))

tbl = table(month.year, sampleInfo$ethnicity)
print(tbl)

x = rowSums(tbl != 0)
mean(x > 1)

library(genefilter)
library(qvalue)

sampleInfo$year = as.factor(format(sampleInfo$date,"%y"))

ind = which(sampleInfo$year %in% c("02","03") & sampleInfo$ethnicity == "CEU")
f = droplevels(sampleInfo$year[ind])
pvals = rowttests(geneExpression[ , ind], f)$p.value
qvals = qvalue(pvals)$qvalues
sum(qvals < 0.05)

ind = which(sampleInfo$year %in% c("03","04") & sampleInfo$ethnicity == "CEU")
f = droplevels(sampleInfo$year[ind])
pvals = rowttests(geneExpression[ , ind], f)$p.value
qvals = qvalue(pvals)$qvalues
sum(qvals < 0.05)

ind = which(sampleInfo$ethnicity %in% c("ASN", "CEU"))
f = droplevels(sampleInfo$ethnicity[ind])
pvals = rowttests(geneExpression[ , ind], f)$p.value
qvals = qvalue(pvals)$qvalues
sum(qvals < 0.05)

ind = which(sampleInfo$ethnicity %in% c("ASN", "CEU") & sampleInfo$year == "05")
f = droplevels(sampleInfo$ethnicity[ind])
pvals = rowttests(geneExpression[ , ind], f)$p.value
qvals = qvalue(pvals)$qvalues
sum(qvals < 0.05)
table(f)

set.seed(3)
ind1 = sample(which(sampleInfo$year=="02" & sampleInfo$ethnicity=="CEU"), size=3)
ind2 = which(sampleInfo$year=="05" & sampleInfo$ethnicity=="ASN")
ind = c(ind1, ind2)
f = droplevels(sampleInfo$ethnicity[ind])
pvals = rowttests(geneExpression[ , ind], f)$p.value
qvals = qvalue(pvals)$qvalues
sum(qvals < 0.05)



