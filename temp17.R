# install.packages("devtools")
library(devtools)
# install_github("ririzarr/rafalib")
library(rafalib)
library(Biobase)
# install_github("genomicsclass/GSE5859")
library(GSE5859)
data(GSE5859)

cors <- cor(exprs(e))
Pairs=which(abs(cors)>0.9999,arr.ind=TRUE)
out = Pairs[which(Pairs[,1]<Pairs[,2]),,drop=FALSE]
if(length(out[,2])>0) e=e[,-out[2]]

out <- grep("AFFX",featureNames(e))
e <- e[-out,]

y <- exprs(e)-rowMeans(exprs(e))
dates <- pData(e)$date
eth <- pData(e)$ethnicity

annotation(e)

library(hgfocus.db)
annot <- select(hgfocus.db, keys=featureNames(e), keytype="PROBEID",columns=c("CHR"))
##for genes with multiples, pick on
annot <-annot[match(featureNames(e),annot$PROBEID),]
annot$CHR <- ifelse(is.na(annot$CHR),NA,paste0("chr",annot$CHR))
chryexp<- colMeans(y[which(annot$CHR=="chrY"),])

mypar2()
hist(chryexp)

sex <- factor(ifelse(chryexp<0,"F","M"))


dim(y)

head(eth)
head(sex)

s <- svd(y) # n.b., y has been "demeaned/detrended", i.e., the mean has been removed from each value

pc <- prcomp(y)

# when values in y has been "demeaned/detrended" then svd and prcomp should yield the same principal components

plot(s$v[,4], pc$rotation[,4])

# how many principal components in our data explain much of the variance in the data?

plot(s$d^2/sum(s$d^2))

# what hidden factors/batch effects do these principal components derive from?

# first, we inspect ethnicity

cols = as.numeric(eth)

plot(s$v[,1], s$v[,2], col=cols, pch=16)

# second, we inspect date => year
year = as.factor(format(dates, "%y"))

cols = as.numeric(year)

plot(s$v[,1], s$v[,2], col=cols, pch=16)

yearmonth = as.factor(format(dates, "%y%m"))

boxplot(split(s$v[,1], yearmonth))
boxplot(split(s$v[,2], yearmonth))

# plot the correlation between date => yearmonth and the principal components
plot(t(cor(as.numeric(yearmonth), s$v)), pch=16)

# plot the correlation between sex and the principal components
plot(t(cor(as.numeric(sex), s$v)), pch=16)

boxplot(split(s$v[,7], sex))
boxplot(split(s$v[,8], sex))
