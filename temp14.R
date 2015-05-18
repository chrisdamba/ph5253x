# library(devtools)
# install_github("genomicsclass/dagdata")
library(dagdata)
data(admissions)

print(admissions)

index = which(admissions$Gender==1)
accepted = sum(admissions$Number[index] * admissions$Percent[index]/100)
applied = sum(admissions$Number[index])
accepted/applied

index = which(admissions$Gender==0)
accepted = sum(admissions$Number[index] * admissions$Percent[index]/100)
applied = sum(admissions$Number[index])
accepted/applied

index = admissions$Gender==1
men = admissions[index,]
women = admissions[!index,]
men.accepted = sum(men$Number*men$Percent/100)
men.not.accepted = sum(men$Number*(1-men$Percent/100))
women.accepted = sum(women$Number*women$Percent/100)
women.not.accepted = sum(women$Number*(1-women$Percent/100))
tab = matrix(c(men.accepted, women.accepted, men.not.accepted, women.not.accepted),2,2)
chisq.test(tab)$p.value

index = admissions$Gender==1
men = admissions[index,]
women = admissions[!index,]
print(data.frame(Major=admissions[1:6, 1], Men=men[, 3], Women=women[, 3]))
print(admissions)

library(dplyr)
library(magrittr)

H = admissions %>% group_by(Major) %>% summarise(Percent = sum(Percent*Number)/sum(Number)) %$% Percent
c('A', 'B', 'C', 'D', 'E', 'F')[which.min(H)]
H[which.min(H)]/100

cor(admissions %>% filter(Gender == 1) %>% select(Number), H)
cor(admissions %>% filter(Gender == 0) %>% select(Number), H)
