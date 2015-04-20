library(dplyr)
library(magrittr)

# library(devtools)
# install_github("genomicsclass/GSE5859Subset")

library(GSE5859Subset)
data(GSE5859Subset)

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