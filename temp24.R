tmpfile = tempfile()
tmpdir = tempdir()
download.file("http://seanlahman.com/files/database/lahman-csv_2014-02-14.zip", tmpfile)

filenames = unzip(tmpfile, list=TRUE)
players = read.csv(unzip(tmpfile, files="Batting.csv", exdir=tmpdir), as.is=TRUE)
unlink(tmpdir)
file.remove(tmpfile)

# Here we use dplyr to obtain the necessary information to perform a hierarchical model

library(dplyr)
library(magrittr)

head(players)

filter(players, yearID==2012) %>%
  mutate(AVG=H/AB) %>%
  filter(AB>=500) %>%
  select(AVG)

dat = filter(players, yearID %in% c(2010, 2011, 2012)) %>%
  mutate(AVG=H/AB) %>%
  filter(AB>=500) %>%
  select(AVG)

mean(dat$AVG)  
sd(dat$AVG)

qqnorm(dat$AVG)
qqline(dat$AVG)

# It is April and after 20 at bats, Jose Iglesias is batting .450 (this is very 
# good). We can think of this as a binomial distribution with 20 trials with 
# probability of success p. Our sample estimate of p is .450. What is our 
# estimate of standard deviation? Hint: This AVG is a sum of Bernoulli trials, that 
# is binomial, divided by 20. The sum of Bernoulli trials (the numerator of AVG) is binomial so it has
# SD sqrt(Np(1−p)). The SD of a random variable times a constant is the SD
# of the random variable times that constant so for the AVG we divide by N get
# sqrt(p(1−p)/N)

sqrt(.45*(1-.45)/20)

sqrt(.45*(1-.45)/20)

# The Binomial is approximated by normal with N = 20, so our sampling distribution is
# approximately normal with mean θ and SD σ=0.11. Earlier we used a baseball
# database to determine that our prior distribution for θ is Normal with mean
# μ=0.275 and SD τ=0.027
# 
# E(θ|Y) = Bμ + (1−B)Y = μ + (1−B)(Y−μ)
# where B = σ^2/(σ^2 + τ^2)

Y = 0.450
sigma = 0.11
mu = 0.275
tau = 0.027

B = sigma^2/(sigma^2 + tau^2)

(E.of.theta.given.Y = B*mu + (1-B)*Y)
