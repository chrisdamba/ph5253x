dbinom(x = 2, size = 4, prob = 0.49)

dbinom(x = 4, size = 10, prob = 0.49)

pbinom(q = 10, size = 20, prob = 0.4, lower.tail = FALSE)

# The probability of winning the lottery is 1 in 175,223,510. If 189,000,000
# randomly generated (with replacement) tickets are sold, what is the
# probability that at least one winning tickets is sold? (give your answer as a
# proportion not percentage)

pbinom(q = 0, size = 189000000, prob = 1/175223510, lower.tail = FALSE)

# or at least two

pbinom(q = 1, size = 189000000, prob = 1/175223510, lower.tail = FALSE)

###############################################################################

size = 20
prob = 0.4
q.lo = 0.35*size
q.hi = 0.45*size

p.lo = pbinom(q = q.lo, size = size, prob = prob)
p.hi = pbinom(q = q.hi, size = size, prob = prob)
p.hi-p.lo

p.lo = pnorm(q = q.lo, mean = size*prob, sd = sqrt(size*prob*(1-prob)))
p.hi = pnorm(q = q.hi, mean = size*prob, sd = sqrt(size*prob*(1-prob)))
p.hi-p.lo

size = 1000
prob = 0.4
q.lo = 0.35*size
q.hi = 0.45*size

p.lo = pbinom(q = q.lo, size = size, prob = prob)
p.hi = pbinom(q = q.hi, size = size, prob = prob)
exact = p.hi-p.lo

p.lo = pnorm(q = q.lo, mean = size*prob, sd = sqrt(size*prob*(1-prob)))
p.hi = pnorm(q = q.hi, mean = size*prob, sd = sqrt(size*prob*(1-prob)))
approx = p.hi-p.lo

abs(exact-approx)

###############################################################################

Ns <- c(5,10,30,100)
ps <- c(0.01,0.10,0.5,0.9,0.99)

par(mfrow = c(length(Ns), length(ps)))
for (N in Ns) {
  k=1:N-1
  for (p in ps) {
    # binomial exact
    exact = dbinom(k,N,p)

    # normal approx
    a = (k+0.5 - N*p)/sqrt(N*p*(1-p))
    b = (k-0.5 - N*p)/sqrt(N*p*(1-p))
    approx = pnorm(a) - pnorm(b)
    
    lim = range(c(approx,exact))
    
    plot(exact, approx, xlim=lim, ylim=lim, col=1, pch=16, main = paste0("N=",N," p=",p))
    abline(0,1)
  }
}

###############################################################################

# exact binomial
N <- 189000000
p <- 1/175223510
dbinom(2,N,p)

# normal approx
a <- (2+0.5 - N*p)/sqrt(N*p*(1-p))
b <- (2-0.5 - N*p)/sqrt(N*p*(1-p))
pnorm(a) - pnorm(b)

# pois approx
dpois(2,N*p)

ppois(1, N*p, lower.tail = FALSE)
