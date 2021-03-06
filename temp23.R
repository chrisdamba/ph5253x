# Pr(+|D) = 0.99
# Pr(−|no D) = 0.99
# Pr(D) = 0.00025
# 
# Pr(D|+) = [Pr(+|D)*Pr(D)]/Pr(+)
# same as
# Pr(D|+) = [Pr(+|D)*P(D)]/[Pr(+|D)*P(D) + Pr(+|no D)*Pr(no D)]
# same as
# Pr(D|+) = [Pr(+|D)*P(D)]/[Pr(+|D)*P(D) + [1-Pr(-|no D)]*[1- Pr(D)]]

(0.99*0.00025)/((0.99*0.00025)+(1-0.99)*(1-0.00025))
