# STUFFFFFFFFFFFFFFFFFFFFED

# We are using the same tree as in tute 1 (3 leaf tree). We want the joint posterior distribution of the branch lengths. So we want P(t(3),t(2)) where t(i) is the time to a coalescent event with i lineages.
# We're going to use the ABC formula.

counter = 0
sample_size = 1000

while (counter < sample_size) {

# We need a prior distribution Pi(t(3),t(2)). How can we get this?
# We use the fact that a) we have observed 300 mutations; and b) the mutation rate is 10^2 mutations/10^3 sites/unit of scaled time (10^5 generations)/individual
# From this, we expect that t_total = 300/10^2 = 3, but since t_total = 3*t(3)+2*t(2), we have 3 = 3*t(3)+2*t(2).
# We also know that t(i) > 0, so we can effectively limit 0 < t(i) < 3.
# We can use a uniform distribution with parameters [0,3] and sample t(3) from there.

t3 = runif(1,min=0,max=3)

# Since 3 = 3*t(3)+2*t(2), this implies t(2) = (3-3*t(3))/2

t2 = (3-3*t3)/
  
# So t3, t2 are our sampled parameters from the prior distribution.
  
}

# START AGAIN
# let's just have as our prior, based on expert opinion of mathematicians, that T(3)~Exp(3) and T(2)~Exp(2)

ScaMutRate = 10^(-6) * 100000 * 1000
counter = 1
sample_size = 1000
X = c(rep(0,1000))
Y = c(rep(0,1000))

while (counter < sample_size) {
  CoalTime = rexp(2, rate = c(3,1)) # Sample theta* parameters from our prior.
  TotalTime = 3*CoalTime[1] + 2*CoalTime[2]
  Mutations = rpois(1, lambda = TotalTime * ScaMutRate)
  if (Mutations == 300) {
    X[counter] = CoalTime[1]
    Y[counter] = CoalTime[2]
    CoalTime
    counter = counter + 1
  }
}

mean(X)
mean(Y)

############## stuffed

ScaMutRate = 10^(-6) * 100000 * 1000
counter = 1
sample_size = 1000
X = c(rep(0,1000))
Y = c(rep(0,1000))
newrate1 = 3
newrate2 = 1

while (counter < sample_size) {
  CoalTime = rexp(2, rate = c(newrate1,newrate2)) # Sample theta* parameters from our prior.
  TotalTime = 3*CoalTime[1] + 2*CoalTime[2]
  Mutations = rpois(1, lambda = TotalTime * ScaMutRate)
  if (Mutations == 300) {
    X[counter] = CoalTime[1]
    Y[counter] = CoalTime[2]
    CoalTime
    counter = counter + 1
  }
}

newrate1 = 1/mean(X)
newrate2 = 1/mean(Y)

counter = 1
while (counter < sample_size) {
  CoalTime = rexp(2, rate = c(newrate1,newrate2)) # Sample theta* parameters from our prior.
  TotalTime = 3*CoalTime[1] + 2*CoalTime[2]
  Mutations = rpois(1, lambda = TotalTime * ScaMutRate)
  if (Mutations == 300) {
    X[counter] = CoalTime[1]
    Y[counter] = CoalTime[2]
    CoalTime
    counter = counter + 1
  }
}


