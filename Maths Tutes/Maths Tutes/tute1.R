# Simulate a tree with three leaves and calculate the total tree length.

# first coalescent event exponential with rate 3 choose 2 (scaled units of time, N generations)
# second coalescent event exponential with rate 2 choose 2
# sum to get overall tree length
CoalTimes = rexp(2, rate = c(3,1))
T_3 = sum(CoalTimes)
T_3
# now need to convert to normal time (in generations), which means multiply T_3 by 100000 (N, population)
T_3_unscaled = 100000*T_3
T_3_unscaled
# ^ the simulated coalescent time in generations

###########

# What is expected number of mutations, with rate 10^(-6) mutations/site/generation/individual and 1000 bp/individual.
# It's a Poisson process! Rate = 10^(-6), but 1000 bp per indiv, so * 1000. Convert to scaled time by * 10000. Multiply by appropriate time (scaled) to get expectated number of mutations.
# Need to do twice once the tree splits into two, then three times after the next coalescent event, to get total number of mutations.

ScaMutRate = 10^(-6) * 100000 * 1000
ScaMutRate

ExpMutations = ScaMutRate*3*CoalTimes[1] + ScaMutRate*2*CoalTimes[2]
ExpMutations

# Okay now need to do this 1000 times, then find mean and variance.
Y = c(rep(0,1000))
for (i in 1:1000 ) {
  NewCoalTimes = rexp(2, rate = c(3,1))
  Time = 3*NewCoalTimes[1] + 2*NewCoalTimes[2] # times 3 because 3 lots of T(3) branches, times 2 because 2 lots of T(2) branches (look at a tree)
  Y[i] = rpois(1, lambda = Time * ScaMutRate)
}

mean(Y)
var(Y)

