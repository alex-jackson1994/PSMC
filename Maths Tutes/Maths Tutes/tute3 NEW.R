# OBJECTIVE: Using the code from tute 1, use ABC to obtain a sample of size 1000 from the joint exact posterior distribution of the branch lengths for a three node (actually, I'm going to assume three leaf not node) tree which observed 300 mutations, given that mutation rate mu = 10^(-6) (mutations per site per gen per indiv), number of sites l = 10^3 (sites per sequence), and effective population size Ne = 10^5 (individuals).

Mu = 10^(-6)
l = 10^3
EffPopSize = 10^5

SampleSize = 10^3

# 1 unit of scaled time is Ne generations passing.
ScaMuRate = Mu * l * EffPopSize # (mutations per sequence per scaled unit) ??????
ScaMuRate

ObservedMutations = 300 # observation for 3 leaf tree is 300 mutations

###########

# Use ABC:
#   1) Sample theta* from prior P(theta)
#   2) Simulate data x* from f(x|theta*). In this case, x* will be number of mutations
#   3) Usually, we use a summary statistic, but since our observations/data is a scalar we can just compare directly. If obs/data match exactly, accept theta* as coming from the true posterior, else reject.
#   4) Rinse and repeat...

###########

# STEP 1: SAMPLE FROM PRIOR

# The first question is, what is the prior going to be? We want to obtain the joint exact posterior posterior distribution for the branch lengths, so I guess we need a v2-component vector describing the rate parameters of the exponentials. I think (from lecture 1 on coalescent theory) that the joint exact posterior should theoretically be BranchLengths = (Exp(3), Exp(1)). So for for the prior distribution for the rate parameters, what about we have the prior as a bivariate independent normal with mean PrMean = (3,1) and 2x2 variance matrix (1, 0, \\ 0, 1/9). We also need to make sure none of the values there are zero.

# Initialise the posterior sample matrix.
Posterior = matrix(rep(0,2*SampleSize), ncol=2)
Counter = 1

# testing the "while" function
#dummy_max = 100
#tmp = 1
#vec = rep(0,100)
#while (tmp <= dummy_max) {
#  realisation = runif(1)
#  if (realisation > 0.5) {
#    vec[tmp] = realisation
#    tmp = tmp + 1
#  }
#}
#vec

# Defining tree topology and stuff...
Leaves = 3
TimeMultVec = c(Leaves:2)

while (Counter <= SampleSize ) {
  ThetaStar = rnorm(2, mean = c(3,1), sd = c(1,1/9) ) # Simulates Theta* from the prior.
  if (min(ThetaStar) > 0) { # Obviously can't have a negative rate parameter.
    # Simulate a tree with mutations.
    TotalTime = TimeMultVec * rexp(2, ThetaStar) # Simulate tree with Theta* as rate parameters for branch lengths.
    ### POISSON PROCESS: Assigning of mutations. The variable total_time keeps track of the amount of time we have to assign mutations, as if there was only one individual (with sequence length 1000). The exponential rates are already in scaled time.
    ### A poisson process has pmf p(x) = a^x*e^(-a)/x!, where a is the Poisson parameter. Remember, for a Poisson process, the rate will be a = lambda * time, where lambda is the exponential rate parameter. In this case, lambda = ScaMuRate.
    SimulatedMutations = rpois(1, lambda = ScaMuRate * TotalTime)
    if (SimulatedMutations == ObservedMutations) { # If our data matches...
      Posterior[Counter,] = ThetaStar # Accept Theta* as coming from the true posterior
      Counter = Counter + 1 # Onto the next sample...
      }
  }
}
Posterior
Rate_Estimates = c(mean(Posterior[,1]), mean(Posterior[,2]))
Rate_Estimates

hist(Posterior[,1])
hist(Posterior[,2])

# On first go, got (2.455130, 1.005982), which is off from the expected (3, 1). But I guess that makes sense as we're restricting this tree by the number of mutations. Next go, (2.470451, 1.007326).