# Doing the MCMC method (Metropolis-Hastings algorithn):

# 1) For t=1, choose starting value Theta_1.
# 2) While t =< Tmax, do
#   3) Sample Theta' from Q(Theta'|Theta_t)
#   4) Calculate Alpha = P(Theta'|X)/P(Theta_t|X)
#   5) If Alpha >= 1
#     6) Theta_t+1 = Theta'
#   7) Else
#     8) Set Theta_t+1 = Theta' with probability Alpha (and if rejected, set Theta_t+1 = Theta_t)
#   9) End if
# 10) t=t+1
# 11) End while

# Doing this for a three leaf tree.

# The genetic data. Remember A,C,G,T=1,2,3,4 respectively.
R1 = c(1,3,4)
R2 = c(1,3,2) 
R3 = c(2,3,2)
Gen = rbind(R1,R2,R3)

#############################
# 1) Initial value
# From earlier lectures, we have shown the branch lengths for the first branch should be distributed Exp(3 choose 2) and the second should be Exp(2 choose 2). So taking those expectations should be a good starting point.
Theta = c(3,1)
ThetaMatrix = matrix(Theta,nrow=1) # we will append further Theta_t guesses onto this matrix.

# 2) Tmax
# Choose it small at first, e.g. 10.
Tmax = 10

# 3) Proposal density Q
# Make it normal. Need to make it non-negative though
SigmaMatrix = matrix(c(1,0,1,0),nrow=2)

# Multivariate normal
# mvrnorm(n = 1, mu, Sigma)

# Set up variables
Negative = FALSE

for (t in 1:Tmax) {
  # 4) Sample ThetaPrime from the proposal density
  ThetaPrime = mvnorm(n=1, Theta, SigmaMatrix) # remember that this is a normal and thus can take negative values (bad)
  if (min(ThetaPrime < 0)) {
    Negative = TRUE
  }
  
  # 5) Calculate Alpha, the ratio of likelihoods.
  
}
