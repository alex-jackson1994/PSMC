simHmm<-function(N,M,A,B,initialProb=rep(1/N,N),maxTime=10000){
  #Simulates a Hidden Markov Model
  #
  #INPUTS
  ##N is the number of possible states
  ##M is the number of possible observations for each state (assumed to be the same for all N)
  ##A is the state transition probability matrix where
  ### A={a_ij}, a_ij=P(at t+1 being at state j | given at t was at state i)
  ##B is the observation symbol probability matrix B={b_j(k)} where
  ### bj(k)=P(observe symbol k at time t | given in state j at time t)
  ##initialProb is the initial probability distribution for each state, set to equally likely for all states by default
  ##maxTime is the number of transitions to simulate, 100 by default
  #
  #OUTPUTS
  ##A list containing the simulated states in a vector as element [[1]] and the simulated observations in a vector as element [[2]]
  q<-vector()
  o<-vector()
  V<-1:M
  states<-1:N
  q[1]<-sample(states,1,prob=initialProb)
  for(t in 1:maxTime){
    o[t]<-sample(V,size=1,replace=TRUE,prob=B[q[t],])
    q[t+1]<-sample(states,size=1,replace=TRUE,prob=A[q[t],])
  }
  return(list(q,o))
}

#####SCRIPTS#####
##Casino Example
N<-2
M<-6
A<-matrix(nrow=2,ncol=2)
A[,1]<-c(0.9,0.1)
A[,2]<-c(0.1,0.9)
B<-matrix(nrow=2,ncol=6)
B[1,]<-rep(1/6,6)
B[2,]<-c(rep(1/10,5),1/2)
qoList<-simHmm(N,M,A,B)
qStates<-qoList[[1]]
obs<-qoList[[2]]
obsFreq<-as.matrix(table(obs))
expectedFreq<-c(rep(2/15,5),1/3)
chisq.test(obsFreq,p=expectedFreq)

##DNA Coding
N<-2 #coding vs noncoding
M<-4 #ACTG
A<-matrix(nrow=N,ncol=N)
A[1,]<-c(0.5,0.5)
A[2,]<-c(0.2,0.8)
B<-matrix(nrow=N,ncol=M)
B[1,]<-rep(1/M,M)
B[2,]<-c(1/3,1/3,1/6,1/6)
qoList<-simHmm(N,M,A,B)
qStates<-qoList[[1]]
obs<-qoList[[2]]
obsFreq<-as.matrix(table(obs))
expectedFreq<-c(13/42,13/42,8/42,8/42)
chisq.test(obsFreq,p=expectedFreq)
print(obsFreq)
