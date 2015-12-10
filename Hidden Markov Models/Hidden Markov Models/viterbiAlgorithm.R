viterbi<-function(A,B,initialProb=rep(1/N,N),obs){
  #Viterbi implements the Viterbi algorithm for a Hidden Markov Model with
  #INPUTS
  ##A is the transition probability matrix
  ##B is the observation probabiliity matrix
  ##initialProb is the vector of initial probabilities for every state, usually denoted \pi
  ##obs is a vector of observed values for any amount of time steps
  #OUTPUTS
  ##the output is a vector containing the most probable sequence of states
  
  #Check out the readme.pdf doc on HMM stuff, this falls under 4 Viterbi Algorithm
  ##Viterbi find a sequence of states that are the most likely to occur given a set of observations and the parameter space \lambda=(A,B,\pi) which are described above
  
  #take the dimensions we're working with and set up matrices right
  if(sum(initialProb)!=1){
    print('Invalid Initial probs, doesnt sum to 1')
    return(0)
  }
  
  N<-dim(B)[1]
  M<-dim(B)[2]
  maxTime<-length(obs)
  delta<-matrix(nrow=N,ncol=maxTime)
  psi<-matrix(nrow=N,ncol=maxTime)
  #set initial values
  delta[,1]<-initialProb*B[,obs[1]] #delta is the prob of seeing the sequence up to a specific point
  psi[,1]<-0 #psi is the state that maximises the probability of getting to that state at that current time step
  #do the damn recursive loop
  for(t in 2:maxTime){
    for(j in 1:N){
      #print(c(t,j))
      maxVal<-max(delta[,t-1]*A[,j])#finds the maximum probability of the previous deltas, gonna work forward from this one
      argMaxPos<-which((delta[,t-1]*A[,j])==maxVal) #This crashes if the max prob repeats
      if(length(argMaxPos)!=1){
        print('Problem: There were repeating maximum probabilities, I dont know what to do with that so I gave up and gave you this nice error message, Sorry!')
        return(0)
      }
      psi[j,t]<-argMaxPos #store this position (which is the state) in psi
      delta[j,t]<-maxVal*B[argMaxPos,obs[t]] #record the probs for this time step
      #print(psi)
      #print(delta)
    }
  }
  #terminate, pick the last values
  maxVal<-max(delta[,maxTime]) #pretty much repeat what happens above except we don't save psi because it's stuck in the returning vector qVec
  argMaxPos<-which((delta[,maxTime])==maxVal) #find position of the maximum prob for the last time step
  p=max(delta[,maxTime]) #this is the probability of the sequence with the highest probability
  qVec<-vector(length=maxTime) #allocate the vector
  qVec[maxTime]<-argMaxPos #store the last state in the vector
  #print(qVec)
  #print(delta)
  
  #backtrack this sequence
  #we gotta work backwards from the only state we know is right (the alst one) back to the beginning
  for(t in (maxTime-1):1){
    qVec[t]<-psi[qVec[t+1],t+1] #yank out the value from the previous (reverse time) (thinking time backwards, this is actually the one that occurs at the time step afterwards in a forward chronological sense...ah fuck thats probably worse, times backwards when we do this stuff so maybe its better to think of it as being the value of the next (forward/chronological) time step)
  }
  return(qVec)
}