#likelihoodTree("AGT","AGG","CGC",0.388,1.287)

likelihoodTree<-function(node1,node2,node3,v1,v3){
  #Calculates the likelihood of a tree in the following form
  #     _|_
  #v2 _|_  | v3
  #v1|  |  | v3
  #  N1 N2 N3
  #
  # REQUIRES:
  #     node1, node2, node3 as strings of their base pairs
  #     v1 and v3 as tree lengths

  pi0<-1/4 #set by equilibrium assumption of jukes-cantor prob model
  v2=v3-v1 #set by tree topology above
  
  #If any of the branch lengths non-positive return likelihood of 0 for hte impossible tree.
  if(v2<=0 || v3<=0 || v1<=0){return(0)}
  
  #Convert strings to character arrays
  node1Array<-strsplit(node1,NULL)[[1]]
  node2Array<-strsplit(node2,NULL)[[1]]
  node3Array<-strsplit(node3,NULL)[[1]]
  dna<-c("A","C","G","T")
  
  #Preallocate vector
  likelihood<-vector(length=length(node1Array))
  
  #loop through each site
  for(l in 1:length(node1Array)){
    
    #double loop through possible dna letters
    for(i in 1:4){
      for(j in 1:4){
        likelihood[l]<-likelihood[l]+jukesCantorProb(dna[i],dna[j],v2)*jukesCantorProb(dna[j],node1Array[l],v1)*jukesCantorProb(dna[j],node2Array[l],v1)*jukesCantorProb(dna[i],node3Array[l],v3)
      }
    }
    likelihood[l]<-likelihood[l]*pi0 #Pi0 times double sum above
  }
  treeLikelihood<-prod(likelihood) #likelihood of whole tree is the product of the likelihoods of each site in the tree
  return(treeLikelihood)
}

jukesCantorProb<-function(node1,node2,v){
  #quick implementation of the probability distribution from the Jukes-Cantor model
  if(node1==node2){
    return(1/4+3/4*exp(-4/3*v))
  }else{
    return(1/4-1/4*exp(-4/3*v))
  }
}