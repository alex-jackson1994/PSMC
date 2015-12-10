LikelihoodCalculations = function(SequencesMatrix, BranchLengthsMatrix)


# If I've interpreted the question correctly, what I'm doing here is writing a script which can calculate the likelihood of a set of three sequences, for this particular tree topology. The sequences can be of any (consistent) length.

# I'm naming the nodes 1-5, bottom to top, left to right. That means 1-3 are the leaf nodes, 4 is an internal node, and 5 is the root node.
  #
  #
  #
  #
  #  |   |   |
  #  1   2   3

##################################

# Define the mutation probability function Mut(I,F,T) (with I as initial NT, F as final NT, and T as time), according to the Jukes-Cantor model.

Mut = function(I,F,T){
  # Two options:
  # If we "mutate" to the same nucleotide (I==F), we use P=0.25+0.75*exp(-4/3*T)
  # If we mutate to a different nucleotide (I=/=F), we use P=0.25-0.25*exp(-4/3*T)
  if (I==F){
    return(0.25+0.75*exp(-4/3*T))
  } else {
    return(0.25-0.25*exp(-4/3*T))
  }
}

##################################

# Store the genetic information. I'm using A,C,G,T = 1,2,3,4. Define the Gen matrix to be the known genetic sequences in the rows. Leaf 1 (AGT=1,3,4) is row 1, leaf 2 is row 2, leaf 3 is row 3.

#R1 = c(1,3,4)
#R2 = c(1,3,2)
#R3 = c(2,3,2)

# R1 = c(1,1,1)
# R2 = c(1,1,1)
# R3 = c(1,1,1)

Gen = SequencesMatrix

# Store the branch lengths (times). Define the V matrix to contain the branch lengths.
# We have a 5-node tree.
nodes = 5
V = matrix(c(rep(0,nodes^2)),nodes,nodes)
# I've just made an empty 5x5 matrix. Now I'll add the relevant values in manually. Position V[i,j] is the time between node i and j.
V[5,4] = 0.899
V[4,1] = 0.388
V[4,2] = 0.388
V[5,3] = 1.287

##################################

# Start the clock!
ptm <- proc.time()




# Got the preliminaries done, now need to write the function that gets the likelihood. Because of the independence between sites assumption, we can consider each site separately. For example, first consider the site 1 tree with A,A,C at nodes 1,2,3 and find its likelihood. Then consider the site 2 tree with G,G,G at nodes 1,2,3. Then we multiply the individual likelihoods together to get the overall tree likelihood.

# First, we need to figure out how long the sequences are (trivially 3 in this case).
seq.length = dim(Gen)[2]

# Second way of doing it. FASTER!

ptm <- proc.time()


L2 = c(rep(0,seq.length))
# This loop in j considers each of the sites in turn.
for (j in 1:seq.length) {
  for (i in 1:4) { # This loop in i considers each of the possibilities (A,C,G,T = 1,2,3,4) for the root node (number 5).
    # This does things a different way. The counter A is our sum over the internal node 4, which increases each time we consider a different nucleotide possibilty. Hopefully this makes sense...
    A = 0
    (for (k in 1:4){ # This loop in k considers each of the possibilities (A,C,G,T = 1,2,3,4) for the internal node (number 4). This doesn't seem like the best way to do it, but I'm not sure how to do it better >_<
      A = A + # 
        Mut(i,k,V[5,4]) * # Probability of going from node 5 to node 4 over time V[5,4] = 0.899.
        Mut(k,Gen[1,j],V[4,1]) * # Probability of going from node 4 to node 1 over time V[4,1] = 0.388.
        Mut(k,Gen[2,j],V[4,2]) # Probability of going from node 4 to node 2 over time V[4,2] = 0.388.
    })
    L2[j] = L2[j] + 0.25 * A * Mut(i,Gen[3,j],V[5,3]) # 0.25 is the equilibrium probability of each of A,C,G,T under Jukes-Cantor assumptions.
    
  }
}

L2
TreeLikelihood2 = prod(L)
TreeLikelihood2

proc.time() - ptm
