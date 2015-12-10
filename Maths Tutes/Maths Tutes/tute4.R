# TUTORIAL 4 - MARKOV CHAIN MONTE CARLO

# OBJECTIVE: Given the sequences in the sequences.txt file, use an MCMC approach to find the joint posterior distribution for the branch lengths of the three leaf tree (remember to rescale the branch lengths).

# First step: how to read in the data.
# Um... just copy/paste it in I guess!
s_1 = "ttgggacatatcttttcagtaaaagctttactgatcttccactgcgccctgcttttttgacatagcagaaattgtataacgattccacccggaagattcgacgacctatcgaagaggagtgactcgaacggcgcaggtaagcctagctactaatcgtcaatgatagctaacattctcaggcctccaggtgtcccttctcgaataagcttatgatgcaggcatcggtgcgagcattcgggctgaaagacgtttgacgatccccacttccaatggtttcgacaaggggtatgactaacggcgtatggatcgttaattctgtcgtctctacgcatatcaacggtttctggatggacgggtggcggggcactcgttagggggtgaggctacaggcgcgctgagaatgctgtatataatcgcctgaaccgtcgtgggcggccgtctggccgccccccgcaggctttgacttgcggggacttgagtccggcagcggtgaggaggcgcagccactagccctggatcagttcttgcattagggcagaacacaccctcttgtttaatccgatcagtcttgccctctggcggactgcaatagccattgtctggataaaacgaaagattggaggatataaccttctcattcctgtcgactgattcggcgttaacgaattgtagcgcactcgcctacgcgtgcaaatacaaccggttacgaacgtcttacccaatcactcaacctgaacctcccggtttttgcctttcacaccctgattacgttgacaaagagggaaggccccaggtagactctcatactgaatctcgccggttaaaccagcacggtcccctttaactcctcccccctcgtgagacgaggctgagtgcttgaaatgatgcaccattggaaagtccaccaaaacccatgtggtacccatactctaggagatttctacgcccgctgggatctaagaattggagcagtaccgctaccactattttcagtcagtga"
s_2 = "ttgggacatatcgtttcagtaaaagctttactaatcttgcactgcgcccggcttatttgacatagcagaaattgtataacgattccaccaggaatattcgacgacctatcgaagatgagtgactcgaacggcgcaggtaagcctagctactaatcgtccatgatagctaacattctcaggcctccaggtgtcccttctcgaataagcttatgatgcaggcatcggtgggagcattggggctgaaagacgtttgacgatccccatttcctatggtttcgacaaggggtatgactaacggcgcatggatcgttaattctgtcgtctctccgcatatcaacggagtccggatggacgggtgccggggtactcgttagggggtaaggctacaggcgcgctgagaatgctgtatataatcgcctgaaccgtcgtgggcggccgtctggccgccccccgcagcctttgacttgcggggacttgagtccggcagcggtgaggaggcgcagccactagccctggatcagttcttgcattaggggagaacacaccctcttgtttaatccgatcagtcttgtcctctggcggactgcgatagcgattgtctggataaaacgaaagattggaggatataactttttcattgctgtcgactgattgggctttaacgaattgtagcgcactcgcctactcgtgcaaaaacaaccggttacgaacgtctaacccaatcactcaacctcaacctcccggtttttgcctttcactccctgattacgttgacaaagagggaaggtcccaggtagactctcatactgaatctcgccagttaaaccagcccggtccactgtaactcctcccccctcgtgagacgaggctgcgtgcttgaaattatgaaccattagaaagtccactaaaacccatgtggtacccatactctaggagatttctacgcccggtgggatctatgaattggagcagtaccgctaccactattttcagtccgtga"
s_3 = "ttaagacatatctattttgcaaaagctttactcatattccactgcgcgcgggctttttgcaacaacagaaactttataaaaattcgaccaggaacatttaccgacctatcgacgaggagtgactcgaacggcgccggtaagccaggctactcagcaacaatgatagctaacatgcgcagggctcttggtgtccgtcatcgtgaaagcatatgatgccgtgatcagtgcgagcattggggctgagagacgtatgacgatccccatttccaagggcttcgacttcgggtattactgacggcgcatggaacgttagtgctatcgtctctacgcatatcaacgtattctggatggacgggtgccgggatactcgtcaagggttaagcctacaggcgcgctgagaatgctgtatataatcgcctgaactgtcgtgggctgccgcttggcacccccccgcaggctctgaattgcggggacttgagtccgcaagcggtgaggaggcgccgccactagccctggacgagatcccttatttggggagcacacaacctgttgttgaatccggtccgtcttgccctctgggggaccgcgatagccattctctgaataaaacgaaagattgggggatatagccttatcattgctgtctagtgatcgtgctttaacaaatggtagcgcttgcgctaactcggacaacaacaaccggtgagaaacgtctaaaccactcaatcaaccttaacctcccgttttatgcctttcacaccctgattacgttgacaaagagggacggcccaaggtcgactttcatactgaatcaggcaacctaaactagaaggaacccctgtaacccctcccccctcgtgagacgaggctgagtgctagaaattatgaaaccttgggtagttcactaaaacccaagtcatagccatactctaggggattactaccccggctgggatctaggaatgggagcagtatcgctaccaatattttccttcagtga"


#################

Leaves = 3
NumBrLengths = Leaves - 1

Tmax = 10 # start low for now
# Also need initial values of parameters.
ParameterSamples = matrix(rep(0,2*Tmax),ncol=2)
ParameterSamples[1,] = c(3,1) # Initial values, based upon my knowledge of coalescent theory. This is a sample for (T(3),T(2)).

# THIS FUNCTION GENERATES NEW PARAMETERS USING THE PROPOSAL DENSITY.
ProposalFunc = function(OldParams){
  # We'll use a independent bivariate normal as the proposal density as it is symmetric. Also need a cut-off so we don't get negative values.
  # Arbitrarily choose (hopefully sensible) variance parameters (1,1/3).
  x = 0 # Dummy
  while ( x = 0 ) {
    NewParams = rnorm(n = NumBrLengths, mean = OldParams, sd = c(1,1/3) )
    if ( min(NewParams > 0) ) { # need to check we haven't got a negative rate parameter
      return(NewParams)
      x = 1 # terminate
    } # end if
  } # end while
}

# THIS FUNCTION GIVES THE DENSITY OF THETA_1 GIVEN THETA_2
ProposalDens = function(ParamOfInterest,ParamGiven) {
  # I'm gonna try using an "apply" function!
  # How do you use sapply and dnorm, but with different mean and s.d.?
  Density = lapply(ParamOfInterest, function(i) dnorm(i,) )
}