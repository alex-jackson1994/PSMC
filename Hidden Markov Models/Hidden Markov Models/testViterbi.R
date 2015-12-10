#####SCRIPTS#####
##Casino Example
A<-matrix(nrow=2,ncol=2)
A[,1]<-c(0.9,0.1)
A[,2]<-c(0.1,0.9)
B<-matrix(nrow=2,ncol=6)
B[1,]<-rep(1/6,6)
B[2,]<-c(rep(1/10,5),1/2)
obs<-c(1,2,3,4,5,6,6,6,6)
viterbi(A,B,c(1/2,1/2),obs)

##Coin Example
A<-matrix(nrow=2,ncol=2)
A[,1]<-c(0.9,0.1)
A[,2]<-c(0.1,0.9)
B<-matrix(nrow=2,ncol=2)
B[1,]<-c(1/2,1/2)
B[2,]<-c(0.1,0.9)
obs<-c(2,1,1,1,2,2,2,2,2,2)
viterbi(A,B,c(1/2,1/2),obs)
