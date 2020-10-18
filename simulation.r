library(evd)
library(truncnorm)
library(MCMCpack)
source('functions.r')

#simulate random bivariate observations from logistic dist.
b1<-rbvevd(400,dep=0.2,model = "log")
b2<-rbvevd(600,dep=0.8,model = "log")
b<-rbind(b1,b2)

#compute pseudo-angles for the two simulated distributions W1 and W2 after having transformed the obs. in Pareto unit dist.

ecdf.b11<-ecdf(b1[,1])
X1<-1/(1-ecdf.b11(b1[,1]-0.0000001))
ecdf.b12<-ecdf(b1[,2])
Y1<-1/(1-ecdf.b12(b1[,2]-0.0000001))

R1<-X1+Y1
WW1<-X1/(R1)

#define extreme values the obs. that exceed the threshold of 90%
threshold1<-as.numeric(quantile(R1,0.90)) 

exceedences1<-which(R1>threshold1)

W1<-WW1[exceedences1]


ecdf.b21<-ecdf(b2[,1])
X2<-1/(1-ecdf.b21(b2[,1]-0.0000001))
ecdf.b22<-ecdf(b2[,2])
Y2<-1/(1-ecdf.b22(b2[,2]-0.0000001))

R2<-X2+Y2
WW2<-X2/(R2)

threshold2<-as.numeric(quantile(R2,0.90))

exceedences2<-which(R2>threshold2)

W2<-WW2[exceedences2]

#combine them as they were from the same dist.
W<-c(W1,W2)


#create the matrix that will be analyzed by the bernstein pol.-based model
x<- cbind(W,1-W) 

n<-nrow(x) 

#starting value of the change-point
tau<-c(as.integer(n/2)) 

#create the two samples divided by the s.v. of the change-point
x1<-x[1:tau[1],]
x2<-x[(tau[1]+1):n,]

n1<-nrow(x1)
n2<-nrow(x2)

p<-ncol(x)

#iterations of MCMC
GI<-1000  

ntotmax<-as.integer(n/2) #good rule of thumb is n 

c<-0.1

#find largest J s.t. ntot <= ntotmax
JJ.output<-JJ(p,ntotmax)
Jmax<-JJ.output$J
ntot<-JJ.output$ntot

#make ntot by p array of Bernstein polynomial indices
ind<-Ber.ind(Jmax,p,ntot)

#put vertex coefficients into first p slots
ind<-Ver(Jmax,p,ind)

#create the density matrix by using function dir: compute Dirichlet densitites
dm1<-array(0,dim=c(n1,ntotmax))
dm2<-array(0,dim=c(n2,ntotmax))
for(i in 1:n1){ #for each component of the first group
  for(j in 1:ntot){ #for each B.P. index
    dm1[i,j]<-dir(p,ind[j,1:p],x1[i,1:p]) #store the dirichlet density of the i-th obs with the j-th index in the density matrix
  }
}

#same for second group
for(i in 1:n2){
  for(j in 1:ntot){
    dm2[i,j]<-dir(p,ind[j,1:p],x2[i,1:p])
  }
}

#initialize parameters for ber.pol.

w1<-c()
w2<-c()
for(i in 1:ntot){
  w1[i]<-1/ntot
}
v1<- rep(0,ntot-p)


for(i in 1:ntot){
  w2[i]<-1/ntot
}
v2<- rep(0,ntot-p)

#initialize the log-likelihood considering the s.v. of tau as change-point
ll.output<-ll.cp(c,p,Jmax,ntot,w1,w2,v1,v2,dm1,dm2,n1,n2)
llold<-ll.output$ll
w1<-ll.output$w1
w2<-ll.output$w2


#run change-point Bernstein polynomials-based model
MCMC.output<-MCMC.cp.det(GI,ntot,p,Jmax,c,tau,llold,w1,w2,v1,v2)

chainw1<-MCMC.output$chainw1
chainw2<-MCMC.output$chainw2
var1<-MCMC.output$var1
var2<-MCMC.output$var2
vart<-MCMC.output$vart
tau<-MCMC.output$tau

par(mfrow=c(1,2))
plot(tau,type="l")
hist(tau,col="cornflowerblue",breaks=10)

#compute the mode of the chain of tau to obtain the estimated change-point
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
changepoint<-Mode(tau)


s<-seq(401,1000,20)  #thinning


#compute the estimated densities of the two distributions divided by the estimated changepoint by
#computing the mean of the estimated densities by the paramters chosen with thinning

d<-seq(0.01,0.99,length.out = 1000)
d<-cbind(d,rep(1,1000)-d)
dmm<-array(dim=c(1000,ntot))
for(i in 1:1000){
  for(j in 1:ntot){
    dmm[i,j]<-dir(p,ind[j,1:p],d[i,1:p])
  }
}

estimated.density11<-array(dim=c(length(s),1000))
estimated.density1<-c()

estimated.density22<-array(dim=c(length(s),1000))
estimated.density2<-c()

for(i in 1:length((s))){
  estimated.density11[i,1:1000]<-dmm[1:1000,1:ntot]%*%chainw1[s[i],1:ntot]
}


for(i in 1:length((s))){
  estimated.density22[i,1:1000]<-dmm[1:1000,1:ntot]%*%chainw2[s[i],1:ntot]
}

for(i in 1:1000){
  estimated.density1[i]<-mean(estimated.density11[,i])
}

for(i in 1:1000){
  estimated.density2[i]<-mean(estimated.density22[,i])
}


#plot the histograms of the pseudo-angles of the two distributions and the respective estimated densities
hist(x[1:changepoint,1],freq=F,breaks=10,xlab="W",col="lightblue")
lines(d[,1],estimated.density1,lwd=2,col="firebrick")
lines(seq(0,1,length.out=100),
      hbvevd(seq(0,1,length.out=100),dep=0.1,half=T),lwd=2,col="darkblue")
legend("topright",legend=c("estimated.density","real.density"),lwd=2,
       col=c("firebrick","darkblue"),cex=0.85)




hist(x[(changepoint+1):n,1],freq=F,breaks=10,xlab="W",col="lightblue")
lines(d[,1],estimated.density2,lwd=2,col="firebrick")
lines(seq(0,1,length.out=100),
      hbvevd(seq(0,1,length.out=100),dep=0.9,half=T),lwd=2,col="darkblue")
legend("topright",y=4,legend=c("estimated.density","real.density"),lwd=2,
       col=c("firebrick","darkblue"),cex=0.85)













