#JJ: find largest possible J and number of parameters ntot
#ntotmax is upper bound on number of basis functions (parameters of bernstein polynomial); ntotmax=n is good rule of thumb
#p=number of dimension (2 in our case)  
JJ<-function(p,ntotmax){    
  Jmax<-p+2  #J order of the Bernstein polynomials, initialize at p+2 since it has to be greater than p
  ntot<-choose(Jmax-1,p-1) #ntot number of basis function (depends on J)
  #find the larger J such that the number of basis function is less than the upper bound of basis functions
  while(ntot<=ntotmax){
    Jmax<-Jmax+1
    ntot<-choose(Jmax-1,p-1)
  }
  Jmax<-Jmax-1 #J
  ntot<-choose(Jmax-1,p-1)#Number of basis function 
  JJ.output<-list("J"=Jmax,"ntot"=ntot)
  return(JJ.output)
}



#Ber.ind: make ntot by p array of Bernstein polynomial indices (alpha in the paper)
#indices are vector of length p whose component sum is equal to 5
Ber.ind<-function(Jmax,p,ntot){
  r<-rep(1,p)
  r[1]<-Jmax-p+1  #r firt index; vector of all ones except first element equal to J-p+1 
  t<-Jmax-p+1 # t and h then used to create indices
  h<-0
  ind<-array(0,dim=c(ntot,p)) # initialize array of Bernstein polynomial indices
  ind[1,1:p]<-r[1:p] #put the first indices in the first row of the array which contains indices
  for(i in 2:ntot){ #for each row of the array of ber. pol. indices
    if(t!=2){
      h<-0
    }
    h<-h+1
    t<-r[h]
    r[h]<-1                 #create index of length p whose component sum equals J
    r[1]<-t-1
    r[h+1]<-r[h+1]+1
    ind[i,1:p]<-r[1:p]
  }
  return(ind)
}





#Ver: put vertex coefficients into first p slots
#vertex coefficient = J-dimensional vector of all ones except element j is J 

Ver<-function(Jmax,p,ind){
  i<-1
  r<-rep(1,p)
  for(j in 2:p){
    i <- i + choose(Jmax-p-2+j,Jmax-p-1) #position of the indices with Jmax - p +1 in position j (vertex coefficients in the paper)
    r[1:p]<-ind[j,1:p]  #store the index in position j to change its position with the vertex coefficient in position i
    ind[j,1:p]<-ind[i,1:p]  #exchange the position of the two indices
    ind[i,1:p]<-r[1:p]
  } 
  return(ind)
}



#dir: compute Dirichlet densitites
#alpha corresponds to the indices, xin to the data
dir<-function(p,alpha,xin){
  #dirichlet density formula in ch.2 of the paper
  dir<- lgamma(sum(alpha))
  for(i in 1:p){
    dir<- dir- lgamma(alpha[i])+(alpha[i]-1)*log(xin[i])    
  }
  dir<-exp(dir)
}



#ll.cp: log-likelihood of distribution with change-point; equals zero if constraints not satisfied
#before computing the log-likelihood, it updates the parameters w by using the v
#n1, n2: sample sizes of the samples divided by the change-point
#dm1,dm2=density matrices (computed with function dir)
#w1, w2, v1, v2: parameters of Bernstein polynomial as in the paper
#c= precision parameter (0.1 recommended)



ll.cp<-function(c,p,Jmax,ntot,w1,w2,v1,v2,dm1,dm2,n1,n2){
  s1<-p+sum(exp(v1[1:(ntot-p)]))   #denominator of formula (8) in the paper
  for( i in (p+1):ntot){   #for each B.P. parameter in positions greater than p
    w1[i]<-exp(v1[i-p])/s1  #compute the parameters w depending on v; (8) in the paper
  }                          
  ll1<-0 #initialize ll for first group
  
  #same for the second group
  s2<-p+sum(exp(v2[1:(ntot-p)]))   
  for( i in (p+1):ntot){
    w2[i]<-exp(v2[i-p])/s2
  }                       
 
  ll2<-0
  
  flag<-0 #used to assure that parameters are between 0 and 1
  w1[1:p]<-1/p  #initialize first p parameters of the two groups
  w2[1:p]<-1/p
  for(i in (p+1):ntot){ #for each parameter in position greater than p
    for(j in 1:p){  #for each of the first p parameters
      w1[j]<-w1[j]-w1[i]*(ind[i,j]-1)/(Jmax-p)  #compute the first p parameters as in formula (7) in the paper
    }
  }
  #same for second group
  for(i in (p+1):ntot){
    for(j in 1:p){
      w2[j]<-w2[j]-w2[i]*(ind[i,j]-1)/(Jmax-p)    
    }
  }
  
  #check that parameters are between 0 and 1
  for(i in 1:p){   
    if(w1[i]<0|w1[i]>1){    
      flag<-1
    }
    if(w2[i]<0|w2[i]>1){
      flag<-1
    }
  }
  if(flag==0){ #compute log-posterior of the model as in the paper ch.3
    ll1<-sum(c*log(w1[1:ntot])) 
    ll1<-ll1+sum(log(dm1[1:n1,1:ntot]%*%w1[1:ntot]))
    ll2<-sum(c*log(w2[1:ntot])) 
    ll2<-ll2+sum(log(dm2[1:n2,1:ntot]%*%w2[1:ntot])) 
    ll<-ll1+ll2
  }else(ll<-0)
  ll.output<-list("w1"=w1,"ll"=ll,"w2"=w2)
  return(ll.output)
}



#MCMC.cp.det: MCMC algorithm
#GI is total number of Gibbs iterates
#tau:vector of length GI with the starting value as first entry
#llold: log likelihood of the distribution with the starting value as change-point (computed with ll.cp)


MCMC.cp.det<-function(GI,ntot,p,Jmax,c,tau,llold,w1,w2,v1,v2){
  chainw1<-array(dim=c(GI,ntot)) #initialize tha markov chains of the two distributions
  chainw2<-array(dim=c(GI,ntot))
  l<-0 #initialize the index used for the R&R variance method
  acc1<-rep(0,ntot-p)     #to compute acceptance ratio every 50 iterations for each parameter of the first distribution
  acc2<-rep(0,ntot-p)     #to compute acceptance ratio every 50 iterations for each parameter of the second distribution
  acct<-0                 #to compute acceptance ratio every 50 iterations for the change-point
  var1<-array(100,dim=c(GI/50+1,ntot-p))   #variance for every parameter of distribution 1
  var2<-array(100,dim=c(GI/50+1,ntot-p))   #variance for every parameter of distribution 2
  vart<-rep(2,GI/50+1)                     #variance of tau
  llold.tau<-0        #initialize log-likelihood for the distribution with the change-point estimated at the current iteration   
  llnew.tau<-0        #initialize log-likelihood for the distribution with the new proposed change-point    
  
  for(gibbs in 1:GI){ #GI iterations of the gibbs sampler
    
    
    taulik<- round(rtruncnorm(1,a= 2, b= n-2 , mean=tau[gibbs], sd=vart[l+1])) #propose a new tau
    
    x1.prop<-x[1:taulik,]        #divide data depending on the proposed tau to compute log likelihood
    x2.prop<-x[(taulik+1):n,]
    n1.prop<-nrow(x1.prop)       #store n for each group
    n2.prop<-nrow(x2.prop)
    dm1.prop<-array(0,dim=c(n1.prop,ntot))    #initialize density matrix (nrow= n of pseudo angles, ncol= n of parameters of BP) used to compute ll
    dm2.prop<-array(0,dim=c(n2.prop,ntot))
    for(i in 1:n1.prop){ #for each pseudo-angle
      for(j in 1:ntot){  #for eac parameter of the B.P.
        dm1.prop[i,j]<-dir(p,ind[j,1:p],x1.prop[i,1:p])       #store dirichlet densities in the density matrices (use function dir)
        
      }
    }
    
    #same for second group
    for(i in 1:n2.prop){
      for(j in 1:ntot){
        dm2.prop[i,j]<-dir(p,ind[j,1:p],x2.prop[i,1:p])
      }
    }
    
    #compute log likelihood of the distribution with changepoint estimated at the current iteration
    ll.output<-ll.cp(c,p,Jmax,ntot,w1,w2,v1,v2,dm1,dm2,n1,n2)
    llold.tau<-ll.output$ll
    
    #compute log likelihood of the distribution with the proposed change-point
    ll.output<-ll.cp(c,p,Jmax,ntot,w1,w2,v1,v2,dm1.prop,dm2.prop,n1.prop,n2.prop)
    llnew.tau<-ll.output$ll
    
    u<-runif(1) #random uniform for the acceptance ratio
    
    #acceptance ratio
    if(log(u)<(llnew.tau -llold.tau + 
               log(dtruncnorm(tau[gibbs],a=2 ,b= n-2, mean=taulik, sd=vart[l+1]))-
               log(dtruncnorm(taulik,a=2 ,b= n-2, mean=tau[gibbs], sd=vart[l+1])))
    )
     #if new tau accepted 
    { tau[gibbs+1]<-taulik  #store the proposed tau as estimated tau for the next iteration
    x1<-x1.prop     #store the properties of divided data with the proposed tau
    x2<-x2.prop
    n1<-n1.prop
    n2<-n2.prop
    dm1<-dm1.prop
    dm2<-dm2.prop
    acct<-acct+1  #update the index used for R&r variance 
    
    #if proposed tau not accepted
    }else{
      tau[gibbs+1]<-tau[gibbs]  #store the current tau as estimated changepoint for next iteration
    }
    
    
    #bernstein polynomial parameters
    
    for(k in 1:(ntot-p)){ #for each parameter v
      vold1<-v1[k] #store the current v 
      
      v1[k]<-v1[k]+rnorm(1)*sqrt(var1[l+1,k])     #Propose new v
      ll.output<-ll.cp(c,p,Jmax,ntot,w1,w2,v1,v2,dm1,dm2,n1,n2)    #log-likelihood of the model with the proposed v
      llnew<-ll.output$ll
      
      
      
      if(llnew!=0){ 
        u<-runif(1) #random uniform for acceptance ratio
        if(log(u)<(llnew-llold)){ #acceptance ratio
          #if new v accepted
          w1<-ll.output$w1   #store the parameters w estimated with the new v as estimated parameters for the next iteration
          llold<-llnew      #store the likelihood of the new model as likelihood of the estimated model for the next iteration
          acc1[k]<-acc1[k]+1  #update number of accepted proposals for R&R variance
          
          #if v not accepted
        }else(v1[k]<-vold1)  #put as v the old parameter stored before
        
      }else(v1[k]<-vold1)
    }
    
    #same for second group
    for(k in 1:(ntot-p)){
      vold2<-v2[k]
      v2[k]<-v2[k]+rnorm(1)*sqrt(var2[l+1,k])     #Propose v??* ~ N(v??,a??)
      ll.output<-ll.cp(c,p,Jmax,ntot,w1,w2,v1,v2,dm1,dm2,n1,n2)    #log-likelihood
      llnew<-ll.output$ll
      
      
      if(llnew!=0){
        u<-runif(1)
        if(log(u)<(llnew-llold)){ 
          w2<-ll.output$w2#probability of acceptance
          llold<-llnew
          acc2[k]<-acc2[k]+1
        }else(v2[k]<-vold2)
        
      }else(v2[k]<-vold2)
    }
    
    #R&R method of variance
    if(gibbs %% 50 == 0){ #every 50 iteration
      vl<-min(0.01,(l+1)^(-0.5))  #compute the weight used to update the variance
      for(k in 1:(ntot-p)){    #for each parameter v of the two distributions
        #if the proposed paramter has been accepted less than 22 times (22/50) decreases the variance; if the opposite happened, increases it.
        if(acc1[k]/50 < 0.44){var1[l+2,k] <- exp(log(var1[l+1,k]) - vl)} else {var1[l+2,k] <- exp(log(var1[l+1,k]) + vl)}
        if(acc2[k]/50 < 0.44){var2[l+2,k] <- exp(log(var2[l+1,k]) - vl)} else {var2[l+2,k] <- exp(log(var2[l+1,k]) + vl)}
      }
      #update variance of tau with the same method
      if(acct/50<0.44){vart[l+2]<-exp(log(vart[l+1]) - vl)} else {vart[l+2]<- exp(log(vart[l+1]) + vl)}
      
      l<-l+1 #update the parameter used to update the weight vl
      #reset the acceptance statistics
      acc1<-rep(0,ntot-p) 
      acc2<-rep(0,ntot-p)
      acct<-0
    }
    chainw1[gibbs,]<-w1   #store the estimated parameters of the B.P. in the markov chain
    chainw2[gibbs,]<-w2
    
    cat('iteration: ',gibbs)
    cat('\n')
    cat('estimated changepoint: ',tau[gibbs+1])
    cat('\n')
    cat('\n')
  }
  MCMC.output<-list("vart"=vart,"tau"=tau,"var1"=var1,"var2"=var2,"chainw1"=chainw1,"chainw2"=chainw2)
  return(MCMC.output)
}
#vart=chain of the variance of tau
#tau=chain of the change-point at each iteration
#var1=chain of the variance of the Bernstein polynomials parameter for first distribution
#var2=chain of the variance of the Bernstein polynomials parameter for second distribution
#chainw1 paramaters of Bernstein polynomial estimated at each iteration for first distribution
#chainw2 paramaters of Bernstein polynomial estimated at each iteration for second distribution


#log-likelihood and MCMC with just the bernstein polynomial based model (without change-point detection)


ll<-function(c,p,Jmax,ntot,w,v,dm,n){
  s<-p+sum(exp(v[1:(ntot-p)]))   #denominator of formula (8) in the paper
  for( i in (p+1):ntot){ #for each B.P. parameter in positions greater than p
    w[i]<-exp(v[i-p])/s     #compute the parameters w depending on v; (8) in the paper
  }                        
  ll<-0 #initialize log likelihood
  flag<-0  #used to assure that parameters are between 0 and 1
  w[1:p]<-1/p   #initialize first p parameters
  for(i in (p+1):ntot){   #for each parameter in position greater than p
    for(j in 1:p){    #for each of the first p parameter
      w[j]<-w[j]-w[i]*(ind[i,j]-1)/(Jmax-p)  #compute the first p parameters as in formula (7) in the paper
    }
  }
  #check that parameters are between 0 and 1
  for(i in 1:p){
    if(w[i]<0|w[i]>1){  
      flag<-1
    }
  }
  if(flag==0){ #if parameters are between 0 and 1
    #compute log-posterior of the model as in the paper ch.3
    ll<-sum(c*log(w[1:ntot]))  
    ll<-ll+sum(log(dm[1:n,1:ntot]%*%w[1:ntot])) 
    #if not log likelihood equals 0
  }else(ll<-0)
  ll.output<-list("w"=w,"ll"=ll)
  return(ll.output)
}


#GI number of iterations of the gibbs sampler
#ntot nummber of B.P. parameters
#p nummber of dimension (2 in our case)
#Jmax needed to compute ll (order of B.P.)
#v and w parameters of B.P. as in paper
#dm density matrix which contains dirichlet densities
#n number of data
MCMC<-function(GI,ntot,p,Jmax,v,w,dm,n){
  #initialize the markov chains for v and w
  chainw<-array(dim=c(GI,ntot))
  chainv<-array(dim=c(GI,ntot-p))
  l<-0 #initialize parameter needed to update variance with R&R method
  acc<-rep(0,ntot-p) #needed to store number of accepted proposals for each parameter
  var<-array(100,dim=c(GI/50+1,ntot-p)) #variance for each of the parameter initialised at 100 
  
  for(gibbs in 1:GI){ #for each iteration of the M-H
    
    for(k in 1:(ntot-p)){  #for each parameter v
      
      vold<-v[k] #store the current parameter estimated
      v[k]<-v[k]+rnorm(1)*sqrt(var[l+1,k])     #propose a new parameter v
      ll.output<-ll(c,p,Jmax,ntot,w,v,dm,n)    #compute likelihood of the distribution with the proposed parameter
      llnew<-ll.output$ll
      if(llnew!=0){  
        u<-runif(1) #random uniform for acceptance ratio
        if(log(u)<(llnew-llold)){ #acceptance ratio
      #if parameter accepted
          w<-ll.output$w  #store the new parameters w computed with the updated v
          llold<-llnew   #store the log likelihood of the distribution with the accepted parameter as likelihood for the next iteration
          acc[k]<-acc[k]+1 #update number of accepted parameters for R&R variance
          #if not accepted store the old v as current v
        }else(v[k]<-vold)
        
      }else(v[k]<-vold)
    }
    #every 50 iteration
    if(gibbs %% 50 == 0){
      vl<-min(0.01,(l+1)^(-0.5)) #update the weight used to update variances
      for(k in 1:(ntot-p)){ #for each parameter v
        #if proposal accepted less than 22 times in the last 50 iterations for a specific parameter, decrease its variance; if not increases it
        if(acc[k]/50 < 0.44){var[l+2,k] <- exp(log(var[l+1,k]) - vl)} else {var[l+2,k] <- exp(log(var[l+1,k]) + vl)}  
      }
      l<-l+1 #update parameter used to compute the weight vl
      acc<-rep(0,ntot-p) #reset the acceptance vector
    }
    #store the estimated parameters v and w in the markov chains
    chainw[gibbs,]<-w
    chainv[gibbs,]<-v
    print(gibbs)
    
  }
  MCMC.output<-list("chainw"=chainw,"chainv"=chainv,"var"=var)
  return(MCMC.output)
}


