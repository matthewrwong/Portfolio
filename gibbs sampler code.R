attach(faithful)
hist(waiting)


##p
mcmc_p_1step=function(p,m1,m2,s1,s2,data,delta){
  new_p=p+runif(1,-delta,delta)
  
  numerator=rep(0,length(data))
  denominator=rep(0,length(data))
  
  for(i in 1:length(data)){
    numerator[i]=new_p*dnorm(data[i],m1,sqrt(s1))+(1-new_p)*dnorm(data[i],m2,sqrt(s2))
    denominator[i]=p*dnorm(data[i],m1,sqrt(s1))+(1-p)*dnorm(data[i],m2,sqrt(s2))
  }
  
  n=sum(log(numerator))
  d=sum(log(denominator))
  a=exp(n-d)
  
  u3=runif(1)
  
  if(u3<a){
    return(new_p)
  }
  else{
    return(p)
  }
}

mcmc_p=function(p,m1,m2,s1,s2,data,delta,n){
  store=c(p,rep(0,n-1))
  for(i in 1:(n-1)){
    store[i+1]=mcmc_p_1step(store[i],m1,m2,s1,s2,data,delta)
  }
  return(store)
}


##mu1

mcmc_mu1_1step=function(p,m1,m2,s1,s2,data,delta){
  new_m1=m1+runif(1,-delta,delta)
  
  numerator=rep(0,length(data))
  denominator=rep(0,length(data))
  
  for(i in 1:length(data)){
    numerator[i]=p*dnorm(data[i],new_m1,sqrt(s1))+(1-p)*dnorm(data[i],m2,sqrt(s2))
    denominator[i]=p*dnorm(data[i],m1,sqrt(s1))+(1-p)*dnorm(data[i],m2,sqrt(s2))
  }
  
  n=sum(log(numerator))+log(dcauchy(new_m1))
  d=sum(log(denominator))+log(dcauchy(m1))
  a=exp(n-d)
  u3=runif(1)
  
  if(u3<a){
    return(new_m1)
  }
  else{
    return(m1)
  }
}

mcmc_mu1=function(p,m1,m2,s1,s2,data,delta,n){
  store=c(m1,rep(0,n-1))
  for(i in 1:(n-1)){
    store[i+1]=mcmc_mu1_1step(p,store[i],m2,s1,s2,data,delta)
  }
  return(store)
}

##mu2
mcmc_mu2_1step=function(p,m1,m2,s1,s2,data,delta){
  new_m2=m2+runif(1,-delta,delta)
  
  numerator=rep(0,length(data))
  denominator=rep(0,length(data))
  
  for(i in 1:length(data)){
    numerator[i]=p*dnorm(data[i],m1,sqrt(s1))+(1-p)*dnorm(data[i],new_m2,sqrt(s2))
    denominator[i]=p*dnorm(data[i],m1,sqrt(s1))+(1-p)*dnorm(data[i],m2,sqrt(s2))
  }
  
  n=sum(log(numerator))+log(dcauchy(new_m2))
  d=sum(log(denominator))+log(dcauchy(m2))
  a=exp(n-d)
  u3=runif(1)
  
  if(u3<a){
    return(new_m2)
  }
  else{
    return(m2)
  }
}

mcmc_mu2=function(p,m1,m2,s1,s2,data,delta,n){
  store=c(m2,rep(0,n-1))
  for(i in 1:(n-1)){
    store[i+1]=mcmc_mu2_1step(p,m1,store[i],s1,s2,data,delta)
  }
  return(store)
}

##sigma1

mcmc_sigma1_1step=function(p,m1,m2,s1,s2,data,delta){
  new_s1=s1+runif(1,-delta,delta)
  
  numerator=rep(0,length(data))
  denominator=rep(0,length(data))
  
  for(i in 1:length(data)){
    numerator[i]=p*dnorm(data[i],m1,sqrt(new_s1))+(1-p)*dnorm(data[i],m2,sqrt(s2))
    denominator[i]=p*dnorm(data[i],m1,sqrt(s1))+(1-p)*dnorm(data[i],m2,sqrt(s2))
  }
  
  n=sum(log(numerator))+log(dgamma(new_s1,.1,100))
  d=sum(log(denominator))+log(dgamma(s1,.1,100))
  a=exp(n-d)
  u3=runif(1)
  
  if(u3<a){
    return(new_s1)
  }
  else{
    return(s1)
  }
}

mcmc_sigma1=function(p,m1,m2,s1,s2,data,delta,n){
  store=c(s1,rep(0,n-1))
  for(i in 1:(n-1)){
    store[i+1]=mcmc_sigma1_1step(p,m1,m2,store[i],s2,data,delta)
  }
  return(store)
}


##sigma2

mcmc_sigma2_1step=function(p,m1,m2,s1,s2,data,delta){
  new_s2=s2+runif(1,-delta,delta)
  
  numerator=rep(0,length(data))
  denominator=rep(0,length(data))
  
  for(i in 1:length(data)){
    numerator[i]=p*dnorm(data[i],m1,sqrt(s1))+(1-p)*dnorm(data[i],m2,sqrt(new_s2))
    denominator[i]=p*dnorm(data[i],m1,sqrt(s1))+(1-p)*dnorm(data[i],m2,sqrt(s2))
  }
  
  n=sum(log(numerator))+log(dgamma(new_s2,.1,100))
  d=sum(log(denominator))+log(dgamma(s2,.1,100))
  a=exp(n-d)
  u3=runif(1)
  
  if(u3<a){
    return(new_s2)
  }
  else{
    return(s2)
  }
}

mcmc_sigma2=function(p,m1,m2,s1,s2,data,delta,n){
  store=c(s2,rep(0,n-1))
  for(i in 1:(n-1)){
    store[i+1]=mcmc_sigma2_1step(p,m1,m2,s1,store[i],data,delta)
  }
  return(store)
}


##gibbs sampler

gibbs_1step=function(p,m1,m2,s1,s2,data,deltap,deltamu,deltasigma){
  new_p=mcmc_p(p,m1,m2,s1,s2,data,deltap,101)[101]
  new_m1=mcmc_mu1(new_p,m1,m2,s1,s2,data,deltamu,101)[101]
  new_m2=mcmc_mu2(new_p,new_m1,m2,s1,s2,data,deltamu,101)[101]
  new_s1=mcmc_sigma1(new_p,new_m1,new_m2,s1,s2,data,deltasigma,101)[101]
  new_s2=mcmc_sigma1(new_p,new_m1,new_m2,new_s1,s2,data,deltasigma,101)[101]
  
  return(c(new_p,new_m1,new_m2,new_s1,new_s2))
}

gibbs_sampler=function(initial_p,initial_m1,initial_m2,initial_s1,initial_s2,data,deltap,deltamu,deltasigma,n){
  p=initial_p
  m1=initial_m1
  m2=initial_m2
  s1=initial_s1
  s2=initial_s2
  
  
  
  store=list()
  for(i in 1:n){
    store[[i]]=gibbs_1step(p,m1,m2,s1,s2,data,deltap,deltamu,deltasigma)
    p=store[[i]][1]
    m1=store[[i]][2]
    m2=store[[i]][3]
    s1=store[[i]][4]
    s2=store[[i]][5]
    print(i)
  }
  return(store)
}



test=gibbs_sampler(.5,50,80,5,7,waiting,.01,1,.1,1100)

samples_problem2=list()
for(i in 1:1000){
  samples_problem2[[i]]=test[[i+99]]
}

p_samples=rep(0,1000)
for(i in 1:1000){
  p_samples[i]=samples_problem2[[i]][1]
}

m1_samples=rep(0,1000)
for(i in 1:1000){
  m1_samples[i]=samples_problem2[[i]][2]
}

m2_samples=rep(0,1000)
for(i in 1:1000){
  m2_samples[i]=samples_problem2[[i]][3]
}

s1_samples=rep(0,1000)
for(i in 1:1000){
  s1_samples[i]=samples_problem2[[i]][4]
}

s2_samples=rep(0,1000)
for(i in 1:1000){
  s2_samples[i]=samples_problem2[[i]][5]
}



hist(p_samples,main="Posterior of p",xlab="p")
hist(m1_samples,main="Posterior of mu1",xlab="mu1")
hist(m2_samples,main="Posterior of mu2",xlab="mu2")
hist(s1_samples,main="Posterior of sigma1",xlab="sigma1")
hist(s2_samples,main="Posterior of sigma2",xlab="sigma2")



waiting_times=rep(0,300)

for(i in 1:300){
  u=runif(1)
  
  if(u<mean(p_samples)){
    waiting_times[i]=rnorm(1,mean(m1_samples),sqrt(s1_samples))
  }
  else{
    waiting_times[i]=rnorm(1,mean(m2_samples),sqrt(s2_samples))
  }
  
}


hist(waiting_times,main="Histogram of generated waiting times",xlab="Generated waiting times")
hist(waiting,main="Histogram of observed waiting times",xlab="Observed waiting times")