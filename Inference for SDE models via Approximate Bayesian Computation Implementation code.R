set.seed=1124234248174

##create data
beta=5
sigma2=20

t=seq(1,100,by=1)

x_true=exp(beta-.5)*t+rnorm(100,0,1)
y_true=x_true+rnorm(100,0,sqrt(sigma2))

##MCMC algorithm

MCMC=function(beta_initial,sigma2_initial,d_initial,reps){
  x=exp(beta_initial-.5)*t+rnorm(100,0,1)
  y=x+rnorm(100,0,sqrt(sigma2_initial))
  beta=beta_initial
  sigma2=sigma2_initial
  d=d_initial
  
  beta_vector=rep(0,reps)
  sigma2_vector=rep(0,reps)
  d_vector=rep(0,reps)
  
  for(i in 1:reps){
    beta_new=beta+runif(1,-.3,.3)
    sigma2_new=sigma2+runif(1,-.3,.3)
    d_new=d+runif(1,-1,1)
    
    x_new=exp(beta_new-.5)*t+rnorm(100,0,1)
    y_new=x_new+rnorm(100,0,sqrt(sigma2_new))
    
    y_new_bar=mean(y_new)
    y_old_bar=mean(y)
    
    alpha_num=dnorm(abs(y_new_bar-mean(y_true))/d_new)/d_new
    alpha_den=dnorm(abs(y_old_bar-mean(y_true))/d)/d
    
    alpha=alpha_num/alpha_den
    
    u=runif(1)
    if(alpha>u){
      beta=beta_new
      sigma2=sigma2_new
      d=d_new
      x=x_new
      y=y_new
    }
    
    beta_vector[i]=beta
    sigma2_vector[i]=sigma2
    d_vector[i]=d
  }
  
  return(list(beta_vector,sigma2_vector,d_vector))
}



##
set.seed=1143414231
test=MCMC(1,10,200,10000)
test[[1]]=test[[1]][-seq(1,5000)]
test[[2]]=test[[2]][-seq(1,5000)]
test[[3]]=test[[3]][-seq(1,5000)]
hist(test[[1]],main="Posterior for beta",xlab="beta")
hist(test[[2]],main="Posterior for sigma2",xlab="sigma2")
hist(test[[3]],main="Posterior for delta",xlab="delta")