
library(truncnorm)
library(RobustRankAggreg)
library(TopKLists)
sigma2_rho=function(rho){rho^(-2)-1}
source(file="http://www.pitt.edu/~mchikina/BIRRA/runBIRRA.R")

#define function to sample from truncated gamma distribution
qtruncgamma=function(a,b,p,shape,rate=1)
{qgamma(pgamma(a,shape=shape,rate=rate)+p*(pgamma(b,shape=shape,rate=rate)-pgamma(a,shape=shape,rate=rate)),shape=shape,rate=rate)}
rtruncgamma=function(n,a,b,shape,rate=1)
{u = runif(n, min = 0, max = 1)
x = qtruncgamma(a=a,b=b,p=u,shape=shape,rate=rate)
return(x)}

#define function to implement diffuse gamma prior
BiG_gamma=function(r,G,S,M=20000,ds=1,dp=1,W=w,sigma_p10=0.5,sigma_p20=0.5,mu0=numeric(G),
                   kappa10=numeric(G),kappa20=numeric(G),sigma_s0=rep(1,S)){
  Kappa2=Kappa1=Mu=matrix(0,G,M)
  Sigma_s2=matrix(0,S,M)
  Sigma_p12=Sigma_p22=numeric(M)
  Mu[,1]=mu0
  Sigma_s2[,1]=sigma_s0
  Kappa1[,1]=kappa10
  Kappa2[,1]=kappa20
  Sigma_p12[1]=sigma_p10
  Sigma_p22[1]=sigma_p20
  for (m in 2:M)
  {
    for (g in 1:G)
    {
      ## update mu
      c=sum((W[g,is.na(w[g,])==F]-c(rep(Kappa1[g,m-1],S/2),rep(Kappa2[g,m-1],S/2))[is.na(w[g,])==F])/Sigma_s2[is.na(w[g,])==F,m-1])
      d=sum(1/Sigma_s2[is.na(w[g,])==F,m-1])+1
      Mu[g,m]=rnorm(1,c/d,sqrt(1/d))
      ## update kappa
      c=sum((W[g,c(is.na(w[g,1:(S/2)])==F,rep(F,S/2))]-Mu[g,m])/Sigma_s2[c(is.na(w[g,1:(S/2)])==F,rep(F,S/2)),m-1])
      d=sum(1/Sigma_s2[c(is.na(w[g,1:(S/2)])==F,rep(F,S/2)),m-1])+1/Sigma_p12[m-1]
      Kappa1[g,m]=rnorm(1,c/d,sqrt(1/d))
      c=sum((W[g,c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F)]-Mu[g,m])/Sigma_s2[c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F),m-1])
      d=sum(1/Sigma_s2[c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F),m-1])+1/Sigma_p22[m-1]
      Kappa2[g,m]=rnorm(1,c/d,sqrt(1/d))
    }
    ## updata sigma_s2
    for (s in 1:(S/2))
    {
      Sigma_s2[s,m]=1/rgamma(1,shape=sum(is.na(w[,s])==F)/2+ds,
                             rate=sum((W[,s]-Mu[,m]-Kappa1[,m])^2,na.rm=T)/2+ds)
    }
    for (s in (S/2+1):S)
    {
      Sigma_s2[s,m]=1/rgamma(1,shape=sum(is.na(w[,s])==F)/2+ds,
                             rate=sum((W[,s]-Mu[,m]-Kappa2[,m])^2,na.rm=T)/2+ds)
    }
    #update sigma_p2
    P1=as.numeric(names(apply(is.na(w[,1:(S/2)]),1,prod)==0))
    P2=as.numeric(names(apply(is.na(w[,(S/2+1):S]),1,prod)==0))
    Sigma_p12[m]=1/rgamma(1,shape=length(P1)/2+dp,
                          rate=sum(Kappa1[P1,m]^2,na.rm=T)/2+dp)
    Sigma_p22[m]=1/rgamma(1,shape=length(P2)/2+dp,
                          rate=sum(Kappa2[P2,m]^2,na.rm=T)/2+dp)
    #update w
    for (l in 1:10)
    {
      for (s in 1:(S/2))
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {aa=max(W[which(r[,s]==(n_T[s]+1)),s]);bb=W[which(r[,s]==(n_T[s]-1)),s]} 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+Kappa1[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        W[which(r[,s]==(n_T[s]+1)),s]=rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=T),b=W[which(r[,s]==n_T[s]),s],
                                                 mean=Mu[which(r[,s]==n_T[s]+1),m]+Kappa1[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
      }
      for (s in (S/2+1):S)
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {aa=max(W[which(r[,s]==(n_T[s]+1)),s]);bb=W[which(r[,s]==(n_T[s]-1)),s]} 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+Kappa2[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        W[which(r[,s]==(n_T[s]+1)),s]=rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=T),b=W[which(r[,s]==n_T[s]),s],
                                                 mean=Mu[which(r[,s]==n_T[s]+1),m]+Kappa2[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
      }
    }
  }
  post.mean.mu=apply(Mu[,ceiling(M/2):M],1,mean)
  return(post.mean.mu)
}

#define function to implement diffuse uniform prior
BiG_unif=function(r,G,S,M=20000,a=0.0202,b=98.5025,W=w,sigma_p10=0.5,sigma_p20=0.5,mu0=numeric(G),
                  kappa10=numeric(G),kappa20=numeric(G),sigma_s0=rep(1,S)){
  Kappa2=Kappa1=Mu=matrix(0,G,M)
  Sigma_s2=matrix(0,S,M)
  Sigma_p12=Sigma_p22=numeric(M)
  Mu[,1]=mu0
  Sigma_s2[,1]=sigma_s0
  Kappa1[,1]=kappa10
  Kappa2[,1]=kappa20
  Sigma_p12[1]=sigma_p10
  Sigma_p22[1]=sigma_p20
  for (m in 2:M)
  {
    for (g in 1:G)
    {
      ## update mu
      c=sum((W[g,is.na(w[g,])==F]-c(rep(Kappa1[g,m-1],S/2),rep(Kappa2[g,m-1],S/2))[is.na(w[g,])==F])/Sigma_s2[is.na(w[g,])==F,m-1])
      d=sum(1/Sigma_s2[is.na(w[g,])==F,m-1])+1
      Mu[g,m]=rnorm(1,c/d,sqrt(1/d))
      ## update kappa
      c=sum((W[g,c(is.na(w[g,1:(S/2)])==F,rep(F,S/2))]-Mu[g,m])/Sigma_s2[c(is.na(w[g,1:(S/2)])==F,rep(F,S/2)),m-1])
      d=sum(1/Sigma_s2[c(is.na(w[g,1:(S/2)])==F,rep(F,S/2)),m-1])+1/Sigma_p12[m-1]
      Kappa1[g,m]=rnorm(1,c/d,sqrt(1/d))
      c=sum((W[g,c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F)]-Mu[g,m])/Sigma_s2[c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F),m-1])
      d=sum(1/Sigma_s2[c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F),m-1])+1/Sigma_p22[m-1]
      Kappa2[g,m]=rnorm(1,c/d,sqrt(1/d))
    }
    ## updata sigma_s2
    for (s in 1:(S/2))
    {
      Sigma_s2[s,m]=1/rtruncgamma(1,a=a,b=b,shape=(sum(is.na(w[,s])==F)-1)/2,
                                  rate=sum((W[,s]-Mu[,m]-Kappa1[,m])^2,na.rm=T)/2)
    }
    for (s in (S/2+1):S)
    {
      Sigma_s2[s,m]=1/rtruncgamma(1,a=a,b=b,shape=(sum(is.na(w[,s])==F)-1)/2,
                                  rate=sum((W[,s]-Mu[,m]-Kappa2[,m])^2,na.rm=T)/2)
    }
    #update sigma_p2
    P1=as.numeric(names(apply(is.na(w[,1:(S/2)]),1,prod)==0))
    P2=as.numeric(names(apply(is.na(w[,(S/2+1):S]),1,prod)==0))
    Sigma_p12[m]=1/rtruncgamma(1,a=a,b=b,shape=(length(P1)-1)/2,
                               rate=sum(Kappa1[P1,m]^2,na.rm=T)/2)
    Sigma_p22[m]=1/rtruncgamma(1,a=a,b=b,shape=(length(P2)-1)/2,
                               rate=sum(Kappa2[P2,m]^2,na.rm=T)/2)
    #update w
    for (l in 1:10)
    {
      for (s in 1:(S/2))
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {aa=max(W[which(r[,s]==(n_T[s]+1)),s]);bb=W[which(r[,s]==(n_T[s]-1)),s]} 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+Kappa1[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        W[which(r[,s]==(n_T[s]+1)),s]=rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=T),b=W[which(r[,s]==n_T[s]),s],
                                                 mean=Mu[which(r[,s]==n_T[s]+1),m]+Kappa1[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
      }
      for (s in (S/2+1):S)
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {aa=max(W[which(r[,s]==(n_T[s]+1)),s]);bb=W[which(r[,s]==(n_T[s]-1)),s]} 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+Kappa2[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        W[which(r[,s]==(n_T[s]+1)),s]=rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=T),b=W[which(r[,s]==n_T[s]),s],
                                                 mean=Mu[which(r[,s]==n_T[s]+1),m]+Kappa2[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
      }
    }
  }
  post.mean.mu=apply(Mu[,ceiling(M/2):M],1,mean)
  return(post.mean.mu)
}

#define function to implement FB
BiG_FB=function(r,G,S,M=10000,W=w,sigma_p10=0.5,sigma_p20=0.5,mu0=numeric(G),
                kappa10=numeric(G),kappa20=numeric(G),sigma_s0=rep(1,S)){
  Kappa2=Kappa1=Mu=matrix(0,G,M)
  Sigma_s2=matrix(0,S,M)
  Sigma_p12=Sigma_p22=numeric(M)
  mu.s=mu.p=nu.s=nu.p=numeric(M)
  Mu[,1]=mu0
  Sigma_s2[,1]=sigma_s0
  Kappa1[,1]=kappa10
  Kappa2[,1]=kappa20
  Sigma_p12[1]=sigma_p10
  Sigma_p22[1]=sigma_p20
  mu.s[1]=mu.p[1]=1
  nu.s[1]=nu.p[1]=1
  for (m in 2:M)
  {
    for (g in 1:G)
    {
      ## update mu
      c=sum((W[g,is.na(w[g,])==F]-c(rep(Kappa1[g,m-1],S/2),rep(Kappa2[g,m-1],S/2))[is.na(w[g,])==F])/Sigma_s2[is.na(w[g,])==F,m-1])
      d=sum(1/Sigma_s2[is.na(w[g,])==F,m-1])+1
      Mu[g,m]=rnorm(1,c/d,sqrt(1/d))
      ## update kappa
      c=sum((W[g,c(is.na(w[g,1:(S/2)])==F,rep(F,S/2))]-Mu[g,m])/Sigma_s2[c(is.na(w[g,1:(S/2)])==F,rep(F,S/2)),m-1])
      d=sum(1/Sigma_s2[c(is.na(w[g,1:(S/2)])==F,rep(F,S/2)),m-1])+1/Sigma_p12[m-1]
      Kappa1[g,m]=rnorm(1,c/d,sqrt(1/d))
      c=sum((W[g,c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F)]-Mu[g,m])/Sigma_s2[c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F),m-1])
      d=sum(1/Sigma_s2[c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F),m-1])+1/Sigma_p22[m-1]
      Kappa2[g,m]=rnorm(1,c/d,sqrt(1/d))
    }
    ## updata sigma_s2
    for (s in 1:(S/2))
    {
      Sigma_s2[s,m]=1/rgamma(1,shape=sum(is.na(w[,s])==F)/2+nu.s[m-1],
                             rate=sum((W[,s]-Mu[,m]-Kappa1[,m])^2,na.rm=T)/2+nu.s[m-1]/mu.s[m-1])
    }
    for (s in (S/2+1):S)
    {
      Sigma_s2[s,m]=1/rgamma(1,shape=sum(is.na(w[,s])==F)/2+nu.s[m-1],
                             rate=sum((W[,s]-Mu[,m]-Kappa2[,m])^2,na.rm=T)/2+nu.s[m-1]/mu.s[m-1])
    }
    #update sigma_p2
    P1=as.numeric(names(apply(is.na(w[,1:(S/2)]),1,prod)==0))
    P2=as.numeric(names(apply(is.na(w[,(S/2+1):S]),1,prod)==0))
    Sigma_p12[m]=1/rgamma(1,shape=length(P1)/2+nu.p[m-1],
                          rate=sum(Kappa1[P1,m]^2,na.rm=T)/2+nu.p[m-1]/mu.p[m-1])
    Sigma_p22[m]=1/rgamma(1,shape=length(P2)/2+nu.p[m-1],
                          rate=sum(Kappa2[P2,m]^2,na.rm=T)/2+nu.p[m-1]/mu.p[m-1])
    #update hyperparameters mu.p, mu.s, nu.p, nu.s
    h=0
    while (h==0)
    {
      mu_star=1/rgamma(1,shape=S*nu.s[m-1]-1,scale=nu.s[m-1]*sum(1/Sigma_s2[,m])) #J*nu-1>0
      p=min(1,exp(-0.05*(mu_star-mu.s[m-1])))
      mu.s[m]=mu_star
      h=(runif(1)<p)
    }
    h=0
    while (h==0)
    {
      mu_star=1/rgamma(1,shape=S*nu.p[m-1]-1,scale=nu.p[m-1]*(1/Sigma_p12[m]+1/Sigma_p22[m])) #J*nu-1>0
      p=min(1,exp(-0.05*(mu_star-mu.p[m-1])))
      mu.p[m]=mu_star
      h=(runif(1)<p)
    }
    h=0
    while(h==0)
    {
      nu_star=exp(log(nu.s[m-1])+rnorm(1,0,1))
      while(nu_star<(1/S))
      {nu_star=exp(log(nu.s[m-1])+rnorm(1,0,1))}
      p=min(1,((nu_star^-1.17)*exp(-0.65*(1/nu_star-1/nu.s[m-1]))/(nu.s[m-1]^-1.17))*
              exp(sum(nu_star*log(1/Sigma_s2[,m]*nu_star/mu.s[m])-(1/Sigma_s2[,m])*(nu_star-nu.s[m-1])/mu.s[m-1]+lgamma(nu.s[m-1])
                      -lgamma(nu_star)-log((1/Sigma_s2[,m])*nu.s[m-1]/mu.s[m])*nu.s[m-1])))
      nu.s[m]=nu_star
      h=(runif(1)<p)
    }
    h=0
    while(h==0)
    {
      nu_star=exp(log(nu.p[m-1])+rnorm(1,0,0.5))
      while(nu_star<(1/2))
      {nu_star=exp(log(nu.p[m-1])+rnorm(1,0,0.5))}
      p=min(1,((nu_star^-1.17)*exp(-0.65*(1/nu_star-1/nu.p[m-1]))/(nu.p[m-1]^-1.17))*
              exp(nu_star*log(1/Sigma_p12[m]*nu_star/mu.p[m])-(1/Sigma_p12[m])*(nu_star-nu.p[m-1])/mu.p[m-1]+lgamma(nu.p[m-1])
                  -lgamma(nu_star)-log((1/Sigma_p12[m])*nu.p[m-1]/mu.p[m])*nu.p[m-1]+
                    nu_star*log(1/Sigma_p22[m]*nu_star/mu.p[m])-(1/Sigma_p22[m])*(nu_star-nu.p[m-1])/mu.p[m-1]+lgamma(nu.p[m-1])
                  -lgamma(nu_star)-log((1/Sigma_p22[m])*nu.p[m-1]/mu.p[m])*nu.p[m-1]))
      nu.p[m]=nu_star
      h=(runif(1)<p)
    }
    #update w
    for (l in 1:10)
    {
      for (s in 1:(S/2))
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {aa=max(W[which(r[,s]==(n_T[s]+1)),s]);bb=W[which(r[,s]==(n_T[s]-1)),s]} 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+Kappa1[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        W[which(r[,s]==(n_T[s]+1)),s]=rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=T),b=W[which(r[,s]==n_T[s]),s],
                                                 mean=Mu[which(r[,s]==n_T[s]+1),m]+Kappa1[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
      }
      for (s in (S/2+1):S)
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {aa=max(W[which(r[,s]==(n_T[s]+1)),s]);bb=W[which(r[,s]==(n_T[s]-1)),s]} 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+Kappa2[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        W[which(r[,s]==(n_T[s]+1)),s]=rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=T),b=W[which(r[,s]==n_T[s]),s],
                                                 mean=Mu[which(r[,s]==n_T[s]+1),m]+Kappa2[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
      }
    }
  }
  post.mean.mu=apply(Mu[,ceiling(M/2):M],1,mean)
  return(post.mean.mu)
}

#define function to implement DA

BiG_DA=function(r,G,S,M=20000,a=0.0202,b=98.5025,dp=1,W=w,sigma_p10=0.5,sigma_p20=0.5,mu0=numeric(G),
                xi10=numeric(G),xi20=numeric(G),sigma_s0=rep(1,S)){
  xi2=xi1=Mu=matrix(0,G,M)
  Sigma_s2=matrix(0,S,M)
  alpha1=alpha2=Sigma_p12=Sigma_p22=numeric(M)
  Mu[,1]=mu0
  Sigma_s2[,1]=sigma_s0
  xi1[,1]=xi10
  xi2[,1]=xi20
  alpha1[1]=alpha2[1]=1
  Sigma_p12[1]=sigma_p10
  Sigma_p22[1]=sigma_p20
  for (m in 2:M)
  {
    for (g in 1:G)
    {
      ## update mu
      c=sum((W[g,is.na(w[g,])==F]-c(rep(alpha1[m-1]*xi1[g,m-1],S/2),rep(alpha2[m-1]*xi2[g,m-1],S/2))[is.na(w[g,])==F])/Sigma_s2[is.na(w[g,])==F,m-1])
      d=sum(1/Sigma_s2[is.na(w[g,])==F,m-1])+1
      Mu[g,m]=rnorm(1,c/d,sqrt(1/d))
      ## update xi
      c=alpha1[m-1]*sum((W[g,c(is.na(w[g,1:(S/2)])==F,rep(F,S/2))]-Mu[g,m])/Sigma_s2[c(is.na(w[g,1:(S/2)])==F,rep(F,S/2)),m-1])
      d=alpha1[m-1]^2*sum(1/Sigma_s2[c(is.na(w[g,1:(S/2)])==F,rep(F,S/2)),m-1])+1/Sigma_p12[m-1]
      xi1[g,m]=rnorm(1,c/d,sqrt(1/d))
      c=alpha2[m-1]*sum((W[g,c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F)]-Mu[g,m])/Sigma_s2[c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F),m-1])
      d=alpha2[m-1]^2*sum(1/Sigma_s2[c(rep(F,S/2),is.na(w[g,(S/2+1):S])==F),m-1])+1/Sigma_p22[m-1]
      xi2[g,m]=rnorm(1,c/d,sqrt(1/d))
    }
    ##update alpha
    P1=as.numeric(names(apply(is.na(w[,1:(S/2)]),1,prod)==0))
    P2=as.numeric(names(apply(is.na(w[,(S/2+1):S]),1,prod)==0))
    c=sum(t((W[P1,1:(S/2)]-Mu[P1,m])*xi1[P1,m])/Sigma_s2[1:(S/2),m-1],na.rm=T)
    d=sum(xi1[P1,m]^2)*sum(1/Sigma_s2[1:(S/2),m-1])+1
    alpha1[m]=rnorm(1,c/d,sqrt(1/d))
    c=sum(t((W[P2,(S/2+1):S]-Mu[g,m])*xi2[P2,m])/Sigma_s2[(S/2+1):S,m-1],na.rm=T)
    d=sum(xi2[P2,m]^2)*sum(1/Sigma_s2[(S/2+1):S,m-1])+1
    alpha2[m]=rnorm(1,c/d,sqrt(1/d))
    #update sigma_p2
    Sigma_p12[m]=1/rgamma(1,shape=length(P1)/2+1,
                          rate=sum(xi1[P1,m]^2,na.rm=T)/2+1)
    Sigma_p22[m]=1/rgamma(1,shape=length(P2)/2+1,
                          rate=sum(xi2[P2,m]^2,na.rm=T)/2+1)
    ##update sigma_s2
    for (s in 1:(S/2))
    {
      Sigma_s2[s,m]=1/rtruncgamma(1,a=a,b=b,shape=(sum(is.na(w[,s])==F)-1)/2,
                                  rate=sum((W[,s]-Mu[,m]-alpha1[m]*xi1[,m])^2,na.rm=T)/2)
    }
    for (s in (S/2+1):S)
    {
      Sigma_s2[s,m]=1/rtruncgamma(1,a=a,b=b,shape=(sum(is.na(w[,s])==F)-1)/2,
                                  rate=sum((W[,s]-Mu[,m]-alpha2[m]*xi2[,m])^2,na.rm=T)/2)
    }
    
    #update w
    for (l in 1:10)
    {
      for (s in 1:(S/2))
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {aa=max(W[which(r[,s]==(n_T[s]+1)),s]);bb=W[which(r[,s]==(n_T[s]-1)),s]} 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+alpha1[m]*xi1[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        W[which(r[,s]==(n_T[s]+1)),s]=rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=T),b=W[which(r[,s]==n_T[s]),s],
                                                 mean=Mu[which(r[,s]==n_T[s]+1),m]+alpha1[m]*xi1[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
      }
      for (s in (S/2+1):S)
      {
        index=sample(1:n_T[s])
        for (g in index)
        {
          if (g==1) {aa=W[which(r[,s]==2),s];bb=Inf} 
          else if (g==n_T[s]) {aa=max(W[which(r[,s]==(n_T[s]+1)),s]);bb=W[which(r[,s]==(n_T[s]-1)),s]} 
          else {aa=W[which(r[,s]==g+1),s];bb=W[which(r[,s]==g-1),s]}
          W[which(r[,s]==g),s]=rtruncnorm(1,a=aa,b=bb,mean=Mu[which(r[,s]==g),m]+alpha2[m]*xi2[which(r[,s]==g),m],sd=sqrt(Sigma_s2[s,m])) 
        }
        W[which(r[,s]==(n_T[s]+1)),s]=rtruncnorm(sum(r[,s]==n_T[s]+1,na.rm=T),b=W[which(r[,s]==n_T[s]),s],
                                                 mean=Mu[which(r[,s]==n_T[s]+1),m]+alpha2[m]*xi2[which(r[,s]==n_T[s]+1),m],sd=sqrt(Sigma_s2[s,m])) 
      }
    }
  }
  post.mean.mu=apply(Mu[,ceiling(M/2):M],1,mean)
  return(post.mean.mu)
}

##########      Parameter setup      ##########
#no. of genes involved in the studies
G=500
#no.studies
S=6
#no.platforms
#P=2
#percentage of genes involved in each study
lambda=rep(0.8,S)
#precision parameters
rho=runif(S,min=0.3,max=0.9)
sigma_p2=c(sigma2_rho(0.9)*0.6,sigma2_rho(0.9)*0.8) #platform
sigma_s2=sigma2_rho(rho)-rep(sigma_p2,c(S/2,S/2))#study
#sigma_s2+rep(sigma_p2,c(S/2,S/2))==sigma2_rho(rho)
#no. of ranked genes
p_T=rep(0.1,S)
n_T=round(p_T*G)

#no. of methods evaluated
m=13
##############################################
#intitialize matrices to store results
result=matrix(NA,G,m)
dimnames(result)[[1]]=1:G
eval=matrix(NA,4,m)
dimnames(eval)[[1]]=c("cov_10","cov_20","cov_50","cov_100")
dimnames(result)[[2]]=dimnames(eval)[[2]]=
  c("MEAN","GEO","MED","RRA","Stuart","MC1","MC2","MC3","BARD","BIRRA","BAYES","CEMC.k","CEMC.s")
##############################################

#########      Simulate data       ##########
mu_g=rnorm(G)
names(mu_g)=1:G
true_rank=rank(-mu_g)
kappa=cbind(rnorm(G,mean=0,sd=sqrt(sigma_p2[1])),rnorm(G,mean=0,sd=sqrt(sigma_p2[2])))
## simulate data for each study ##
ss=1 ## initialize s to gurantee each gene is included in at least one study
while(ss>0)
{
  w=r=r.k=matrix(NA,G,S)
  rownames(w)=rownames(r)=rownames(r.k)=1:G
  colnames(r)=colnames(r.k)=1:S
  for (s in 1:(S/2))
  {
    index=which(rbinom(G,1,lambda[s])==1)
    n_T[s]=min(n_T[s],length(index))
    w[index,s]=mu_g[index]+rnorm(length(index),mean=0,sd=sqrt(sigma_s2[s]))+kappa[index,1]
    r[index,s]=r.k[index,s]=rank(-w[index,s])
    r[is.na(r[,s])==F&r[,s]>n_T[s],s]=n_T[s]+1
    r.k[is.na(r[,s])==F&r[,s]>n_T[s],s]=NA
  }
  for (s in (S/2+1):S)
  {
    index=which(rbinom(G,1,lambda[s])==1)
    n_T[s]=min(n_T[s],length(index))
    w[index,s]=mu_g[index]+rnorm(length(index),mean=0,sd=sqrt(sigma_s2[s]))+kappa[index,2]
    r[index,s]=r.k[index,s]=rank(-w[index,s])
    r[is.na(r[,s])==F&r[,s]>n_T[s],s]=n_T[s]+1
    r.k[is.na(r[,s])==F&r[,s]>n_T[s],s]=NA
  }
  ss=sum(apply(is.na(w),1,sum)==S)
}  
ind=apply(r.k,1,function(x) all(is.na(x)))
r_k=r.k[!ind, ]
space=lapply(as.list(data.frame(r)),function(x)which(!is.na(x)))
r.k.list=lapply(as.list(data.frame(r.k)),order,na.last=NA)
true_rank.k=rank(true_rank[rownames(r_k)])

##############################################
result[,"BAYES"]=rank(-BiG_gamma(r=r,G=G,S=S))
#############      RRA methods   ###############
sr=t(t(r)/apply(r,2,max,na.rm=T)) ## normalize rank
result[,"RRA"]=order(as.numeric(as.vector(aggregateRanks(rmat = sr, method = "RRA")[,1])))
result[,"Stuart"]=order(as.numeric(as.vector(aggregateRanks(rmat = sr, method = "stuart")[,1])))

## Borda
result[,"MEAN"]=rank(apply(r,1,mean,na.rm=T))
result[,"GEO"]=rank(apply(r,1,geo.mean,na.rm=T))
result[,"MED"]=rank(apply(r,1,median,na.rm=T))

## MC
MCO=MC(input=r.k.list,space=space)
ind=as.numeric(names(true_rank.k))
result[ind,"MC1"]=order(MCO$MC1.TopK)
result[ind,"MC2"]=order(MCO$MC2.TopK)
result[ind,"MC3"]=order(MCO$MC3.TopK)

## CEMC
result[ind,"CEMC.s"]=order(CEMC(input=r.k.list,space=space,dm="s")$TopK)
result[ind,"CEMC.k"]=order(CEMC(input=r.k.list,space=space,dm="k")$TopK)

###############     BIRRA      ################# 
result[,"BIRRA"]=rank(-BIRRA(r))
## BARD
result[,"BARD"]=rank(-read.table(file=paste0("bard_calc_",k,"/RhoBar.txt"),header = T,row.names = 1)[,2]) ## "RhoBar.txt" obtained from C++ program 

##############################################
############## performance evaluation ############
eval["cov_10",]=apply(matrix(true_rank[apply(result,2,order)[1:10,]]<=10,nrow=10),2,mean)
eval["cov_20",]=apply(matrix(true_rank[apply(result,2,order)[1:20,]]<=20,nrow=20),2,mean)
eval["cov_50",]=apply(matrix(true_rank[apply(result,2,order)[1:50,]]<=50,nrow=50),2,mean)
eval["cov_100",]=apply(matrix(true_rank[apply(result,2,order)[1:100,]]<=100,nrow=100),2,mean)
