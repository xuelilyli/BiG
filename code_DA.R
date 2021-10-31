#define function to sample from truncated gamma distribution
qtruncgamma=function(a,b,p,shape,rate=1)
{qgamma(pgamma(a,shape=shape,rate=rate)+p*(pgamma(b,shape=shape,rate=rate)-pgamma(a,shape=shape,rate=rate)),shape=shape,rate=rate)}
rtruncgamma=function(n,a,b,shape,rate=1)
{u = runif(n, min = 0, max = 1)
x = qtruncgamma(a=a,b=b,p=u,shape=shape,rate=rate)
return(x)}

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
BiG_DA(r=r,G=G,S=S,M=20)
