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
BiG_FB(r=r,G=G,S=S,M=20)




