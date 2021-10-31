# install.packages("truncnorm")
# install.packages("RobustRankAggreg")
library(truncnorm)
library(RobustRankAggreg)
tau_rho=function(rho){1/(rho^(-2)-1)}

start=proc.time()
##########      Parameter setup      ##########
#no. of genes involved in the studies
I=100
#no.studies
J=5
#percentage of genes involved in each study
lambda=0.8
#precision parameters
rho=0.8
tau=rep(tau_rho(rho),J)
#no. of ranked genes
p_T=rep(0.2,J)
n_T=round(p_T*I)
#no. of MCMC runs
M=10000
#no. of simulated data sets for the same parameter setting
K=1
#set hyper-parameters
alpha=1.17
beta=0.65
delta=0.1

## initialize matrices to store results 
Bayes_rank=post.median.theta=post.mean.theta=true_rank=theta=matrix(0,I,K)
rownames(true_rank)=1:I
coverage_Bayes_10=coverage_Bayes_20=cor_20_Bayes=cor_10_Bayes=cor_Bayes=matrix(NA,1,K)
Theta=array(0,dim=c(I,M,K))
Tau=array(tau_rho(rho),dim=c(J,M,K))
post.median.tau=post.mean.tau=matrix(0,J,K)

##############################################

##########      Simulate data       ##########
set.seed(113)
for (k in 1: K)
{
  theta[,k]=rnorm(I)
  names(theta[,k])=1:I
  true_rank[,k]=rank(-theta[,k])
  ## simulate data for each study ##
  s=1
  while(s>0)
  {
    Z=z=y=matrix(NA,I,J)
    rownames(z)=rownames(y)=1:I
    colnames(y)=1:J
    for (j in 1:J)
    {
      index=which(rbinom(I,1,lambda)==1)
      n_T[j]=min(n_T[j],length(index))
      z[index,j]=theta[index,k]+rnorm(length(index),mean=0,sd=sqrt(1/tau[j]))
      y[index,j]=rank(-z[index,j])
      y[is.na(y[,j])==F&y[,j]>n_T[j],j]=n_T[j]+1
    }
    s=sum(apply(is.na(z),1,sum)==J)
  }
  ##############################################
  
  ##########      Bayesian model      ##########
  # var.j=apply(z,2,var,na.rm=T)
  # tau.j=1/(var.j-1)
  Z=z
  Theta[,1,k]=theta[,k]
  mu=nu=numeric(M)
  mu[1]=tau_rho(rho)
  nu[1]=tau_rho(rho)^2/0.25
  for (m in 2:M)
  {
    for (i in 1:I)
    {
      c=sum(Z[i,is.na(z[i,])==F]*Tau[is.na(z[i,])==F,m-1,k])
      d=sum(Tau[is.na(z[i,])==F,m-1,k])+1
      Theta[i,m,k]=rnorm(1,c/d,sqrt(1/d))
    }
    for (j in 1:J)
    {
      Tau[j,m,k]=rgamma(1,shape=nu[m-1]+sum(is.na(z[,j])==F)/2,
                        rate=nu[m-1]/mu[m-1]+sum((Z[,j]-Theta[,m,k])^2,na.rm=T)/2)
    }
    h=0
    while (h==0)
    {
      mu_star=1/rgamma(1,shape=J*nu[m-1]-1,rate=nu[m-1]*sum(Tau[,m,k])) #J*nu-1>0
      p=min(1,exp(-delta*(mu_star-mu[m-1])))
      mu[m]=mu_star
      h=(runif(1)<p)
    }
    l=0
    while(l==0)
    {
      nu_star=exp(log(nu[m-1])+rnorm(1,0,1.1))
      while(nu_star<(1/J))
      {nu_star=exp(log(nu[m-1])+rnorm(1,0,1.1))}
      p=min(1,((nu_star^-alpha)*exp(-beta*(1/nu_star-1/nu[m-1]))/(nu[m-1]^-alpha))*
              exp(sum(nu_star*log(Tau[,m,k]*nu_star/mu[m])-Tau[,m,k]*(nu_star-nu[m-1])/mu[m-1]+lgamma(nu[m-1])
                      -lgamma(nu_star)-log(Tau[,m,k]*nu[m-1]/mu[m])*nu[m-1])))
      nu[m]=nu_star
      l=(runif(1)<p)
    }
    for (l in 1:10)
    {
      for (j in 1:J)
      {
        index=sample(1:n_T[j])
        for (i in index)
        {
          if (i==1) {aa=Z[which(y[,j]==2),j];bb=Inf} 
          else if (i==n_T[j]) {aa=-Inf;bb=Z[which(y[,j]==(n_T[j]-1)),j]} 
          else {aa=Z[which(y[,j]==i+1),j];bb=Z[which(y[,j]==i-1),j]}
          Z[which(y[,j]==i),j]=rtruncnorm(1,a=aa,b=bb,mean=Theta[which(y[,j]==i),m,k],sd=sqrt(1/Tau[j,m,k])) 
        }
      }
    }
  }
  ##############################################
  post.mean.theta[,k]=apply(Theta[,ceiling(M/2):M,k],1,mean)
  post.median.theta[,k]=apply(Theta[,ceiling(M/2):M,k],1,median)
  post.mean.tau[,k]=apply(Tau[,ceiling(M/2):M,k],1,mean)
  post.median.tau[,k]=apply(Tau[,ceiling(M/2):M,k],1,median)
  Bayes_rank[,k]=rank(-post.mean.theta[,k])
  coverage_Bayes_10[k]=mean(true_rank[Bayes_rank[,k]<=10,k]<=10)
  coverage_Bayes_20[k]=mean(true_rank[Bayes_rank[,k]<=20,k]<=20)
  cor_10_Bayes[k]=cor(true_rank[true_rank[,k]<=10,k],Bayes_rank[true_rank[,k]<=10,k])
  cor_20_Bayes[k]=cor(true_rank[true_rank[,k]<=20,k],Bayes_rank[true_rank[,k]<=20,k])
  cor_Bayes[k]=cor(true_rank[,k],Bayes_rank[,k])
  ##############################################
}
duration=proc.time()-start
setwd("D:/research/Bayes")
save.image(file = "full.RData")
pdf("full.pdf",width=12)
par(mfrow=c(5,2),mai=c(0.25,0.5,0.15,0.1))
for (j in 1:5) 
{
  plot(Tau[j,,1],type="l")
  plot(1/Tau[j,,1],type="l")
}
par(mfrow=c(3,1),mai=c(0.3,0.6,0.15,0.1))
plot(Theta[1,,1],type="l")
plot(mu,type="l")
plot(nu,type="l")
dev.off()

plot(nu[1:1500],type="l")

load("full.RData")
pdf("D:\\Dropbox\\Research_Lily_Sherry\\Report\\Presentations\\fullTrace_p.pdf",width=8)
par(mfrow=c(5,1),mai=c(0.25,0.6,0.2,0.2))
plot(Theta[1,,1],type="l",main="Trace Plots from Full Model")
plot(1/Tau[1,,1],type="l")
plot(1/Tau[5,,1],type="l")
plot(mu,type="l")
plot(nu,type="l")
dev.off()
