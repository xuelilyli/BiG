library(truncnorm)
library(RobustRankAggreg)
tau_rho=function(rho){1/(rho^(-2)-1)} ######### funtion to calculate tau based on rho

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
n_T=n_T0=round(p_T*I)
#no. of MCMC runs
M=10000
#no. of simulated data sets for the same parameter setting (replicates)
K=1
#intitialize matrices to store results
Bayes_rank=post.median.theta=post.mean.theta=true_rank=theta=matrix(0,I,K)
rownames(true_rank)=1:I
mse.theta=coverage_Bayes_10=coverage_Bayes_20=cor_20_Bayes=cor_10_Bayes=cor_Bayes=matrix(NA,1,K)
Theta=array(0,dim=c(I,M,K))
Tau=array(tau_rho(rho),dim=c(J,M,K))
post.median.tau=post.mean.tau=matrix(0,J,K)
##############################################

# note for some of the parameter named the same, if it starts with a lower letter it is a true value; if it starts with an uppler letter, it is simulated value from the MCMC

##########      Simulate data       ##########
set.seed(113)
for (k in 1: K)
{
  theta[,k]=rnorm(I)
  names(theta[,k])=1:I
  true_rank[,k]=rank(-theta[,k])
  ## simulate data for each study ##
  s=1 ## initialize s to make sure each gene is included in at least one study
  while(s>0)
  {
    z=y=matrix(NA,I,J)
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
#   var.j=apply(z,2,var,na.rm=T)
#   tau.j=1/(var.j-1)
# initialize theta and tau
Theta[,1,k]=theta[,k]
Tau[,1,k]=tau
# initialize Z
Z=z
  for (m in 2:M)
  {
    ## update Theta
    for (i in 1:I)
    {
      c=sum(Z[i,is.na(z[i,])==F]*Tau[is.na(z[i,])==F,m-1,k])
      d=sum(Tau[is.na(z[i,])==F,m-1,k])+1
      Theta[i,m,k]=rnorm(1,c/d,sqrt(1/d))
    }
    ## update Tau
    for (j in 1:J)
    {
      Tau[j,m,k]=rgamma(1,shape=0.01+sum(is.na(z[,j])==F)/2,
                        rate=0.01+sum((Z[,j]-Theta[,m,k])^2,na.rm=T)/2)
    }
    ## update Z
    for (l in 1:10) ## repeat 10 times
    {
      for (j in 1:J)
      {
        index=sample(1:n_T[j])
        for (i in index)
        {
          if (i==1) {aa=Z[which(y[,j]==2),j];bb=Inf} 
          else if (i==n_T[j]) {aa=max(Z[which(y[,j]==(n_T[j]+1)),j]);bb=Z[which(y[,j]==(n_T[j]-1)),j]} 
          else {aa=Z[which(y[,j]==i+1),j];bb=Z[which(y[,j]==i-1),j]}
          Z[which(y[,j]==i),j]=rtruncnorm(1,a=aa,b=bb,mean=Theta[which(y[,j]==i),m,k],sd=sqrt(1/Tau[j,m,k])) 
        }
        Z[which(y[,j]==(n_T[j]+1)),j]=rtruncnorm(sum(y[,j]==n_T[j]+1,na.rm=T),b=Z[which(y[,j]==n_T[j]),j],
                                                 mean=Theta[which(y[,j]==n_T[j]+1),m,k],sd=sqrt(1/Tau[j,m,k])) 
      }
    }
  }
  ##############################################
  ## calculate Bayesian estimates and measurements for performance evaluation
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
  mse.theta[k]=mean((theta[,k]-post.mean.theta[,k])^2)
  ##############################################
}
duration=proc.time()-start
setwd("D:/research/Bayes")
save.image(file = "D:/research/Bayes/gamma.RData")
load("D:/research/Bayes/gamma.RData")
## trace plosts
pdf("gammaTrace.pdf",width=12)
par(mfrow=c(5,2),mai=c(0.25,0.5,0.15,0.1))
for (j in 1:5) 
{
  plot(Tau[j,,1],type="l")
  plot(1/Tau[j,,1],type="l")
}
par(mfrow=c(1,1),mai=rep(0.8,4))
plot(Theta[1,,1],type="l")
dev.off()
pdf("D:\\Dropbox\\Research_Lily_Sherry\\Report\\Presentations\\gammaTrace_p.pdf",width=8)
par(mfrow=c(3,1),mai=c(0.3,0.6,0.35,0.2))
plot(Theta[1,,1],type="l",main="Trace Plots from Gamma Prior")
plot(1/Tau[1,,1],type="l")
plot(1/Tau[5,,1],type="l")
dev.off()

