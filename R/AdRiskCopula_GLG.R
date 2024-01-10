###
#
# Quick check on 
#
###
inputDir<-file.path("H:","Rutter","Orient","20230705","AdenomaRisk")
headDir<-file.path("H:","CRC-SPIN","NatHistMod","AdRisk")
setwd(headDir)

library("rmutil") #needed for rggamma()
library("MASS") #needed for mvrnorm()
library("flexsurv") #needed for GLG distribution; called by numAdProbs_GLG_nh_poisson in range_check.R
load(file.path(inputDir,"model_parms_2.4.2.rda"))
mu<-model_parms_2.4.2[["ar.mean"]]
sigma<-(model_parms_2.4.2[["ar.sd"]])

#sink(file.path(headDir,"GLGmoments_log.txt"))
print("#####")
print("Using Mean and Variance from previous calibration's posterior mean for Normal misture distribution")
print(paste0("Mean= ",round(mu,3),"  SD= ",round(sigma,3)))

#for(lambda in c(seq(from=0,to=sigma,by=0.2),sigma,seq(from=(floor(sigma/0.2)+1)*0.2,to=2,by=0.2))){
calcGLGmean<-function(lambda,mu,sigma){
   ifelse(lambda==0,mu,
      mu + sigma*(digamma(lambda^(-2)) + log(lambda)^(-2))/abs(lambda)
   )
}

calcGLGvar<-function(lambda,mu,sigma){
   ifelse(lambda==0,sigma^2,
      sigma^2*(trigamma(lambda^(-2)))/(lambda^2)
   )
}

#1 dimensional
nSmpl<-4000
lambda<-2
mu<-4
sigma<-1.25
rho<-0.2
simDist1D<-function(nSmpl,lambda,mu,sigma,rho){
   print("###")
   print(paste0("#Shape:",lambda))
   if(lambda==0){print("Special Case: Normal")}
   if(lambda==sigma){print("Special Case: Negative Binomial")}
   #sample
   set.seed(20230816)
   
   Z<-rnorm(nSmpl, mean = mu, sd = sigma)
   U<-pnorm(Z,mean=mu,sd=sigma)
   hist(U,breaks=seq(from=min(U),to=max(U),length=200),freq=F)
   G<-log(qgengamma(U,mu=mu,sigma=sigma,Q=lambda))
   hist(G,breaks=seq(from=min(G),to=max(G),length=200),freq=F)
   G2<-log(rgengamma(nSmpl,mu=mu,sigma=sigma,Q=lambda))
   hist(G2,breaks=seq(from=min(G2),to=max(G2),length=200),freq=F)
   
   Z2d<-mvrnorm(nSmpl,mu=c(0,0), Sigma=cbind(c(1,rho),c(rho,1)))
   Z<-Z2d[,1]
   hist(Z,breaks=seq(from=min(Z),to=max(Z),length=200),freq=F)
   Z2<-Z2d[,2]
   hist(Z2,breaks=seq(from=min(Z2),to=max(Z2),length=200),freq=F)
   
   U<-pnorm(Z,mean=0,sd=1)
   hist(U,breaks=seq(from=min(U),to=max(U),length=200),freq=F)
   U2<-pnorm(Z2,mean=0,sd=1)
   hist(U2,breaks=seq(from=min(U2),to=max(U2),length=200),freq=F)
   
   G<-qnorm(U,mean=mu,sd=sigma)
   hist(G,breaks=seq(from=min(G),to=max(G),length=200),freq=F)
   c(mean(G),sd(G))
   G2<-log(qgengamma(U2,mu=mu,sigma=sigma,Q=lambda))
   hist(G2,breaks=seq(from=min(G2),to=max(G2),length=200),freq=F)
   c(mean(G2),sd(G2))
   c(calcGLGmean(lambda,mu,sigma),calcGLGvar(lambda,mu,sigma))
   
   cor.test(Z,Z2,method="spearman")
   cor.test(U,U2,method="spearman")
   cor.test(G,G2,method="spearman")
   
   qqnorm(G, main='Normal?')
   qqline(G)
   print(shapiro.test(G))
   qqnorm(G2, main='Normal?')
   qqline(G2)
   print(shapiro.test(G2))
   
   print("apply log-transform")
   logY<-log(rgengamma(nSmpl, mu = mu, sigma = sigma, Q = lambda))
   print(paste0("Mean exp: ",round(calcGLGmean(lambda,mu,sigma),3)," vs obs: ",round(mean(logY),3)))
   print(paste0("Var exp: ",round(calcGLGvar(lambda,mu,sigma),3)," vs obs: ",round(var(logY),3)))
  
   #eval
   hist(logY,breaks=seq(from=min(logY),to=max(logY),length=200)) #Q can be a vector of shape parameters
   qqnorm(logY, main='Normal?')
   qqline(logY)
   #print(shapiro.test(logY)) #test for reject normality: H0: ~N(0,1)
   #Error in shapiro.test(logY) : sample size must be between 3 and 5000
   
   # print(data.frame(mean=mean(logY),sd=sd(logY),
   #    med=round(log(qgengamma(p=0.5, mu = mu, sigma = sigma, Q=lambda, lower.tail = TRUE, log.p = FALSE)),2),
   #    Fmu=round(pgengamma(q=exp(mu), mu = mu, sigma = sigma, Q=lambda, lower.tail = TRUE, log.p = FALSE),2)))
   print("###")
}
#sink()
simDist(nSmpl=10000,lambda=2,mu=4,sigma=1)

logY<-rggamma(n=5000,s=200,m=1,f=1) #shape s=2; scale m; family f=1 is ~Beta
hist(logY,breaks=seq(from=min(logY),to=max(logY),length=200)) #Q can be a vector of shape parameters
c(mean(logY),var(logY))
c(mean(log(logY)),var(log(logY)))
c(digamma(200)+log(1/200),trigamma(200))
c(calcGLGmean(lambda=(1/200)^2,mu,sigma),var(log(logY)))