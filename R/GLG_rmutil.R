###
#
# Investigate rmulti package for generalized gamma
#
###
inputDir<-file.path("H:","Rutter","Orient","20230705","AdenomaRisk")
headDir<-file.path("H:","CRC-SPIN","NatHistMod","AdRisk")
setwd(headDir)

library("rmutil") #needed for rggamma()
library("MASS") #needed for mvrnorm()

#sink(file.path(headDir,"GLGmoments_log.txt"))
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
quickHist<-function(X){
   hist(X,breaks=seq(from=min(X),to=max(X),length=200)) 
}
nSmpl<-1000
# shape: s; sigma (doc); k (Stacy 1962)
# scale: m; mu (doc); sigma/mu=a (stacy 1962)
# family: f (f=1 ~Beta); nu (doc); b (Stacy 1962)
k<-40
m<-1
b<-2
Y<-rggamma(n=1000,s=k,m=m,f=b)   
#quickHist(Y)
c(mean(Y),var(Y))
c(mean(log(Y)),var(log(Y)))
c(digamma(k)/b+log(m/k),trigamma(k)/b^2)

alpha<-log()

