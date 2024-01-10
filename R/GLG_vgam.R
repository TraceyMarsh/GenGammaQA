####
# VGAM library:
# lgamma 
# prentice74 - estimation procedure, per param (3) Prentice 74

#ggamma {usefr}
#GenGamma {flexsurv}
#{ggamma}
###
library("VGAM")
nSmpl<-1000
mu<-0
sigma<-1
lambda<-200

quickHist<-function(X){
   hist(X,breaks=seq(from=min(X),to=max(X),length=200)) 
}
Z<-rlgamma(nSmpl, location=mu, scale=sigma, shape=lambda)
plot(Z,log(Y))
cor.test(Z,Y,method="spearman")

?VGAM
