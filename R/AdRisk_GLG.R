###
#
# Reference: Fabio et al. arxiv
#     Reparam of X~generalized gamma rv, log(x) scale, that has Normal (mu, sigma^2) correspond to lambda=0
#             "shape":lambda, "location":mu, "scale":sigma
#     Classic param have X -> log-Normal as k -> infty
#     20240104: suspect reparam uses k = 1/(lambda^2); i.e. abs(lambda) = 1/sqrt(k) and mu=0
#               and then achieve symmetry about this reference distribution by letting X->sign(lambda)*X
#               and then centering at mu: X->X+mu         
###
inputDir<-file.path("H:","Rutter","Orient","20230705","AdenomaRisk")
headDir<-file.path("H:","CRC-SPIN","NatHistMod","AdRisk","GenGammaQA")
#setwd(headDir)

library("flexsurv") #needed for GLG distribution; called by numAdProbs_GLG_nh_poisson in range_check.R
# calculate moments per Fabio et al.
calcGLGmean<-function(lambda,mu,sigma){
   mu+ifelse(lambda==0,0,
      -1*sign(lambda)*
         sigma*(digamma(lambda^(-2)) + log(lambda^(-2)))/abs(lambda)
   )
}
calcGLGvar<-function(lambda,mu,sigma){
   ifelse(lambda==0,sigma^2,
      sigma^2*(trigamma(lambda^(-2)))/(lambda^2)
   )
}

simDist<-function(nSmpl,lambda,mu,sigma,plotsOn=F){
   if(plotsOn){print(paste0("#Shape:",lambda," #mu: ",mu, " #sigma: ",sigma))}
   if(lambda==0){print("Special Case: Normal")}
   if(lambda==sigma){print("Special Case: Negative Binomial")}
   #sample
   set.seed(20240101)
  
   #print("Apply log-transform of gengamma")
   logY<-log(rgengamma(nSmpl, mu = mu, sigma = sigma, Q = lambda))
   #print(paste0("Mean calc: ",round(calcGLGmean(lambda,mu,sigma),3)," vs obs: ",round(mean(logY),3)))
   #print(paste0("Var calc: ",round(calcGLGvar(lambda,mu,sigma),3)," vs obs: ",round(var(logY),3)))
   
   if(plotsOn){
      #eval
      hist(logY,breaks=seq(from=min(logY),to=max(logY),length=200)) #Q can be a vector of shape parameters
      qqnorm(logY, main='Normal?')
      qqline(logY)
      #print(shapiro.test(logY)) #test for reject normality: H0: ~N(0,1); N in [3,5000], else error
   }
   
   data.frame(
      lambda=lambda, mu=mu, sigma=sigma,
      meanO=mean(logY),varO=var(logY),
      meanE=round(calcGLGmean(lambda,mu,sigma),3),
      varE=round(calcGLGvar(lambda,mu,sigma),3),
      medE=round(log(qgengamma(p=0.5, mu = mu, sigma = sigma, Q=lambda, lower.tail = TRUE, log.p = FALSE)),2),
      FmuE=round(pgengamma(q=mu, mu = mu, sigma = sigma, Q=lambda, lower.tail = TRUE, log.p = FALSE),2))
}

#sink(file.path(headDir,"Logs",GLGmoments_log.txt"))
mu<-0
sigma<-1
resAll<-NULL
#for(lambda in c(seq(from=-sigma,to=sigma,by=0.2),sigma,seq(from=(floor(sigma/0.2)+1)*0.2,to=2,by=0.2))){
for(lambda in seq(from=-2,to=2,by=0.2)){
   resAll<-rbind(resAll,simDist(nSmpl=5000,lambda=lambda,mu=mu,sigma=sigma))
}
resAll
with(resAll,summary(varO-varE))
with(resAll,summary(meanO-meanE))

#sink()


