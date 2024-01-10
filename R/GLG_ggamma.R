
library("ggamma")
quickHist<-function(X){
   hist(X,breaks=seq(from=min(X),to=max(X),length=200)) 
}

k<-400
b<-5
a<-2
X<-rggamma(n=1000,a=a,b=b,k=k)
#quickHist(X)
#c(mean(X),var(X))
Y<-log(X)
#quickHist(Y)
c(mean(Y),var(Y))
c(digamma(k)/b+log(a),trigamma(k)/b^2)

library("flexsurv")
Y<-log(rgengamma(5000, mu = log(a)+digamma(k)/b, sigma = 1/sqrt(b^2*k), Q = 1/sqrt(k)))
c(mean(Y),var(Y))
c(digamma(k)/b+log(a),trigamma(k)/b^2)

Y<-log(rgengamma(5000, mu = log(a)+digamma(k)/b, sigma = sqrt(trigamma(k)/b^2), Q = 1/sqrt(k)))
c(mean(Y),var(Y))
c(digamma(k)/b+log(a),trigamma(k)/b^2)

k<-2*k
Y<-log(rgengamma(5000, mu = log(a)+digamma(k)/b, sigma = sqrt(trigamma(k)/b^2), Q = 1/sqrt(k)))
c(mean(Y),var(Y))
c(digamma(k)/b+log(a),trigamma(k)/b^2)

k<-2*k
Y<-log(rgengamma(5000, mu = log(a)+digamma(k)/b, sigma = sqrt(trigamma(k)/b^2), Q = 1/sqrt(k)))
c(mean(Y),var(Y))
c(digamma(k)/b+log(a),trigamma(k)/b^2)

#### try setting mu & sigma and see if formulas
Y<-log(rgengamma(5000, mu = log(a)+digamma(k)/b, sigma = sqrt(b^2*k), Q = 1/sqrt(k)))
c(mean(Y),var(Y))
c(digamma(k)/b+log(a),trigamma(k)/b^2)

Y<-log(rgengamma(5000, mu = log(a)+digamma(k)/b, sigma = sqrt(trigamma(k)/b^2), Q = 1/sqrt(k)))
c(mean(Y),var(Y))
c(digamma(k)/b+log(a),trigamma(k)/b^2)

k<-2*k
Y<-log(rgengamma(5000, mu = log(a)+digamma(k)/b, sigma = sqrt(trigamma(k)/b^2), Q = 1/sqrt(k)))
c(mean(Y),var(Y))
c(digamma(k)/b+log(a),trigamma(k)/b^2)

k<-2*k
Y<-log(rgengamma(5000, mu = log(a)+digamma(k)/b, sigma = sqrt(trigamma(k)/b^2), Q = 1/sqrt(k)))
c(mean(Y),var(Y))
c(digamma(k)/b+log(a),trigamma(k)/b^2)


resAll<-NULL
for(k in 10^c(0:5)){
   X<-rggamma(n=1000,a=2,b=5,k=k)
   Y<-log(X)
   resAll<-rbind(resAll,c(mean(Y),digamma(k)/b+log(a),var(Y),trigamma(k)/b^2))
}
resAll

# q<-1/sqrt(k)
# alpha<-log(a)+digamma(k)/b
# sigma<-q/b
# Z<-(Y-alpha)/sigma*q + digamma(q^(-2))
# quickHist(Z)

