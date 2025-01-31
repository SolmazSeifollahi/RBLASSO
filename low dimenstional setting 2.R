
library(MASS)
library(statmod)
library(tmvmixnorm)
library(glmnet)
library(coda)
source("criteria.R")

n<- 150        # total number of observations
ntr<- 50      # number of observations for training
nts<- 100       # number of observations for testing

p<- 15        #the number of covariates p=100, 200
rho<- 0      #the correlation value rho=0, 0.5

#characters of MCMC
nsample<- 1e4   # MCMC path length
burnin<- 1e3    # MCMC burn-in number
thin<- 10        # MCMC thining number

nrep<- 50       # the number of replications

# real value of parameters
tbeta<- c(rep(3, 5), rep(0, 10))
rsigma2<- 4

# Restrictions on parameters
H<- matrix(0, nrow=3, ncol=p)
H[1,1]<- 1
H[2,2]<- 1
H[3,3]<- 1
U<- rep(3.3, 3)         #upper bound
L<- rep(2.6, 3)   #lower bound

# the covariance matrix of covariates
Sb<- toeplitz(rho^(0:(p-1)))

#--------------storing data----------------
#storing beta
rebeta<- matrix(nrow = nrep, ncol=p)
unbeta<- matrix(nrow = nrep, ncol=p)

#storing sigma2
resigma2<- c(); unsigma2<- c()

#------------storing outputs--------------
MSEy<- matrix(nrow = nrep, ncol=2); colnames(MSEy)<- c("RE", "UN")
MAEy<- matrix(nrow = nrep, ncol=2); colnames(MAEy)<- c("RE", "UN")
MSEb<- matrix(nrow = nrep, ncol=2); colnames(MSEb)<- c("RE", "UN")
MAEb<- matrix(nrow = nrep, ncol=2); colnames(MAEb)<- c("RE", "UN")
recerit<- matrix(nrow = nrep, ncol=3)
uncerit<- matrix(nrow = nrep, ncol=3)

for (m in 1:nrep) {
  
  set.seed(110+m)
  # generate observations
  Xr<- mvrnorm(n , mu = rep(0, p), Sigma = Sb)
  X<- scale(Xr, center = TRUE, scale = FALSE)
  E<- rnorm(n , 0, rsigma2)
  y<- X%*%tbeta+ E
  
  # dividing training and testing sets
  nv<- sample(1:n, ntr, replace = F)
  Xtr<- X[c(nv),]
  Xts<- X[-c(nv),]
  
  ytr<- y[c(nv)]
  yts<- y[-c(nv)]
  
  #-----------------------restricted model------------------------#
  #initializing values
  beta1<- c(rep(3.2, times=3), 2, 2, rep(0.2, 10))
  Xtr1 <- scale(Xtr, center = TRUE, scale = FALSE)
  cv_fit<- cv.glmnet(Xtr1, ytr, alpha = 1, nfolds = 10, intercept = FALSE)
  lam1<- cv_fit$lambda.min 
  sigma21<- 1
  
  beta<- beta1
  lam<- lam1
  sigma2<- sigma21
  
  # storing vector or matrix
  rebetap<- matrix(0, nrow = nsample, ncol = p)
  relamp<- c()
  resigma2p<- c()
  set.seed(110+m)
  
  for (k in 1:nsample) {
    
    # generating the wi's
    rewp<- c()
    for (j in 1:p){
      rewp[j]<- ifelse(beta[j]!=0, 1/rinvgauss(1,lam*sqrt(sigma2)/abs(beta[j]),lam^2), rgamma(1,0.5,rate=(lam^2)*0.5))
    }
    
    # generating the lambda
    relamp[k]<- sqrt(rgamma(1, p+1, 1+0.5*sum(rewp)))
    
    # generating the sigma2
    srate<- crossprod(ytr-Xtr%*%beta)+sum(beta^2/rewp)
    resigma2p[k]<- 1/rgamma(1, 0.5*(ntr+p), 0.5*srate)
    
    # generating beta's
    reSS<- crossprod(Xtr)+ diag(1/rewp)
    resigmap<- solve(reSS)
    remup<- resigmap%*% crossprod(Xtr, ytr)
    rebetap[k,]<- rtmvn(n=1, Mean=remup, as.numeric(resigma2p[k])*resigmap,
                        D=H, lower=L, upper=U, int=beta, burn=10)
    
    beta<- rebetap[k,]
    lam<- relamp[k]
    sigma2<- resigma2p[k]
    
  }
  
  
  #-----------------------unrestricted model------------------------#
  #initializing values
  beta<- beta1
  lam<- lam1
  sigma2<- sigma21
  
  # storing vector or matrix
  unbetap<- matrix(0, nrow = nsample, ncol = p)
  unlamp<- c()
  unsigma2p<- c()
  set.seed(110+m)
  for (k in 1:nsample) {
    
    # generating the wi's
    unwp<- c()
    for (j in 1:p){
      unwp[j]<- ifelse(beta[j]!=0, 1/rinvgauss(1,lam*sqrt(sigma2)/abs(beta[j]),lam^2), rgamma(1,0.5,rate=(lam^2)*0.5))
    }
    
    # generating the lambda
    unlamp[k]<- sqrt(rgamma(1, p+1, 1+0.5*sum(unwp)))
    
    # generating the sigma2
    srate<- crossprod(ytr-Xtr%*%beta)+sum(beta^2/unwp)
    unsigma2p[k]<- 1/rgamma(1, 0.5*(ntr+p), 0.5*srate)
    
    # generating beta's
    unSS<- crossprod(Xtr)+ diag(1/unwp)
    unsigmap<- solve(unSS)
    unmup<- unsigmap%*% crossprod(Xtr, ytr)
    unbetap[k,]<- mvrnorm(1 , mu = unmup, Sigma = as.numeric(unsigma2p[k])*unsigmap)
    
    beta<- unbetap[k,]
    lam<- unlamp[k]
    sigma2<- unsigma2p[k]
    
  }
  re_chain <- rebetap[seq((burnin+1),nsample,by = thin),]
  un_chain <- unbetap[seq((burnin+1),nsample,by = thin),]
  
  
  for (i in 1:p) {
    
    if(quantile(re_chain[,i], probs = 0.025)<= 0 && quantile(re_chain[,i], probs = 0.975) >= 0){
      rebeta[m,i]<- 0
    }else {
      rebeta[m,i]<- median(re_chain[,i])
    }
    
    if(quantile(un_chain[,i], probs = 0.025)<= 0 && quantile(un_chain[,i], probs = 0.975) >= 0){
      unbeta[m,i]<- 0
    } else {
      unbeta[m,i]<- median(un_chain[,i])
    }
  }
  
  #predicted values
  reypred<- Xts%*%rebeta[m,]
  unypred<- Xts%*%unbeta[m,]
  
  # calculate MSEy
  MSEy[m,1]<- mean((yts-reypred)^2)
  MSEy[m,2]<- mean((yts-unypred)^2)  
 
  #calculate MAEy
  MAEy[m,1]<- mean(abs(yts-reypred))
  MAEy[m,2]<- mean(abs(yts-unypred)) 
  
  # calculate MSEb
  MSEb[m,1]<- mean((tbeta-rebeta[m,])^2)
  MSEb[m,2]<- mean((tbeta-unbeta[m,])^2)  
  
  #calculate MAEb
  MAEb[m,1]<- mean(abs(tbeta-rebeta[m,]))
  MAEb[m,2]<- mean(abs(tbeta-unbeta[m,])) 
  # calculate total_miss, fdr, fn
  recerit[m,]<- criteria(rebeta[m,], tbeta)
  uncerit[m,]<- criteria(unbeta[m,], tbeta)
}

output2<- matrix(nrow=7, ncol=2)
rownames(output2)<- c("MSEy", "MAEy", "MSEb", "MAEb",  "total_miss", "fdr", "fnr")
colnames(output2)<- c("REmedian", "UNmedian")
output2[1,]<- colMeans(MSEy, na.rm = T)
output2[2,]<- colMeans(MAEy, na.rm = T)
output2[3,]<- colMeans(MSEb, na.rm = T)
output2[4,]<- colMeans(MAEb, na.rm = T)
output2[5:7,1]<- colMeans(recerit, na.rm = T)
output2[5:7,2]<- colMeans(uncerit, na.rm = T)


round(output2, 4)


