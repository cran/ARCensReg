
###matriz de covariancia (MatArp * sig2)

MatArp<-function(phi,n) {
  p = length(phi)
  if (n==1) Rn = 1
  else Rn= toeplitz(ARMAacf(ar=phi, ma=0, lag.max = n-1))
  rhos = ARMAacf(ar=phi, ma=0, lag.max = p)[(1:p)+1]
  return(Rn/(1-sum(rhos*phi)))
}

###############phi ----- pi
estphit = function(pit) {
  p = length(pit)
  Phi = matrix(0,ncol=p,nrow=p)
  if (p>1) {
    diag(Phi) = pit
    for (j in 2:p) {
      for (k in 1:(j-1)) {
        Phi[j,k] = Phi[j-1,k] - pit[j]*Phi[j-1,j-k]
      }
    }
    return(Phi[p,])
  }
  else return(pit)
}


############estimar pi - caso com censura
lc = function(pi,D,n) {
  phi =estphit(pi)
  p = length(phi)
  lambda = matrix(c(-1,phi))
  spi = t(lambda)%*%D%*%lambda
  gp = 1
  for (i in 1:p) gp = gp*((1-pi[i]^2)^(-i))
  sig2hat = spi/n
  l = as.numeric(-n/2*log(sig2hat) -1/2*log(gp))
  return(-l)
}

############estimar pi - caso sem censura
lcc = function(pi,y,x) {
  n = length(y)
  phi =estphit(pi)
  p = length(phi)
  lambda = matrix(c(-1,phi))
  betahat = solve(t(x)%*%solve(MatArp(phi,n))%*%x)%*%t(x)%*%solve(MatArp(phi,n))%*%y
  spi = t(lambda)%*%Dbeta(betahat,y,x,p)%*%lambda
  gp = 1
  for (i in 1:p) gp = gp*((1-pi[i]^2)^(-i))
  l = as.numeric(-n/2*log(spi) -1/2*log(gp))
  return(-l)
}

################################################################################
## Log-likelihood ##
################################################################################

LogVerosCensLeft<-function(cc,y,media,Psi){
  m=length(cc)
  gammai=media
  
  if(sum(cc)==0){
    ver<-dmvnorm(as.vector(y),as.vector(gammai),Psi)
  }
  if(sum(cc)>0){
    if(sum(cc)==m){
      auxupper<-y-gammai
      ver<- pmvnorm(upper=c(auxupper),mean=rep(0,sum(cc)),sigma=Psi)
    }
    else{
      muc<- gammai[cc==1,]+Psi[cc==1,cc==0]%*%solve(Psi[cc==0,cc==0])%*%(y[cc==0]-gammai[cc==0,])
      Sc<- Psi[cc==1,cc==1]-Psi[cc==1,cc==0]%*%solve(Psi[cc==0,cc==0])%*%Psi[cc==0,cc==1]
      auxupper<- y[cc==1]-muc
      ver<-dmvnorm(y[cc==0],gammai[cc==0,],Psi[cc==0,cc==0])*
        (pmvnorm(upper=c(auxupper),sigma=Sc))
    }
  }
  
  obj.out <- list(ver = ver)
  return(obj.out)
  
}

LogVerosCensRig<-function(cc,y,media,Psi){  
  gammai=as.matrix(media)
  m=length(cc)
  
  if(sum(cc)==0){
    ver<-dmvnorm(as.vector(y),as.vector(gammai),Psi)
  }
  if(sum(cc)>0){
    if(sum(cc)==m){
      auxupper<-y-gammai
      ver<- 1-pmvnorm(upper=c(auxupper),mean=rep(0,sum(cc)),sigma=Psi)
    }
    else{
      muc<- gammai[cc==1,]+Psi[cc==1,cc==0]%*%solve(Psi[cc==0,cc==0])%*%(y[cc==0]-gammai[cc==0,])
      Sc<- Psi[cc==1,cc==1]-Psi[cc==1,cc==0]%*%solve(Psi[cc==0,cc==0])%*%Psi[cc==0,cc==1]
      auxupper<- y[cc==1]-muc
      ver<-dmvnorm(y[cc==0],gammai[cc==0,],Psi[cc==0,cc==0])*
        (1-pmvnorm(upper=c(auxupper),sigma=Sc))
    }
  }    
  
  obj.out <- list(ver = ver)
  return(obj.out)
  
}

################################################################################
## Amostrador Gibbs
################################################################################


amostradordegibbs <- function(M,M0,nj,t1,cc1,y1,media,Gama,miss,cens){
  
  draws <- matrix(NA,nrow=M,ncol=nj)
  draws[1,] <- t1
  
  if (length(miss)>0) {
    g= media
    t1[-miss] <- y1[-miss]
    muc <- as.vector(g[miss]+Gama[miss,-miss]%*%solve(Gama[-miss,-miss])%*%
                       (y1[-miss]-g[-miss]))
    #Sc <- Gama[miss,miss]-Gama[miss,-miss]%*%solve(Gama[-miss,-miss])%*%
    #  Gama[-miss,miss]
    #ym <- rtmvnorm(1, mean = muc, sigma = Sc, lower = rep(-Inf, length = length(muc)),
    #               upper = rep(Inf, length = length(muc)), algorithm="gibbs", thinning=2)
    y1[miss] = muc
  }
  
  if(sum(cc1)==0){  
    for(i in 2:M){  
      draws[i,] <- y1
    }
  }
  if(sum(cc1)>0 & sum(cc1)==nj){  
    for(i in 2:M){
      g <- as.vector(media)
      if (cens=='left') {
        t1 <- as.vector(rtmvnorm(1, mean = g, sigma = Gama, lower = rep(-Inf,
                        length = length(g)), upper = y1, algorithm="gibbs", thinning=2))
      }
      else {
        t1 = as.vector(rtmvnorm(1, mean = g, sigma = Gama, upper = rep(Inf, 
                        length = length(g)), lower = y1, algorithm="gibbs", thinning=2)) 
      }
      draws[i,] <- t1
    }
  }	
  if(sum(cc1)>0 & sum(cc1)<nj){
    if(sum(cc1)==1){
      for(i in 2:M){
        g <- media
        t1[cc1==0] <- y1[cc1==0]
        muc <- as.numeric(g[cc1==1]+Gama[cc1==1,cc1==0]%*%solve(Gama[cc1==0,cc1==0])%*%(y1[cc1==0]-g[cc1==0]))
        Sc <- as.numeric(Gama[cc1==1,cc1==1]-Gama[cc1==1,cc1==0]%*%solve(Gama[cc1==0,cc1==0])%*%Gama[cc1==0,cc1==1])
        yr <-  ifelse(cens=='left',rtnorm(1, mean = muc, sd=sqrt(Sc), lower=-Inf, upper=y1[cc1==1]),
                      rtnorm(1, mean = muc, sd=sqrt(Sc), upper=Inf, lower=y1[cc1==1]))
        t1[cc1==1] <- yr
        draws[i,] <- t1
      }
    }
    else{
      for(i in 2:M){
        g <- media
        t1[cc1==0] <- y1[cc1==0]
        muc <- as.vector(g[cc1==1]+Gama[cc1==1,cc1==0]%*%solve(Gama[cc1==0,cc1==0])%*%(y1[cc1==0]-g[cc1==0]))
        Sc <- Gama[cc1==1,cc1==1]-Gama[cc1==1,cc1==0]%*%solve(Gama[cc1==0,cc1==0])%*%Gama[cc1==0,cc1==1]
        if (cens == 'left') {
          yr <- rtmvnorm(1, mean = muc, sigma = Sc, lower = rep(-Inf, length = length(muc)),
                                upper = y1[cc1==1], algorithm="gibbs", thinning=2)
        }
        else {
          yr <- rtmvnorm(1, mean = muc, sigma = Sc, upper = rep(Inf, length = length(muc)),
                                lower = y1[cc1==1], algorithm="gibbs", thinning=2)
        }
        t1[cc1==1] <- yr
        draws[i,] <- t1
      }
    }
  }		
  
  # Amostra com burnin (M0)
  amostragibbs <- draws[(M0+1):M,]
  
  obj.out <- list(amostragibbs = amostragibbs)
  return(obj.out)
  
}


###############################################################################
##Derivadas
###############################################################################
aphi = function(phi) ifelse(length(phi)==1,log(MatArp(phi,length(phi))),
                            log(det(MatArp(phi,length(phi)))))
Dbeta = function(beta,y,x,p) {
  n = length(y)
  D = matrix(0,p+1,p+1)
  for (ii in 1:(p+1)) {
    for (jj in 1:(p+1)) {
      D[ii,jj] = sum((y-x%*%beta)[ii:(n+1-jj)]*(y-x%*%beta)[jj:(n+1-ii)])
    }
  }
  return(D)
}
Dphi1 = function(beta,y,xx,p) matrix(Dbeta(beta,y,xx,p)[2:(p+1),1])
Dphiphi2 = function(beta,phi,y,xx,p) (Dbeta(beta,y,xx,p)[2:(p+1),2:(p+1)])%*%phi

Jt = function(theta,y,x) {
  l=ncol(x)
  n =length(y)
  beta = matrix(theta[1:l])
  sig2 = theta[l+1]
  phi = theta[(l+2):length(theta)]
  p=length(phi)
  Mn = MatArp(phi,n)
  lambda = matrix(c(-1,phi))
  spi = t(lambda)%*%Dbeta(beta,y,x,p)%*%lambda
  dbeta = 1/sig2*(t(x)%*%solve(Mn)%*%y - t(x)%*%solve(Mn)%*%x%*%beta)
  dsig2 = -n/2/sig2 +1/2/sig2^2*spi
  da = matrix(jacobian(aphi,phi))
  dphi = -1/sig2*(-Dphi1(beta,y,x,p) + Dphiphi2(beta,phi,y,x,p))-1/2*da
  return(rbind(dbeta,dsig2,dphi))
}

Ht = function(theta,y,x) {
  l=ncol(x)
  n =length(y)
  r = length(theta)
  beta = matrix(theta[1:l])
  sig2 = theta[l+1]
  phi = theta[(l+2):r]
  p = length(phi)
  Mn = MatArp(phi,n)
  lambda = matrix(c(-1,phi))
  spi = t(lambda)%*%Dbeta(beta,y,x,p)%*%lambda
  invMn = solve(Mn)
  dbetabeta = -1/sig2*(t(x)%*%invMn%*%x)
  dsig2sig2 = n/2/sig2^2 - 1/sig2^3*spi
  daa = (hessian(aphi,phi))
  dphiphi = -1/sig2*Dbeta(beta,y,x,p)[2:(p+1),2:(p+1)] - 1/2*daa
  dbetasig2 = -1/sig2^2*(t(x)%*%invMn%*%y - t(x)%*%invMn%*%x%*%beta )
  dD1beta = jacobian(Dphi1,beta,y=y,xx=x,p=p)
  dDphibeta = jacobian(Dphiphi2,beta,phi=phi,y=y,xx=x,p=p)
  dbetaphi = 1/sig2*(dD1beta - dDphibeta)
  dphisig = 1/sig2^2*(-Dphi1(beta,y,x,p)+ Dphiphi2(beta,phi,y,x,p))
  H = matrix(0,r,r)
  H[1:l,1:l] = dbetabeta
  H[l+1,l+1] = dsig2sig2
  H[(l+2):r,(l+2):r] = dphiphi
  H[l+1,1:l] = H[1:l,l+1] = dbetasig2
  H[(l+2):r,1:l] = dbetaphi
  H[1:l,(l+2):r] = t(dbetaphi)
  H[l+1,(l+2):r] = H[(l+2):r,l+1] = dphisig
  return(H)
}
