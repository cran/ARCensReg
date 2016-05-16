# library(mvtnorm)
# library(mnormt)
# require("msm")
# library(tmvtnorm)
# require("MASS")
# library(numDeriv)
#
# source('auxfunctions.R')

################################################################################
## Algoritmo SAEM ##
################################################################################
#cc must be 0 for miss

SAEM<-function(cc,y,x,p=1,M=10,cens='left',perc=0.25,MaxIter=400,pc=0.18,x_pred=NULL,miss=NULL,tol=0.0001){
  pb <- txtProgressBar(min = 0, max = MaxIter, style = 3)
  #valores iniciais
  m<-length(y)
  if (length(miss)==0) {
    beta1 = solve(t(x)%*%x)%*%t(x)%*%y
    l<-length(beta1)
    pi1 = as.numeric(pacf((y - x%*%beta1),lag.max = p,plot=F)$acf)
    pi1 = optim(pi1,lcc,y=y,x=x,lower = rep(-.999,p), upper = rep(0.999,p),method='L-BFGS-B')$par
    phi1 = estphit(pi1)
    lambda = matrix(c(-1,phi1))
    beta1 = solve(t(x)%*%solve(MatArp(phi1,m))%*%x)%*%t(x)%*%solve(MatArp(phi1,m))%*%y
    spi = t(lambda)%*%Dbeta(beta1,y,x,p)%*%lambda
    sigmae = as.numeric(spi/m)
  }
  else
  {
    beta1<-solve(t(x[-miss,])%*%x[-miss,])%*%t(x[-miss,])%*%y[-miss]
    l<-length(beta1)
    pi1 = as.numeric(pacf((y - x%*%beta1)[-miss],lag.max=p,plot=F)$acf)
    pi1 = optim(pi1,lcc,y=y[-miss],x=as.matrix(x[-miss,]),lower = rep(-.999,p), upper = rep(0.999,p),method='L-BFGS-B')$par
    phi1 = estphit(pi1)
    lambda = matrix(c(-1,phi1))
    beta1 = solve(t(as.matrix(x[-miss,]))%*%solve(MatArp(phi1,m-length(miss)))%*%as.matrix(x[-miss,]))%*%
      t(as.matrix(x[-miss,]))%*%solve(MatArp(phi1,m-length(miss)))%*%as.matrix(y[-miss])
    spi = t(lambda)%*%Dbeta(beta1,y[-miss],as.matrix(x[-miss,]),p)%*%lambda
    sigmae = as.numeric(spi/(m-length(miss)))
    cc[miss] = rep(0,length(miss))
  }
  teta <- c(beta1,sigmae,phi1)
  r = length(teta)
  tempoi = Sys.time()

  if ((sum(cc)==0)&(length(miss)==0)) {
    teta1 <- teta
    media = x%*%beta1
    V<-MatArp(phi1,m)
    Psi<-sigmae*V
    #
    H = -Ht(teta1,y,x)
    ep_par = sqrt(diag(solve(H)))
    logver <- log(LogVerosCensLeft(cc,y,media,Psi)$ver)
    npar<-length(c(teta1))
    loglik<-logver
    AICc<- -2*loglik +2*npar
    AICcorr<- AICc + ((2*npar*(npar+1))/(m-npar-1))
    BICc <- -2*loglik +log(m)*npar
    tempof=Sys.time()
    dift =tempof-tempoi
    setTxtProgressBar(pb, MaxIter)

    if (!is.null(x_pred)) {
      h = nrow(x_pred)
      sig_pred = MatArp(phi1,m+h) *sigmae
      pred = x_pred%*%beta1 + sig_pred[m+1:h,1:m]%*%solve(sig_pred[1:m,1:m])%*%(y-x%*%beta1)
      obj.out <- list(beta1 = beta1, sigmae = sigmae, phi1 = phi1,pi1=pi1, loglik=loglik,
                      AIC=AICc, BIC=BICc, AICcorr=AICcorr,theta = teta1,ep = ep_par,
                      pred=pred,timediff=dift)
    }
    else {
      obj.out <- list(beta1 = beta1, sigmae = sigmae, phi1 = phi1,pi1=pi1, loglik=loglik,
                      AIC=AICc, BIC=BICc, AICcorr=AICcorr,theta = teta1,ep = ep_par,timediff=dift)
    }

    class(obj.out) <- "NCens"

  }

  else {

    MG <- round(M/(1-perc),0)
    M0 <- MG - M #burn in

    Theta <- matrix(NA,MaxIter,length(teta))

    criterio <- 10000
    count <- 0

    media= x%*%beta1
    V<-MatArp(phi1,m)
    Psi<-sigmae*V

    ## sequencia decrescente de numero positivos para o parametro de suaviazacao do SAEM
    ####pc = ponto de corte
    if(pc==1)
    {
      seqq=rep(1,pc*MaxIter)
    } else
    {
      seqq = c(rep(1,pc*MaxIter),(1/((((pc*MaxIter)+1):MaxIter)-(pc*MaxIter))))
      seqq = c(rep(1,MaxIter-length(seqq)),seqq)
    }

    SAEM_ss <- array(data=0,dim=c(MaxIter+1))
    SAEM_sig = array(data=0,dim=c(MaxIter+1))
    SAEM_xx <- array(data=0,dim=c(MaxIter+1,l,l))
    SAEM_xy <- array(data=0,dim=c(MaxIter+1,l))
    SAEM_y = array(data=0,dim=c(MaxIter+1,sum(cc)))
    SAEM_D = array(data=0,dim=c(MaxIter+1,p+1,p+1))
    SAEM_delta = array(data=0,dim=c(MaxIter+1,r))
    SAEM_G = array(data=0,dim=c(MaxIter+1,r,r))
    SAEM_H = array(data=0,dim=c(MaxIter+1,r,r))

    while(criterio > tol & count<MaxIter){

      count <- count + 1
      setTxtProgressBar(pb, count)

      ## Passo de Simulacao: Gera das distribuicoes condicionis
      t1 <- y
      gibbs <- amostradordegibbs(MG,M0,m,t1,cc,y,media,Psi,miss,cens)
      amostragibbs <- gibbs$amostragibbs
      uyi <- matrix(amostragibbs[,1:m],nrow=M,ncol=m)

      ## Passo de Aproximacao
      somaD = matrix(0,p+1,p+1)
      somass <- 0
      somaxx <- matrix(0,l,l)
      somaxy <- matrix(0,l,1)
      somay <- matrix(0,sum(cc),1)
      somadelta = matrix(0,r,1)
      somaG = matrix(0,r,r)
      for (k in 1:M) {
        yi <- matrix(uyi[k,],nrow=m,ncol=1)
        somass <- somass + (t(yi-media)%*%solve(V)%*%(yi-media))
        somaxx <- somaxx + (t(x)%*%solve(V)%*%x)
        somaxy <- somaxy + t(x)%*%solve(V)%*%(yi)
        somay <- somay + as.matrix(yi[cc==1,])
        J = Jt(teta,yi,x)
        H = Ht(teta,yi,x)
        somadelta = somadelta + J
        somaG = somaG + (-H - (J)%*%t(J))
        somaD = somaD + Dbeta(beta1,yi,x,p)
      }
      E_D = 1/M*somaD
      E_ss = (1/M)*somass
      E_xx = (1/M)*somaxx
      E_xy = (1/M)*somaxy
      E_y = (1/M)*somay
      E_delta = 1/M * somadelta
      E_G = 1/M*somaG

      ## Aproximacao Estocastica
      SAEM_D[count+1,,] = SAEM_D[count,,] + seqq[count]*(E_D - SAEM_D[count,,])
      SAEM_ss[count+1] = SAEM_ss[count] + seqq[count]*(E_ss - SAEM_ss[count])
      SAEM_xx[count+1,,] = SAEM_xx[count,,] + seqq[count]*(E_xx - SAEM_xx[count,,])
      SAEM_xy[count+1,] = SAEM_xy[count,] + seqq[count]*(E_xy - SAEM_xy[count,])
      SAEM_y[count+1,] = SAEM_y[count,] + seqq[count]*(E_y - SAEM_y[count,])
      SAEM_delta[count+1,] = SAEM_delta[count,] + seqq[count]*(E_delta - SAEM_delta[count,])
      SAEM_G[count+1,,] = SAEM_G[count,,] + seqq[count]*(E_G - SAEM_G[count,,])
      SAEM_H[count+1,,] = SAEM_G[count+1,,] - (SAEM_delta[count+1,])%*%t(SAEM_delta[count+1,])
      ## Passo M
      beta1<- solve(SAEM_xx[count+1,,])%*%SAEM_xy[count+1,]
      media<-x%*%beta1

      sigmae<- (1/m)*(SAEM_ss[count+1])
      sigmae<-as.numeric(sigmae)

      pi1 = optim(pi1,lc,lower = rep(-.999,p), upper = rep(0.999,p), n=m, D = SAEM_D[count+1,,],method='L-BFGS-B')$par
      phi1 = estphit(pi1)
      teta1 <- c(beta1,sigmae,phi1)

      V<-MatArp(phi1,m)
      Psi<-sigmae*V

      if (count>1){
        criterio <- sqrt((teta1/teta-1)%*%(teta1/teta-1))
      }

      Theta[count,] <- teta1
      teta<-teta1
    }
    setTxtProgressBar(pb, MaxIter)
    Theta=Theta[1:count,]

    if (cens=='left') {
      logver <- ifelse(length(miss)>0,log(LogVerosCensLeft(as.matrix(cc[-miss]),as.matrix(y[-miss]),as.matrix(media[-miss,]),
                                                           Psi[-miss,-miss])$ver),log(LogVerosCensLeft(cc,y,media,Psi)$ver))
    }
    else {
      logver <- ifelse(length(miss)>0,log(LogVerosCensRig(cc[-miss],y[-miss],media[-miss,],
                                                          Psi[-miss,-miss])$ver),log(LogVerosCensRig(cc,y,media,Psi)$ver))
    }
    tempof=Sys.time()
    dift = tempof - tempoi#difftime(tempof,tempoi,units='mins')
    npar<-length(c(teta1))
    loglik<-logver
    AICc<- -2*loglik +2*npar
    AICcorr<- AICc + ((2*npar*(npar+1))/(m-npar-1))
    BICc <- -2*loglik +log(m)*npar
    ep_par = sqrt(diag(solve(SAEM_H[count+1,,])))
    ###
    yest = numeric(m)
    yest[cc==0] = y[cc==0]
    yest[cc==1] = SAEM_y[count+1,]
    ###
    if (length(miss)>0) {
      yest[miss] = yi[miss]#apply(uyi,2,mean)[miss]
    }
    ####previsao
    if (!is.null(x_pred)) {
      h = nrow(x_pred)
      sig_pred = MatArp(phi1,m+h) *sigmae
      pred = x_pred%*%beta1 + sig_pred[m+1:h,1:m]%*%solve(sig_pred[1:m,1:m])%*%(yest-x%*%beta1)
      obj.out <- list(beta1 = beta1, sigmae = sigmae, phi1 = phi1,pi1=pi1, Theta=Theta, loglik=loglik,
                      AIC=AICc, BIC=BICc, AICcorr=AICcorr,theta = teta1,iter=count,ep = ep_par,
                      pred=pred,yest=yest,timediff = dift,criteria=criterio)
    }
    else {
      obj.out <- list(beta1 = beta1, sigmae = sigmae, phi1 = phi1,pi1=pi1, Theta=Theta, loglik=loglik,
                      AIC=AICc, BIC=BICc, AICcorr=AICcorr,theta = teta1,iter=count,ep = ep_par,
                      yest=yest,timediff = dift,criteria=criterio)
    }

    class(obj.out) <- "SAEM_Cens"

  }
  return(obj.out)
}