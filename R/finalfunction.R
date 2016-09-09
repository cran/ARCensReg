#source('princfunction.R')

ARCensReg = function(cc,y,x,p=1,cens='left',x_pred=NULL,miss=NULL,tol=0.0001,show.convergence=TRUE,M=10,perc=0.25,MaxIter=400,pc=0.18,show_se=TRUE)
{
  m = length(y)

  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.numeric(x)) stop("x must be a numeric matrix")
  if (!is.matrix(x)) x=as.matrix(x)
  if (det(t(x)%*%x)==0) stop("the columns of x must be linearly independent")
  ## Verify error at parameters specification

  if (cens!='left'& cens!='right') stop('cens must be left or right')
  #No data
  if( (length(x) == 0) | (length(y) == 0) | (length(cc) == 0)) stop("All parameters must be provided.")

  #Validating if exists NA's
  if (length(miss)>0) {
    cc[miss] = 0
    if (sum(is.na(y[-miss]))>0) stop("NA values in y must be specified in the argument miss")
    if (sum(is.na(y[miss]))!=length(miss)) cat("\n","WARNING: y[miss] is not NA","\n","\n")
  }
  else {
    if (sum(is.na(y))>0) stop("NA values in y must be specified in the argument miss")
  }
  if (sum(cc%in%c(0,1))< length(cc)) stop("The elements of the vector cc must be 0 or 1")
  if(sum(is.na(x)) > 0) stop("There are some NA values in x")
  if (sum(is.na(cc)) > 0) stop("There are some NA values in cc")


  #Validating dims data set
  if (ncol(as.matrix(y)) > 1) stop("y must have just one column")
  if (ncol(as.matrix(cc)) > 1) stop("cc must have just one column")
  if( length(y) != nrow(as.matrix(x)) ) stop("x does not have the same number of lines than y")
  if( length(y) != length(cc) ) stop("cc does not have the same length than y")

  if (!is.null(x_pred)) {
    x_pred = as.matrix(x_pred)
    if (ncol(x_pred)!=ncol(as.matrix(x))) stop("x_pred must have the same number of columns than x")
    if (sum(is.na(x_pred))>0) stop("There are some NA values in x_pred")
    if (!is.numeric(x_pred)) stop("x_pred must be a numeric matrix")
  }

  if (sum(miss %in% 1:m)<length(miss) ) stop("miss must indicate the index of the missing data on y")

  #Validating supports
  if(length(p)!=1) stop("p must be a positive integer value")
  if(!is.numeric(p)) stop("p must be a positive integer value")
  if(p!= round(p)|p<=0) stop("p must be a positive integer value")
  if(tol <= 0) stop("tolerance must be a positive value (suggested to be small)")
  if(!is.numeric(MaxIter)) stop("MaxIter must be a positive integer value")
  if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
  if(!is.numeric(M)) stop("M must be a positive integer value")
  if(M <= 1 |M%%1!=0) stop("M must be a positive integer value (greater than 1)")
  if(!is.numeric(pc)) stop("pc must be a real number in [0,1]")
  if(pc > 1 | pc < 0) stop("pc must be a real number in [0,1]")
  if(!is.numeric(perc)) stop("perc must be a real number in [0,1)")
  if(perc >= 1 | perc < 0) stop("perc must be a real number in [0,1)")
  if(!is.logical(show.convergence)) stop("show.convergence must be TRUE or FALSE.")
  if(!is.logical(show_se)) stop("show_se must be TRUE or FALSE.")

  #Load required libraries

  #Running the algorithm
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')

  out <-suppressWarnings(SAEM(cc,y,x,p,M=M,cens=cens,perc=perc,MaxIter=MaxIter,pc = pc,x_pred = x_pred,miss = miss,tol=tol,show_ep=show_se))

  cat('\n\n')
  cat('---------------------------------------------------\n')
  cat('  Censored Linear Regression Model with AR Errors \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat("p =",p)
  cat('\n')
  cat('---------\n')
  cat('Estimates\n')
  cat('---------\n')
  cat('\n')
  l = ncol(x)
  lab = numeric(p +l +1)
  for (i in 1:ncol(x)) lab[i] = paste('beta',i-1,sep='')
  lab[l+1] = 'sigma2'
  for (i in ((l+2):length(lab))) lab[i] = paste('phi',i-l-1,sep='')
  if (show_se) {
    tab = round(rbind(out$theta,out$ep),4)
    colnames(tab) = lab
    rownames(tab) = c("","s.e.")
  }
  else {
    tab = round(rbind(out$theta),4)
    colnames(tab) = lab
    rownames(tab) = c("")
  }
  print(tab)
  cat('\n')
  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  critFin <- c(out$loglik, out$AIC, out$BIC, out$AICcorr)
  critFin <- round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","AICcorr"))
  print(critFin)
  cat('\n')
  cat('-------\n')
  cat('Details\n')
  cat('-------\n')
  cat('\n')
  cat('Type of censoring =',cens)
  cat('\n')
  cat('Number of missing values =',ifelse(is.null(miss),0,length(miss)))
  cat('\n')
  if (sum(cc)>0) {
  cat("Convergence reached? =",(out$iter < MaxIter))
  cat('\n')
  cat('Iterations =',out$iter,"/",MaxIter)
  cat('\n')
  cat('MC sample =',M)
  cat('\n')
  cat('Cut point =',pc)
  cat('\n')
  }
  cat("Processing time =",out$timediff,units(out$timediff))
  cat('\n','\n')

  if (sum(cc)==0) show.convergence = FALSE
  if(show.convergence)
  {
    cpl = pc*MaxIter
    npar   = l+1+p
    labels = list()
    for(i in 1:l){labels[[i]] = bquote(beta[.(i-1)])}
    labels[[l+1]] = bquote(sigma^2)
    for(i in 1:p){labels[[i+l+1]] = bquote(phi[.(i)])}

    par(mar=c(4, 4.5, 1, 0.5))
    op <- suppressWarnings(par(mfrow=c(ifelse(npar%%3==0,npar%/%3,(npar%/%3)+1),3)))

    for(i in 1:npar)
    {
      suppressWarnings(plot.ts(out$Theta[,i],xlab="Iteration",ylab=labels[[i]]))
      abline(v=cpl,lty=2)
    }
    par(mfrow=c(1,1))
  }


  if (!is.null(x_pred)) res = list(beta = out$beta,sigma2= out$sigmae,phi = out$phi1,pi1=out$pi1,theta =out$theta,SE=out$ep,
                                   loglik=out$loglik,AIC=out$AIC,BIC=out$BIC,AICcorr=out$AICcorr,time = out$timediff,pred=out$pred,criteria = out$criteria)
  else res = list(beta = out$beta,sigma= out$sigmae,phi = out$phi1,pi1=out$pi1,theta =out$theta,SE=out$ep,
            loglik=out$loglik,AIC=out$AIC,BIC=out$BIC,AICcorr=out$AICcorr,time = out$timediff,criteria = out$criteria)
  if (sum(cc)==0) {
    obj.out = list(res = res)
  }
  else obj.out = list(res = res,yest=out$yest,yyest=out$yyest,iter = out$iter)#,conv=out$Theta
  class(obj.out)  = ifelse(sum(cc)==0,'ARp-LRM','ARp-CRM')
  return(obj.out)
}
