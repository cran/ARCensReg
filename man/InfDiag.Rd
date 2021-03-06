\name{InfDiag}
\alias{InfDiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Influence Diagnostic in Censored Linear Regression Model with Autoregressive Errors
}
\description{
It performs influence diagnostic by a local influence approach (Cook, 1986) with three possible perturbations schemes: response perturbation (y), scale matrix perturbation (Sigma) or explanatory variable perturbation (x).
A benchmark value is calculated that depends on k.
}
\usage{
InfDiag(theta,yest,yyest,x,k=3,plots=T,indpar=rep(1,length(theta)),
                   perturbation ='y',indcolx = rep(1,ncol(x)))
}

\arguments{
 \item{theta}{
Vector of estimated parameters.
}
\item{yest}{
Vector of responses of length \code{n} with agmented data. Should be the value yest of the ARCensReg function in the case that at least one observation is censored.
}
\item{yyest}{
Should be the value yyest of the ARCensReg function in the case that at least one observation is censored. Otherwise, must be \code{y\%*\%t(y)}.
}
\item{x}{
Matrix of covariates of dimension \code{n x l}, where \code{l} is the number of fixed effects including the intercept, if considered (in models which include an intercept \code{x} should contain a column of ones).
}
\item{k}{
Constant to be used in the benchmark calculation: \code{M0+k*sd(M0)}.
}
\item{plots}{
TRUE or FALSE. Indicates if a graph should be plotted.
}
\item{indpar}{
Vector of length equal to the number of parameters, with each element 0 or 1 indicating if the respective parameter should be taking into account in the influence calculation.
}
\item{perturbation}{
Perturbation scheme. Possible values: "y" for response perturbation, "Sigma" for scale matrix perturbation or "x" for explanatory variable perturbation.
}
\item{indcolx}{
If \code{perturbation="x"}, \code{indcolx} must be a vector of length equal to the number of columns of x, with each element 0 or 1 indicating if the respective column of x should be perturbed. All columns are perturbed by default.
}
}
\details{
The function returns a vector of length n with the aggregated contribution (M0) of all eigenvectors of the matrix associated with the normal curvature.
For details see (Schumacher et. al., 2016).
}
\value{
M0
}
\author{Fernanda L. Schumacher <\email{fernandalschumacher@gmail.com}>,
Victor H. Lachos
<\email{hlachos@ime.unicamp.br}> and Christian E. Galarza <\email{cgalarza88@gmail.com}>

Maintainer: Fernanda L. Schumacher <\email{fernandalschumacher@gmail.com}>
}

\examples{
\dontrun{
#generating the data
set.seed(12341)
x = cbind(1,runif(100))
dat = rARCens(n=100,beta = c(1,-1),pit = c(.4,-.2),sig2=.5,
            x=x,cens='left',pcens=.05)
#creating an outlier
dat$data$y[40] = 5
plot.ts(dat$data$y)

#fitting the model
fit = ARCensReg(cc=dat$data$cc,y=dat$data$y,x,p=2,cens='left',
      tol=0.001,show_se=F)

#influence diagnostic
M0y = InfDiag(theta=fit$res$theta, yest=fit$yest, yyest=fit$yyest,
        x=x, k = 3.5, perturbation = "y")
M0Sigma = InfDiag(theta=fit$res$theta, yest=fit$yest, yyest=fit$yyest,
        x=x, k = 3.5, perturbation = "Sigma")
M0x = InfDiag(theta=fit$res$theta, yest=fit$yest, yyest=fit$yyest,
        x=x, k = 3.5, perturbation = "x",indcolx =c(0,1))

#perturbation on a subset of parameters
M0y1 = InfDiag(theta=fit$res$theta, yest=fit$yest, yyest=fit$yyest,
        x=x, k = 3.5, perturbation = "y",indpar=c(1,1,0,0,0))
M0y2 = InfDiag(theta=fit$res$theta, yest=fit$yest, yyest=fit$yyest,
        x=x, k = 3.5, perturbation = "y",indpar=c(0,0,1,1,1))
plot(M0y1,M0y2)
abline(v = mean(M0y1)+3.5*sd(M0y1),h = mean(M0y2)+3.5*sd(M0y2),lty=2)
}
}


\references{
Cook, R. D. (1986). Assessment of local influence. Journal of the Royal Statistical Society, Series B, 48, 133-169.

Schumacher, F. L., Lachos, V. H. & Vilca-Labra, F. E. (2016) Influence diagnostics for censored regression models with autoregressive errors. Submitted.

Zhu, H. & Lee, S. (2001). Local influence for incomplete-data models. Journal of the Royal Statistical Society, Series B, 63, 111-126.
}

\keyword{package}
\keyword{censored regression}
\keyword{autoregressive errors}
\keyword{influence diagnostic}

\seealso{\code{\link[ARCensReg]{ARCensReg}}}
