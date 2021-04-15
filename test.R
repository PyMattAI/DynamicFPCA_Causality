library(far)
library(fda)
library(vars)
library(freqdom)
library(freqdom.fda)
library(pcdpca)


######### Function DFPCA.test(X, Y, num, ndx, ndy, L, q)
#
# Non-causality test in the sense of Granger between functional observations X and Y if X cause Y.
#
# Parameters:
# num: number of observations by curves 
# X: data matrix num*n with n sample size
# Y: data matrix num*n with n sample size
# ndx: number of Dynamic functional principal components used for X
# ndy: number of Dynamic functional principal components used for Y
# L: non negative integer, DFPCA filters coefficients at lags L. By default L=30.
# q: positive integer, window size for the kernel estimator of the spectral density estimator.By default, q=max(1,floor(n^0.33)). 

DFPCA.test <- function(X, Y, num, ndx, ndy, L = 30,
                        q = ceiling((dim(X)[2])^{0.33})){
  data.fd.x = Data2fd(argvals=0:(num-1)/num, X)
  filters.x = fts.dpca.filters(fts.spectral.density(data.fd.x, q = q), Ndpc = ndx, q = L)
  data.fd.y = Data2fd(argvals=0:(num-1)/num, Y)
  filters.y = fts.dpca.filters(fts.spectral.density(data.fd.y, q = q), Ndpc = ndy, q = L)
  score.x=fts.dpca.scores(data.fd.x,dpcs = filters.x)
  score.y=fts.dpca.scores(data.fd.y,dpcs = filters.y)
  colnames(score.y) <- paste("a",seq(1,ndy),sep="")
  colnames(score.x) <- paste("b",seq(1,ndx),sep="")
  var.2c <- VAR(y =cbind(score.y,score.x))
  test <- causality(var.2c,cause = paste("b",seq(1,ndx),sep=""))
  print(test$Granger)
}


##### Example

data<-  simul.farx(m=25,n=100,base=base.simul.far(24,5),
                   base.exo=base.simul.far(24,5),
                   d.a=matrix(c(0.5,0),nrow=1,ncol=2),
                   alpha.conj=matrix(c(0.2,0),nrow=1,ncol=2),
                   d.rho=diag(c(0.45,0.90,0.34,0.45)),
                   alpha=diag(c(0.5,0.23,0.018)),
                   d.rho.exo=diag(c(0.45,0.90,0.34,0.45)),
                   cst1=0.05)
y <- select.fdata(data,name="X")
x <- select.fdata(data,name="Z")
data<-as.fdata(list("X"=x[[1]],"Y"=y[[1]]))
#plot of the two whole series
multplot(data,xval=seq(from=0,to=100-0.04,by=0.04),whole=T,legend = T)
# test with 2 DFPCA
DFPCA.test(data$X,data$Y,num=25,ndx = 2, ndy = 2)
