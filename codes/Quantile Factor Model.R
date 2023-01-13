loss.fun <-function(tau, vec){
  return(mean(sapply(vec,function(x){return(ifelse(x>0, tau*x, (1-tau)*abs(x)))})))
}

tau.vec <- seq(from=0.5, to=0.95, by=0.05)

cdf.est <- function(x, pred.quant, tau.vec=tau.vec){
  # pred.quant : (T x S) x K, K : number of quantile levels. #tau.vec 별로 quantile 값
  # in this function, pred.quant is K-dim vector. 
  # Then this function would be applied to the ((ST)XK) object pred.quant. 
  if(x < pred.quant[1]){
    return(0)
  }else if(x > max(pred.quant)){
    return(1)
  }else{
    index <- rank(c(x, pred.quant), ties="min")[1]
    # if(index==50){
     return(index)
    # }else{
    #return(tau.vec[index-1]+(x-pred.quant[index-1])*(tau.vec[index]-tau.vec[index-1])/(pred.quant[index]-pred.quant[index-1]))
    # }
  }
}
library(matlib)
#######################################################################################
## when AB=0 : find estimators under the condition $\frac{1}{T}F_\tau^TF_\tau = I_{\hat{r}}$ and $\Lambda_\tau^T\Lambda_\tau$ is a diagonal matrix with non-increasing elements. 
##             which was suggested in Chen(2019) "Quantile Factor Models"
## when AB=1 : find estimators under the restriction according to bai(2013) PC2 or, Ando and Bai(2020) Quantile Co-Movement in FinancialMarkets \dots.
##             roughly, $\frac{1}{T}F_\tau^TF_\tau = I_{\hat{r}}$ and the upper part of Lamda is lower-triangular matrix.
#######################################################################################

do_qfm <- function(X,r.vec,AB=0,save=1,s.ind=0,fname=""){
  N = ncol(X)
  Te = nrow(X)
  temp <- svd(X)
  pred.quant.i = list()
  for(ind in 1:length(tau.vec)){
    print(paste("###########", ind, "th QFM estimation###########", sep=""))
    tau <- tau.vec[ind]
    # To minimize the target function of QFM, we use fully iterative approach. 
    # L.hat: current loading
    # L.hat.new : updated loading
    # F.hat : current factor
    # F.hat.new : updated factor
    # maximum iteration number: 200 
    # If spanned space by updated factors is near to that by previous ones enough, stop the iteration.
    #L.hat <- temp$loadings[,1:r]
    #F.hat <- temp$scores[,1:r]
    r = ifelse(length(r.vec)==1,r.vec,r.vec[ind])
    L.hat = sqrt(N)*(temp$v)[,1:r]
    
    L.hat.new <- matrix(0, nrow=N, ncol=r)
    F.hat.new <- matrix(0, nrow=Te, ncol=r)
    
    iter <- 1
    
    while(AB){
      for(i in 1:N){
        L.hat.new[i,] <- rq(X[,i] ~ F.hat-1, tau=tau)$coeff
      }
      for(t in 1:nrow(F.hat)){
        F.hat.new[t,] <- rq(X[t,] ~ L.hat.new-1, tau=tau)$coeff
      }
      q1 = qr(F.hat.new)
      q2 = qr(qr.R(q1)%*%t(L.hat.new))
      F.hat.new = sqrt(Te) * qr.Q(q1) %*% qr.Q(q2)
      L.hat.new = t(qr.R(q2))/sqrt(Te)
      if(iter>1 ){
        if(abs(loss.fun(tau,as.vector(X-F.hat.new %*% t(L.hat.new))) - loss.fun(tau,as.vector(X-F.hat %*% t(L.hat))))< 1e-10 ){
          print(loss.fun(tau,as.vector(X-F.hat.new %*% t(L.hat.new))) - loss.fun(tau,as.vector( X-F.hat %*% t(L.hat))))
          print(iter)
          break
        }
      }
      iter <- iter+1
      F.hat <- F.hat.new
      L.hat <- L.hat.new
      if(iter==200){
        print(paste("Max iteration at tau=", tau))
        break;
      }
    }
    while(!AB){
      for(t in 1:Te){
        F.hat.new[t,] <- rq(X[t,] ~ L.hat-1, tau=tau)$coeff
      }
      for(i in 1:N){
        L.hat.new[i,] <- rq(X[,i] ~ F.hat.new-1, tau=tau)$coeff
      }
      
      if(iter>1 ){
        if(abs(loss.fun(tau,as.vector(X-F.hat.new %*% t(L.hat.new))) - loss.fun(tau,as.vector(X-F.hat %*% t(L.hat))))< 1e-10 ){
          print(loss.fun(tau,as.vector(X-F.hat.new %*% t(L.hat.new))) - loss.fun(tau,as.vector( X-F.hat %*% t(L.hat))))
          print(iter)
          break
          #save_q1 = list(F = F.hat.new, L = L.hat.new)
          L.hat = svd(L.hat.new)
          q1 = F.hat.new %*% L.hat$v %*% diag(L.hat$d) %*% t( L.hat$v)
          q2 = svd(q1)
          L.hat.new = F.hat.new %*% L.hat$v %*% diag(1/abs(L.hat$d)) %*% t( L.hat$v) %*% q2$v * sqrt(N)
          F.hat.new = F.hat.new %*% t(q2$v) %*% L.hat$v %*% diag(abs(L.hat$d)) %*% t( L.hat$v) / sqrt(N)
          
          break
        }
      }
      iter <- iter+1
      F.hat <- F.hat.new
      L.hat <- L.hat.new
      if(iter==200){
        print(paste("Max iteration at tau=", tau))
        break;
      }
    }
    pred.quant.i [[ind]] = list(F=F.hat.new, L=L.hat.new)
  }
  if(save){
    filename.f <- paste("factor(r=",r,",index=",s.ind, ").RDS", sep="")
    saveRDS(pred.quant.i, file=ifelse(fname=="",filename.f,fname))
  }
  return(pred.quant.i)
}
#######################################################################################
##find_r is a function that estimates the number of factors according to Chen(2019) "Quantile Factor Models"
#######################################################################################

find_r<-function(X,r=80){
  N = ncol(X)
  Te = nrow(X)
  pred.quant.i = do_qfm(X,r,AB=0,save=1,0)
  r.vec = rep(0,length(tau.vec))
  for(ind in 1:length(tau.vec)){
    L.hat.new=pred.quant.i[[ind]]$L
    pnt = diag(t(L.hat.new)%*%L.hat.new )
    r.vec[ind]=sum(pnt >pnt[1]/N**(1/3))
    print(r.vec[ind])
  }
  return(r.vec)
}
r.vec = find_r(list.em$X,r=20)

#Example
library(quantreg)
r=3
pred.quant.i =do_qfm(list.em$X, r, AB=0, save=1,fname="QFM.RDS") 
