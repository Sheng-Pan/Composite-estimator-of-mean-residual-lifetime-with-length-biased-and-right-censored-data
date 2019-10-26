S_a <- function(a, v, Id, t)   # V = Y - A   t <- Y
{
  Tt <- cbind(time = t, address = c(1: length(t)))
  Aa <- cbind(time = a, address = rep(0, length(a)))
  Vv <- cbind(time = v, address = rep(0, length(v)))
  u <- rbind(Tt, Aa, Vv)
  u <- u[order(u[,1]), ]
  Q <- colSums(outer(a, u[,1], function(x,y) as.numeric(x <= y)) + Id * outer(v, u[,1], function(x,y) as.numeric(x <= y)))
  dQ <- Q - c(0, Q[-length(Q)])
  K <- colSums(outer(a, u[,1], function(x,y) as.numeric(x >= y)) + outer(v, u[,1], function(x,y) as.numeric(x >= y)))
  Sa_avt <- cumprod(1 - divide(dQ, K))        #value at points a_i, v_i and t_i
  u <- cbind(u, Sa_avt)[order(u[,2]),]
  Sa_hat <- u[which(u[,2] > 0),3]            # Sa_hat is not ordered from 0 to 1, but with the same order of the input Y
  return(Sa_hat)
}

S.HQ <- function(a, y, Id, t, is.order)    ####   is.order = 1 represent that Fhat is order from 0 to 1
{
  n <- length(a)
  Tt <- cbind(time = t, address = c(1: length(t)))
  Yy <- cbind(time = y, address = rep(0, length(y)))
  u <- rbind(Tt, Yy)
  u <- u[order(u[,1]), ]
  R <- colSums(outer(y, u[,1], function(x,y) as.numeric(x >= y ))) - n * S_a(a, y-a, Id, u[,1])
  N <- colSums(Id * outer(y, u[,1], function(x,y) as.numeric(x <= y )))
  dN <- N - c(0, N[-length(N)])
  F_ty <- 1 - cumprod(1 - divide(dN, R))
  if(is.order == 1)
  {
    u <- cbind(u, F_ty)
    F_t <- u[which(u[,2] > 0),3]
  }
  else
  {
    u <- cbind(u, F_ty)[order(u[,2]),]
    F_t <- u[which(u[,2] > 0),3]
  }
  return(1-F_t)
}


S.new <- function(a, y, Id, t, is.order)           ####   is.order = 1 represent Fhat is order from 0 to 1
{
  Tt <- cbind(time = t, address = c(1: length(t)))
  Yy <- cbind(time = y, address = rep(0, length(y)))
  u <- rbind(Tt, Yy)
  u <- u[order(u[,1]), ]
  N <- colSums(Id * outer(y, u[,1], function(x,y) as.numeric(x <= y )))
  dN <- N - c(0, N[-length(N)])
  K1 <- outer(a, u[,1], function(x,y) as.numeric(x <= y )) * outer(y, u[,1], function(x,y) as.numeric(x >= y ))
  K2 <- outer((y - a), u[,1], function(x,y) as.numeric(x <= y )) * outer(y, u[,1], function(x,y) as.numeric(x >= y ))
  K <- colSums(K1 + Id * K2) / 2 
  F_ty <- 1 - cumprod(1 - divide(dN, K))
  if(is.order == 1)
  {
    u <- cbind(u, F_ty)
    F_t <- u[which(u[,2] > 0),3]         ###  reserve Fhat at t and get rid of those Yy 
  }
  else
  {
    u <- cbind(u, F_ty)[order(u[,2]),]
    F_t <- u[which(u[,2] > 0),3]
  }
  return(1-F_t)
}
S.v <- function(y, Id, t, is.order)           ####   is.order = 1 represent Fhat is order from 0 to 1
{
  Tt <- cbind(time = t, address = c(1: length(t)))
  Yy <- cbind(time = y, address = rep(0, length(y)))
  u <- rbind(Tt, Yy)
  u <- u[order(u[,1]), ]
  N <- colSums(Id * outer(y, u[,1], function(x,y) as.numeric(x <= y )))
  dN <- N - c(0, N[-length(N)])
  K <- colSums(outer(y, u[,1], function(x,y) as.numeric(x >= y )))
  F_ty <- 1 - cumprod(1 - divide(dN, K))
  if(is.order == 1)
  {
    u <- cbind(u, F_ty)
    F_t <- u[which(u[,2] > 0),3]         ###  reserve Fhat at t and get rid of those Yy 
  }
  else
  {
    u <- cbind(u, F_ty)[order(u[,2]),]
    F_t <- u[which(u[,2] > 0),3]
  }
  return(1-F_t)
}
###C的生存函数K-M estimation
S.c <- function(y, Id, t, is.order)           ####   is.order = 1 represent Fhat is order from 0 to 1
{
  Tt <- cbind(time = t, address = c(1: length(t)))
  Yy <- cbind(time = y, address = rep(0, length(y)))
  u <- rbind(Tt, Yy)
  u <- u[order(u[,1]), ]
  N <- colSums((1-Id) * outer(y, u[,1], function(x,y) as.numeric(x <= y )))
  dN <- N - c(0, N[-length(N)])
  K <- colSums(outer(y, u[,1], function(x,y) as.numeric(x >= y )))
  F_ty <- 1 - cumprod(1 - divide(dN, K))
  if(is.order == 1)
  {
    u <- cbind(u, F_ty)
    F_t <- u[which(u[,2] > 0),3]         ###  reserve Fhat at t and get rid of those Yy 
  }
  else
  {
    u <- cbind(u, F_ty)[order(u[,2]),]
    F_t <- u[which(u[,2] > 0),3]
  }
  return(1-F_t)
}
#S.mle(bootA,bootY,bootdelta,xx,1) S.mle(A,Y,delta,xx,1) 
S.mle<- function(a,y, id,x,is.order)
{ pnew <- Vardi(20, y, id)
t <- sort(unique(y))
n_F <- length(t)
nk <- outer(t, t, function(x,y) as.numeric(x <= y ))
col_nk <- apply(nk, 2, function(x) max(which(x == 1)))
ptproduct <- 1 / t * pnew;
oneupper <- matrix(1, n_F, n_F)
oneupper[lower.tri(oneupper)] <- 0
sum <- t(ptproduct) %*% oneupper
hatF <- sum[col_nk] / sum[n_F]                     # estimation of unbiased distribution \hat{F}!!!
sm<-1-hatF
# sm[which(sm==min(sm))]=sm[which(sm==min(sm))]+1e-20
v1<-outer(t, x, function(x,y) as.numeric(x <= y))#找出x的每个元素在Y中的位置
v0<-rep(0,length(x))
for(k in 1:length(x))
{v2<-sm*v1[,k]
if(min(v1[,k])==1){
  v0[k]<-0
}else{
  v0[k]<-min(v2[which(v2>0)],1)
}
}
return(v0)
}