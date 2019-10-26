Vardi <- function(times, y, Id)      
{
  sample_n <- length(y)
  YD <- y[Id == 1]
  YC <- y[Id == 0]
  if ( length(unique(y)) == length(y))
  {
    t <- sort(y)
    n_t <- length(t)
    xi <- colSums(outer(YD, t, function(x,y) as.numeric(x == y)))
    zeta <- colSums(outer(YC, t, function(x,y) as.numeric(x == y)))
    
    pold <- rep(1/n_t, n_t)                   #initial value of the iteration
    one <- matrix(1, n_t, n_t)
    one[upper.tri(one)] <- 0
    onelow <- one                        ###  lower matrix
    
    one <- matrix(1, n_t, n_t)
    one[lower.tri(one)] <- 0
    oneup <- one                         ###  upper matrix
    for( j in 1 : times)                 ###  iteration
    {
      pt <- pold / t
      ptsum <- t(pt) %*% onelow
      zetapro <- zeta * t(1 / ptsum)
      s <- t(zetapro) %*% oneup
      pnew <-  1/sample_n * (xi + pt * t(s))   #This is the estimation of length-bias distribution \hat{G}!!!
      pold <- pnew                                       
    }
    return(pnew)
  }
  else
  {
    warning('ties in the data');
    t <- sort(unique(y))
    n_t <- length(t)
    xi <- colSums(outer(YD, t, function(x,y) as.numeric(x == y)))
    zeta <- colSums(outer(YC, t, function(x,y) as.numeric(x == y)))
    #zeta <- c(rep(0, length(YD1)),rep(1, length(YC1)))    LiuPeng's code
    #xi <- c(rep(1, length(YD1)),rep(0, length(YC1)))
    
    pold <- rep(1/n_t, n_t)                   #initial value of the iteration
    one <- matrix(1, n_t, n_t)
    one[upper.tri(one)] <- 0
    onelow <- one                        ###  lower matrix
    
    one <- matrix(1, n_t, n_t)
    one[lower.tri(one)] <- 0
    oneup <- one                         ###  upper matrix
    for( j in 1 : times)                 ###  iteration
    {
      pt <- pold / t
      ptsum <- t(pt) %*% onelow
      zetapro <- zeta * t(1 / ptsum)
      s <- t(zetapro) %*% oneup
      pnew <-  1/n_t * (xi + pt * t(s))   #This is the estimation of length-bias distribution \hat{G}!!!
      pold <- pnew
      return(pnew)
    }
  }
}
