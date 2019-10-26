generate_data <- function(n,thresh)
 { while( sum(is.select) < n)
  {
    W0 <- runif(52500, min = 0, max = 100)
    T0 <- rweibull(52500, shape = 2, scale = 2)
    is.select <- as.numeric(W0 + T0 >= thresh)
  }
W <- W0[which(is.select == 1)[1 : n]]
T <- T0[which(is.select == 1)[1 : n]]
A <- thresh - W
C <- runif(n, min = 0, max = 5.5)     #3.8 -- CR:30%
Y <- pmin(T, A + C)
delta <- as.numeric(T <= (A + C))
# # CR[i] <- 1 - sum(delta) / n
LBdata <- cbind(Y, delta, Y - A, A, 1 - delta,  C * (1-delta) , T * delta - A * delta)
}