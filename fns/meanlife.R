COMP_simpsons_rule<- function(f,a,b,n)
{
  h=(b-a)/n; xi=seq(a, b, by = h);xi=t(xi)
  aa = xi[seq(3, length(xi)-2, by = 2)];
  bb = xi[seq(2, length(xi), by = 2)];
  I= h/3*(f(xi[1])+2*sum(f(aa))+4*sum(f(bb))+f(xi[length(xi)]));
  return(I)
}
##mean residual lifetime estimation function,
##Arguments 
## s methods choices,  'new', 'mle', 'HQ'
## t time
## a left truncation time
## y survival time
## id 1
meanlife <-function(s,t,a,y,id)
{ 
  g <- get(paste0('S.',s))
  n0 <- 300
  h1 <- rep(0,length(t))
  f1 <- function(x){return(g(a,y,id,x,1))}
  for(i in 1:length(t))
  {if(g(a,y,id,t[i],1)==0){
    h1[i]<-0}else{
      h1[i]<-COMP_simpsons_rule(f1,t[i],6,n0)/g(a,y,id,t[i],1) 
    }
  }
  return(h1)
}
