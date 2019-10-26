resamp <- function(arg1,arg2,arg3)
{
  nsamps <- arg1    #flag --- nsamps not set
  noresamp <- arg3   #flag --- 0 represents sample without replacement
  data <- arg2
  RowCol <- dim(as.matrix(data))
  
  if(noresamp == 0)
  {
    indstotal = sample.int(RowCol[1], RowCol[1])
    if (nsamps > RowCol[1])
    {
      warning('resamp: cannot sample more than length(data) without replacement.')
      warning('using replacement sampling.')
      inds = sample.int(RowCol[1], nsamps, replace = TRUE)
    }
    else
      inds = indstotal[1:nsamps]
  }
  else
  {
    inds = sample.int(RowCol[1], nsamps, replace = TRUE)
  }
  
  if(RowCol[1] == 1 || RowCol[2] == 1)
  {
    res = data[inds]
  }
  else
  {
    res = data[inds,]
  }
  return(res)
}