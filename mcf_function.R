# Compute the MCF for a given individual
# ---------------------------------------------------------------
# For questions, please contact:
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# https://sites.google.com/site/liuxiaosite1/
# ---------------------------------------------------------------
# ---------------------------------------------------------------

mcf = function(Data){
  
  Data = as.matrix(Data)
  # number of systems
  n = nrow(Data)
  # extract censoring times
  t.c = apply(Data, 1, max, na.rm=TRUE)
  # set the censoring time to -10
  tmp = lapply(1:n, mcf.sub.censor, T.c=t.c, Data.tmp=Data)
  Data1 = do.call(rbind, tmp) # the same as Data but with censoring time being set to -10
  
  # sort data
  tmp = as.vector(Data1)
  t.sort = sort(tmp)
  t.sort = t.sort[t.sort>0]
 
  # delta matrix
  tmp = lapply(1:n, mcf.sub.delta, T.c=t.c, T.sort=t.sort)
  delta = do.call(cbind,tmp)
  
  # d matrix
  tmp = lapply(1:n, mcf.sub.d, T.c=t.c, T.sort=t.sort, Data1.tmp=Data1)
  d = do.call(cbind,tmp)
  
  # complete the MCF table
  delta.dot = rowSums(delta)
  d.dot = rowSums(d)
  d.bar = d.dot/delta.dot
  mu = cumsum(d.bar)
  
  # compute the variance
  tmp = lapply(1:n, mcf.sub.var, Delta.dot=delta.dot,D.bar=d.bar,Delta=delta,D=d)
  tmp = do.call(cbind,tmp)
  v = rowSums(tmp)
  
  # output table
  output = cbind(t.sort, delta, d, delta.dot, d.dot, d.bar, mu, v)
  output = data.frame(output)
  colnames(output) = c("time",
                      paste("delta",c(1:n),sep=""),
                      paste("d",c(1:n),sep=""),
                      "delta.dot","d.dot","d.bar",
                      "mcf","variance")
  return(output)
}

mcf.sub.censor = function(i, T.c, Data.tmp){ # used to remove the censoring time
  case = which(Data.tmp[i,]==T.c[i])
  output = Data.tmp[i,]
  output[case] = -10
  output = matrix(output,nrow=1)
  return(output)
}
mcf.sub.delta = function(i, T.c, T.sort){  # used to construct the delta function
  case = which(T.sort<=T.c[i])
  output = array(0, dim=c(length(T.sort),1))
  output[case] = 1
  return(output) 
}
mcf.sub.d = function(i, T.c, T.sort, Data1.tmp){  # used to construct the d function
  case = which(T.sort %in% Data1.tmp[i,])
  output = array(0, dim=c(length(T.sort),1))
  output[case] = 1
  return(output) 
}
mcf.sub.var = function(i, Delta.dot, D.bar, Delta, D){ # used to compute the variance
  tmp = Delta[,i]/Delta.dot*(D[,i]-D.bar)
  tmp = cumsum(tmp)
  tmp = tmp^2
  return(tmp)
}

# ------------------------------------