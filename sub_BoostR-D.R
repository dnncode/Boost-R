# The Boost-R-D main algorithm; called by Boost-R-D.R

# ---------------------------------------------------------------
# For questions, please contact:
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# https://sites.google.com/site/liuxiaosite1/
# ---------------------------------------------------------------
# ---------------------------------------------------------------

BoostR2 = function(data, X, Z, K.value, lambda.value, gamma.value, D.max, u.value, v.value){

  
  K = K.value # number of trees
  lambda = lambda.value
  gamma = gamma.value
  
  source("sub_node_2.R")
  source("mcf_function.R")
  source("tree_2.R")
  clusterEvalQ(cl, source("sub_node_2.R"))
  clusterEvalQ(cl, source("mcf_function.R"))
  clusterEvalQ(cl, source("tree_2.R"))
  clusterEvalQ(cl, library(foreach))
  clusterEvalQ(cl, library(optimParallel))
  clusterEvalQ(cl, library("splines2"))
  clusterEvalQ(cl, library("splines"))
  # --------------------------------
  # Boosting iterations
  # --------------------------------
  t.max = ceiling( max(apply(data, 1, max, na.rm=TRUE), na.rm=TRUE) )
  n.data = nrow(data)
  n.x = ncol(X)
  
  # Initial estimate
  #para0 = c(a,b)
  model.initial = function(t){
    b = 0
    a = 0.001
    H = a * t^(b)
    return(H)
  }
  m.hat = model.initial(c(1:t.max))
  m.hat = matrix(rep(m.hat,each=n.data),nrow=n.data)
  
  m.hat.list = list()
  m.hat.list[[1]] = m.hat
  
  tree.seq = list() # record each tree
  for (k in 1:(K+1)){
    print(k)
    # prepare the current estimates (k-1):
    m.hat.in = m.hat.list[[k]]
    tree.out = tree(m.hat.tree = m.hat.in, gamma, lambda, data.tree = data, D.max=D.max,z.tree=Z, u=u.value, v=v.value)
    tree.seq[[k]] = tree.out
    tree.level = length(tree.out)
    x.low = tree.out[[tree.level]][[1]] 
    x.high = tree.out[[tree.level]][[2]]
    print(ncol(x.low))
    #sum(tree.out[[3]][[3]][,1]^2)
    
    m.hat.new = array(0/0, dim=c(n.data, ncol(m.hat)))
    for (i in 1:n.data){
      # determine the leaf
      case = c(1:ncol(x.low))
      for (j in 1:n.x){
        case1 = which( (X[i,j] - x.low[j,]) * (X[i,j] - x.high[j,])<=0)
        case = intersect(case, case1)
      }  
      case = case[1]
      
      # get m.hat.new
      v = v.value
      u = u.value
      para.opt = tree.out[[tree.level]][[3]][,case]
      para.opt = matrix( para.opt, nrow=u+v)  
      
      f = array(0,dim=c(nrow(Z[[i]]),1))
      n.z = ncol(Z[[1]])
      for (l in 1:n.z){
        B <- bs(x=Z[[i]][,l], degree=v, knots=u, intercept = TRUE, Boundary.knots=c(0, 1))
        
        for (j in 1:ncol(B)){
          f = f + para.opt[j,l] * cumsum(B[,j])
        }
      }
      
      m.hat.new[i,] =  c(f, rep(0/0,t.max-length(f)))
    }
    
    
    
    # update the estimates after tree k has been grown:
    m.hat.list[[k+1]] = m.hat.list[[k]] + learning.rate * m.hat.new
    
  
  }
  
  output = list()
  output[[1]] = tree.seq
  output[[2]] = m.hat.list
  return(output)
}