# Sub-routines for Boost-R-D.R
# ---------------------------------------------------------------
# For questions, please contact:
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# https://sites.google.com/site/liuxiaosite1/
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# intial set up for cumulative intensity function
model=function(para,t){
  b = para[2]
  a = para[1]
  H = a * t^(b)
  return(H)
}

# generate mcf for a group of individuals
mcf.group = function(i, data.mcf.group){
  output = mcf(data.mcf.group[i,])
  return(output)
}


# Objective on a single node
# Input: 
# 1). all data on this node
# 2). the current model on this node 
# Output:
# 1). the objective function for any given parameters
node.2 = function(para.input, t.max.in, g.in, h.in, z.in, gamma.in,u.in, v.in){
  
  # v: order 3
  v = v.in
  u = u.in
  
  n.Id = nrow(g.in)
  n.Id = which( rowSums(g.in, na.rm=TRUE) != 0 )
  n.z = ncol(z.in[[1]])
  
  para.input.2 = matrix( para.input, nrow=u+v) 
  term = array(0,dim=c(nrow(g.in),1))
  for (i in n.Id){
    g.i = g.in[i,]
    h.i = h.in[i,]
    g.max = max( which(!is.na(g.i)) )
    
    # calculate f(t)
    f = array(0,dim=c(nrow(z.in[[i]]),1))
    for (l in 1:n.z){
      B <- bs(x=z.in[[i]][,l], degree=v, knots=u, intercept = TRUE, Boundary.knots=c(0, 1))
      
      for (j in 1:ncol(B)){
        f = f + para.input.2[j,l] * cumsum(B[,j])
      }
    }
    
    tmp = sum(g.i[1:g.max]*f[1:g.max] + 0.5*h.i[1:g.max]*f[1:g.max]^2)
    
    term[i,1] = tmp
    
  }
  
  obj = sum( term )
  
  # add the penalty term here:
  for (l in 1:n.z){
    obj = obj + 1/2 * gamma.in * sqrt(sum(para.input.2[,l]^2))
  }
  
  return(obj)
}


