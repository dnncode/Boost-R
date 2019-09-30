# Sub-routines for Boost-R.R
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
node = function(para.input, t.max.in, g.j.in, h.j.in, gamma.in){
  
  para.input.2 = c(para.input[1],exp(para.input[2])) 
  omega = model(para.input.2, c(1:t.max.in))
  obj = sum( omega * g.j.in + 1/2 * (h.j.in + gamma.in) * omega^2 )
  return(obj)
}




