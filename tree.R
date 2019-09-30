# Tree growing algorithm for Boost-R.R

# ---------------------------------------------------------------
# For questions, please contact:
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# https://sites.google.com/site/liuxiaosite1/
# ---------------------------------------------------------------
# ---------------------------------------------------------------

tree = function(m.hat.tree, gamma, lambda, data.tree, D.max){
  
  output.list = output.tree = list()
  
  # ------------------------------------
  # root node
  # ------------------------------------
  D = 1
  
  # preparation for optim: 
  N.data = nrow(data.tree)
  T.max = ceiling( max(apply(data.tree, 1, max, na.rm=TRUE), na.rm=TRUE) )
  T.max = min(T.max, ncol(m.hat.tree))
  M.hat.node = m.hat.tree[,1:T.max]
  m =  array(0/0, dim=c(N.data, T.max))
  Tmp.mcf.list = lapply(1:N.data, mcf.group, data.tree)
  
  n.min = round(N.data/10) 
  n.min = 50 
  
  # get m_{i,j};  
  # This is the computation bottleneck. 
  for (i in 1:N.data){
    tmp = Tmp.mcf.list[[i]]
    if (nrow(tmp)>0){
      for (j in 1:ceiling(tmp$time[nrow(tmp)])){
        tmp1 = findInterval(j, tmp$time)
        m[i,j] = c(0,tmp$mcf)[tmp1+1]
      }
      #if ((j+1)<=t.max){m[i,(j+1):t.max] = 0/0}
    }else{
      m[i,] = 0/0
    }
  }
  
  # get g_{\cdot,j} and h_{\cdot,j}
  g.j = 2 * colSums(M.hat.node-m, na.rm=TRUE)
  case = which(is.na(g.j))
  h.j = rep(2 * N.data,T.max) 
  if (length(case)>0){h.j[case] = 0/0}
  
  para.opt = optim(par=c(0.01,log(1)), node,
                   g.j.in = g.j, h.j.in = h.j, t.max.in = T.max, gamma.in = gamma,
                   method = c("BFGS"),
                   control = list(trace=6), hessian = FALSE)
  
  obj.before = para.opt$value + lambda * D
  obj.before = round(obj.before,4)
  node.individual = para.opt$value
  node.individual.before =node.individual
  
  # record the tree structure here: 
  node.vector = c(0)  # length equals number of leafs D; 1 means do not further split this leaf
  node.opt = c(para.opt$value) # the optimum objective function on the leaf
  x.low = array(0, dim=c(n.x,D))
  x.high = array(1, dim=c(n.x,D))
  data.id.list = list() # a list of data IDs on the D terminal node
  data.id.list[[1]] = c(1:n.data)
  
  # ------------------------------------
  # Start the tree node splitting iterations
  # ------------------------------------
  flag = 0 # if flag=1, stop the tree
  node.vector.after = array()
  tree.level = 1
  importance = array(0, dim=c(n.x,1))
  
  while (flag==0){
    
    #print(paste("leaf:",D,sep=""))
    
    output.after.list = list() # store the results if this node i is splitted
    # get all terminal node; then, split all of them to see if there is a gain
    for (ii in 1:D){
      data.id = data.id.list[[ii]]
      if (length(data.id)<=n.min){node.vector[ii]=1}
      #output.after.list = list() # store the results if this node i is splitted
      
      if (node.vector[ii]==0){
        # find the optimum splitting value and splitting variable
        count = 1
        obj.gain =  array()
        obj.after.list = para.opt.left.list = para.opt.right.list = list()
        jj.start = 1
        
        search.range.list = list() # a list that contains the search range for individual features
        obj.gain.list = list()
        obj.gain.min = Inf
        for (j in 1:n.x){
          # parallel
          search.step = (x.high[j,ii]-x.low[j,ii])/20
          search.range = seq(x.low[j,ii]+search.step, x.high[j,ii]-search.step, search.step)
          #print(search.range)
          if (max(search.range)-min(search.range)<=0.02){search.range=mean(search.range)}
          search.range.list[[j]] = search.range
          
          dummy = foreach(jj = search.range) %dopar% {
            # get the left daughter
            case.left = which(X[data.id,j]<jj)
            data.left = data.tree[data.id[case.left], ]
            # get the right daughter
            data.right = data.tree[data.id[-case.left], ]
            
            para.opt.left.list = "oh"
            para.opt.right.list = "my"
            obj.after.list = "god"
            obj.gain = Inf
            output.foreach = list()
            
            if ( (nrow(data.left)>n.min) & (nrow(data.right)>n.min) ){
              # get the optimum objective on the left
              # preparation for optim: 
              N.data = nrow(data.left)
              T.max = ceiling( max(apply(data.left, 1, max, na.rm=TRUE), na.rm=TRUE) )
              T.max = min(T.max, ncol(m.hat.tree[data.id[case.left],]))
              M.hat.node = m.hat.tree[data.id[case.left],1:T.max]
              m =  array(0/0, dim=c(N.data, T.max))
              Tmp.mcf.list = lapply(1:N.data, mcf.group, data.left)
              
              # get m_{i,j};  
              # This is the computation bottleneck. 
              for (i in 1:N.data){
                #data.tmp = data.input[i, ]
                tmp = Tmp.mcf.list[[i]]
                if (nrow(tmp)>0){
                  for (j in 1:ceiling(tmp$time[nrow(tmp)])){
                    tmp1 = findInterval(j, tmp$time)
                    m[i,j] = c(0,tmp$mcf)[tmp1+1]
                  }
                  #if ((j+1)<=t.max){m[i,(j+1):t.max] = 0/0}
                }else{
                  m[i,] = 0/0
                }
              }
              
              # get g_{\cdot,j} and h_{\cdot,j}
              g.j = 2 * colSums(M.hat.node-m, na.rm=TRUE)
              case = which(is.na(g.j))
              h.j = rep(2 * N.data,T.max) 
              if (length(case)>0){h.j[case] = 0/0}
              
              para.opt.left.list  = optim(par=c(0.01,log(1)), node,
                                          g.j.in = g.j, h.j.in = h.j, t.max.in = T.max, gamma.in = gamma,
                                          method = c("BFGS"),
                                          control = list(trace=6), hessian = FALSE)
              
              # get the optimum objective on the right
              # preparation for optim: 
              N.data = nrow(data.right)
              T.max = ceiling( max(apply(data.right, 1, max, na.rm=TRUE), na.rm=TRUE) )
              T.max = min(T.max, ncol(m.hat.tree[data.id[-case.left],]))
              M.hat.node = m.hat.tree[data.id[-case.left],1:T.max]
              m =  array(0/0, dim=c(N.data, T.max))
              Tmp.mcf.list = lapply(1:N.data, mcf.group, data.right)
              
              # get m_{i,j};  
              # This is the computation bottleneck. 
              for (i in 1:N.data){
                #data.tmp = data.input[i, ]
                tmp = Tmp.mcf.list[[i]]
                if (nrow(tmp)>0){
                  for (j in 1:ceiling(tmp$time[nrow(tmp)])){
                    tmp1 = findInterval(j, tmp$time)
                    m[i,j] = c(0,tmp$mcf)[tmp1+1]
                  }
                  #if ((j+1)<=t.max){m[i,(j+1):t.max] = 0/0}
                }else{
                  m[i,] = 0/0
                }
              }
              
              # get g_{\cdot,j} and h_{\cdot,j}
              g.j = 2 * colSums(M.hat.node-m, na.rm=TRUE)
              case = which(is.na(g.j))
              h.j = rep(2 * N.data,T.max) 
              if (length(case)>0){h.j[case] = 0/0}
              
              para.opt.right.list  = optim(par=c(0.01,log(1)), node,
                                           g.j.in = g.j, h.j.in = h.j, t.max.in = T.max, gamma.in = gamma,
                                           method = c("BFGS"),
                                           control = list(trace=6), hessian = FALSE)
              
              
              # get the optimum objective after splitting the node
              obj.after.list= para.opt.left.list$value + para.opt.right.list$value 
              obj.after.list= round(obj.after.list,4)
              #obj.gain = (obj.after.list - obj.before + lambda) # This line seems to be incorrect. 
              obj.gain = (obj.after.list - node.individual[ii] + lambda) 
              #count = count + 1
            }
            
            output.foreach[[1]] = para.opt.left.list 
            output.foreach[[2]] = para.opt.right.list 
            output.foreach[[3]] = obj.after.list # This is the sum of loss after a node is splitted
            output.foreach[[4]] = obj.gain
            return(output.foreach)
          }
          
          jjj = 1
          tmp = array()
          for (jj in jj.start:(jj.start-1+length(search.range))){
            para.opt.left.list[[jj]] = dummy[[jjj]][[1]]
            para.opt.right.list[[jj]] = dummy[[jjj]][[2]]
            obj.after.list[[jj]] = dummy[[jjj]][[3]]
            obj.gain[jj] = dummy[[jjj]][[4]]
            tmp[jjj] = dummy[[jjj]][[4]] 
            jjj = jjj + 1
          }
          jj.start = jj+1
          obj.gain.list[[j]] = tmp
          
          if (min(obj.gain)<obj.gain.min){
            obj.gain.min = min(obj.gain)
            x.opt.tmp = j
          }
          #print(c(j,jj))
        }
        
        #print(obj.gain)
        if (min(obj.gain) <0){ # otherwise, do not split the node
          case = which(obj.gain == min(obj.gain))[1]
          obj.optim = obj.after.list[[case]]
          para.optim.left = para.opt.left.list[[case]]$par
          para.optim.right = para.opt.right.list[[case]]$par
          optim.left = para.opt.left.list[[case]]$value
          optim.right = para.opt.right.list[[case]]$value
          
          x.opt = x.opt.tmp
          
          value.opt = search.range.list[[x.opt]][which(obj.gain.list[[x.opt]]==min(obj.gain))[1]]
          #print(which(obj.gain.list[[x.opt]]==min(obj.gain)))
          
          case.left.opt = data.id[which(X[data.id,x.opt]<value.opt)]
          output.after.list[[ii]] = c(obj.optim, x.opt, value.opt, 
                                      para.optim.left, para.optim.right, optim.left, optim.right, case.left.opt)
          node.vector.after[ii] = 0
          if (length(search.range.list[[x.opt]])==1){node.vector.after[ii] = 1}
          
          # update importance
          importance[x.opt,1] = importance[x.opt,1] + min(obj.gain)-lambda
        }else{
          # set node.vector to 1 for this node
          node.vector.after[ii] = 1
        }
        
        
      }else{
        node.vector.after[ii] = node.vector[ii]
      } # if ... else ...
      
    } # for loop
    
    # update node termination
    node.vector = node.vector.after
    
    # update the tree structure
    x.low.list = x.high.list = para.opt.left.list = para.opt.right.list = para.opt.list = list()
    #node.opt.after = obj.before
    node.opt.after = 0
    node.individual = array()
    count = 1
    count2 = 1
    for (i in 1:D){
      if (node.vector[i]==0){
        para.opt.left.list[[i]] = matrix( output.after.list[[i]][4:5], ncol=1 )
        para.opt.right.list[[i]] = matrix( output.after.list[[i]][6:7], ncol=1 )
        para.opt.list[[i]] = cbind(para.opt.left.list[[i]], para.opt.right.list[[i]])
        node.individual[count2:(count2+1)] = output.after.list[[i]][8:9]
        count2 = count2 + 2
        
        x.opt = output.after.list[[i]][2]
        x.low.left = x.low[,i]
        #x.low.left[x.opt,1] = x.low[x.opt,i]
        x.high.left = x.high[,i]
        x.high.left[x.opt] = output.after.list[[i]][3]
        
        x.low.right = x.low[,i]
        x.low.right[x.opt] = output.after.list[[i]][3]
        #x.low.left[x.opt,1] = x.low[x.opt,i]
        x.high.right = x.high[,i]
        
        x.low.list[[i]] = cbind(x.low.left, x.low.right)
        x.high.list[[i]] = cbind(x.high.left, x.high.right)
        
        # update the optimum node.opt
        node.opt.after = node.opt.after + output.after.list[[i]][1] + 2*lambda
        
        # update data list
        data.id.left.after =  output.after.list[[i]][10:length(output.after.list[[i]])]
        tmp = match(data.id.list[[i]],output.after.list[[i]][10:length(output.after.list[[i]])])
        data.id.right.after =  data.id.list[[i]][which(is.na(tmp))]
        
        data.id.list[[count]] = data.id.left.after
        count = count + 1
        data.id.list[[count]] = data.id.right.after
        count = count + 1
      }else{
        if (D==1){
          para.opt.list[[i]] = matrix( para.opt$par, ncol=1 )
        }else{
          para.opt.list[[i]] = matrix( para.opt[,i], ncol=1 )
        }
        node.individual[count2] = node.individual.before[i]
        count2 = count2+1
        node.opt.after = node.opt.after + node.individual.before[i] + lambda
        
        x.low.list[[i]] = matrix(x.low[,i],ncol=1)
        colnames(x.low.list[[i]]) = colnames( x.low)[i]
        x.high.list[[i]] = matrix(x.high[,i],ncol=1)
        colnames(x.high.list[[i]]) = colnames( x.high)[i]
        data.id.list[[count]] = data.id.list[[i]]
        count = count +1
      }
    }
    
    x.low = do.call(cbind, x.low.list)
    x.high = do.call(cbind, x.high.list)
    para.opt = do.call(cbind, para.opt.list)
    node.opt = node.opt.after
    obj.before = node.opt # This is the total loss including penalty
    D.old = D
    D = ncol(x.low)
    node.individual.before = node.individual
    
    if (sum(node.vector)==length(node.vector)){flag=1}
    #print(D)
    #print(node.vector)
    #print(node.opt)
    #print(x.low)
    #print(x.high)
    
    # get the new node.vector
    tmp.list = list()
    for (i in 1:D.old){
      if (node.vector.after[i]==0){
        tmp.list[[i]] = c(0,0)
      }else{
        tmp.list[[i]] = 1
      }
    }
    node.vector = do.call(c, tmp.list)
    
    
    output.list[[1]] = x.low
    output.list[[2]] = x.high
    output.list[[3]] = para.opt
    output.list[[4]] = node.opt
    output.list[[5]] = data.id.list
    output.list[[6]] = importance
    #output.list[[6]] = data.id.list
    
    output.tree[[tree.level]] = output.list
    tree.level = tree.level + 1
    
    #print(node.opt)
    
    if (D>=D.max){flag=1}
  } # while loop
  
  
  return(output.tree)
  
}

