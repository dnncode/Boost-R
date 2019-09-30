# Sub-routines for visualizations; required by Boost-R-D.R
# ---------------------------------------------------------------
# For questions, please contact:
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# https://sites.google.com/site/liuxiaosite1/
# ---------------------------------------------------------------
# ---------------------------------------------------------------
library(spectral)

# binary paritions of individual trees (dim(X)=2 only)
Plot_Partition = function(i, BoostR.obj){
  graphics.off()
  tree.seq = BoostR.obj[[1]]
  
  n.tree = length(tree.seq)
  tree.i = tree.seq[[1]]
  x.low = tree.i[[length( tree.i)]][[1]]
  n.x = nrow(x.low)
  
  if (i>n.tree){
    print("the specified i exceeds the number of trees in the ensemble")
  }else if (n.x != 2){
    print("this 2D visualization is only available when there are only two static features")
  }else{
    par(mfrow=c(1,2))
    par(mar=c(3.2,3.2,2,1))
    
    tree.i = tree.seq[[i]]
    
    # Get the binary partition
    x.low = tree.i[[length( tree.i)]][[1]]
    x.high = tree.i[[length( tree.i)]][[2]]
    
    plot(0,0,xlim=c(-0.01,1.01),ylim=c(-0.01,1.01),pch=20,cex=0.001,col="white",
         xlab="",ylab="")
    mtext(text = expression(x[1]),
          side = 1,line = 2)
    mtext(text = expression(x[2]),
          side = 2,line = 2)
    title(paste("tree:",i,sep=""))
    
    points(X,pch=3,col="cornsilk2",cex=1)
    n.leaf = ncol(x.low)
    for (j in 1:n.leaf){
      rect(xleft=x.low[1,j], ybottom=x.low[2,j], xright=x.high[1,j], ytop=x.high[2,j], 
           density = NULL, angle = 45,
           col = "NA", border = "blue", lty = 1, lwd = 1)
      text(x=(x.low[1,j]+x.high[1,j])/2,  
           y=(x.low[2,j]+x.high[2,j])/2,
           labels=j)
    }
    
    # Get the estimated 
    para = ( tree.i[[length( tree.i)]][[3]] )
    para[2,] = exp( para[2,])
    x = c(0:max(data))
    curve = array(0/0,dim=c(length(x),n.leaf))
    for (j in 1:n.leaf){
      curve[,j] = para[1,j]*x^para[2,j]
    }
    
    j = 1
    plot(x, curve[,j],
         type="l",xlab="",ylab="",
         ylim=c(-0.5,max(curve)*1.01), col="darkgreen")
    mtext(text = "time",
          side = 1,line = 2)
    #mtext(text = paste("contribution from tree", i, sep=" "),
    #side = 2,line = 2)
    mtext(text = "contribution from leafs",
          side = 2,line = 2)
    text(x=max(x)-(n.leaf+1-j)*2,  
         y=max(curve[,j]),
         labels=j)
    if (n.leaf>1){
      for (j in 2:n.leaf){
        lines(x, curve[,j], col="darkgreen")
        text(x=max(x)-(n.leaf+1-j)*2,  
             y=max(curve[,j]),
             labels=j)
      }
    }
  }
  
}

# number of leafs for all trees in the ensemble 
Plot_Leaf = function(BoostR.obj){
  graphics.off()
  tree.seq = BoostR.obj[[1]]
  
  n.tree = length(tree.seq)
  n.leaf = array(0,dim=c(n.tree-1,1))
  for (i in 1:(n.tree-1)){
    tree.i = tree.seq[[i]]
    x.low = tree.i[[length( tree.i)]][[1]]
    n.leaf[i] = ncol(x.low)
    par(mar=c(4,4,1,1))
    plot(c(1:(n.tree-1)),n.leaf,type="h",col="darkgreen",xlab="trees",ylab="# of leafs")
  }
}

# estimated cumulative intensity for selected individuals; i specifies the individual ID
Plot_Individual = function(i, BoostR.obj){
  graphics.off()
  tree.seq = BoostR.obj[[1]]
  m.hat.list = BoostR.obj[[2]]
  n.tree = length(tree.seq)
  
  par(mar=c(4,4,1,1))
  
  case = which(data[i,]>0)-1
  if (length(case)<2){
    
  }else{
    fail.times =c(0, as.numeric( data[i,case+1] ))
    
    y.max = length(fail.times)*1.01
    for (k in 2:n.tree){
      m.hat.now = m.hat.list[[k]]
      y.max = max(y.max, m.hat.now[i,], na.rm=TRUE)
    }
    
    plot( fail.times, c(0:(length(fail.times)-1)) ,type="p", 
          ylim=c(0,y.max),col="blue",pch=3,ylab="cumulative events",xlab="time")
    m.hat.now = m.hat.list[[1]]
    tmp = m.hat.now
    lines(m.hat.now[i,],col="grey",lty=2)
    text(x=max(fail.times,na.rm=TRUE)-10,  
         y=max(m.hat.now[i,]),pos=3,
         labels="initial")
    
    for (k in 2:(n.tree-1)){
      m.hat.now = m.hat.list[[k]]
      lines(m.hat.now[i,],col="darkgreen")
      #text(x=max(fail.times,na.rm=TRUE)-(n.tree+1-k)*4,  
      #y=max(m.hat.now[i,]),pos=1,
      #labels=k-1)
    }
    
    k = n.tree
    m.hat.now = m.hat.list[[k]]
    lines(m.hat.now[i,],col="red",lwd=3)
    
    dist = max(m.hat.now[i,])-max(tmp[i,])
    arrows(x0=max(fail.times,na.rm=TRUE)-10, 
           y0= max(tmp[i,])+dist*0.15, 
           x1 =max(fail.times,na.rm=TRUE)-10, 
           y1 =max(tmp[i,])+dist*0.73,
           length = 0.15, angle = 30, col="blue",lwd=2, code=2)
    
    title(paste("individual",i,sep=" "))
    
    legend("topleft",legend=c("final ensemble","intermidate ensemble","initial"),
           lty=c(1,1,2),col=c("red","darkgreen","grey"),bty="n")
  }
}

# variable important; standardize=TRUE generates the standardized variable importance
Plot_Imp = function(BoostR.obj, standardize){
  graphics.off()
  tree.seq = BoostR.obj[[1]]
  n.tree = length(tree.seq)
  tree.i = tree.seq[[1]]
  x.low = tree.i[[length( tree.i)]][[1]]
  n.x = nrow(x.low)
  
  imp = array(0, dim=c(n.x,1))
  for (i in 1:n.tree){
    imp = imp + tree.seq[[i]][[length(tree.seq[[i]])]][[6]]
  }
  if (standardize){
    imp.std = -(imp - max(imp))/(max(imp)-min(imp))
  }else{
    imp.std = imp
  }
  par(mar=c(4,4,1,1))
  barplot(height=as.numeric(imp.std), names.arg=c(1:n.x),
          ylab="feature importance",
          xlab="features")
}

# interaction between splines coefficients and static features (dim(X)=2 only)
Plot_Interaction(BoostR.obj){
  
  graphics.off()
  tree.seq = BoostR.obj[[1]]
  n.tree = length(tree.seq)
  tree.i = tree.seq[[1]]
  x.low = tree.i[[length( tree.i)]][[1]]
  n.x = nrow(x.low)
  
  if (n.x >2){
    print("this 2D visualization is only available when there are only two static features")
  }else{
    X_1 = X_2 = seq(0,1,0.01)
    output1 = output2 = output3 = output4 = output5 = array(0/0, dim=c(length(X_1),length(X_2)))
    i = 1
    for (x_1 in X_1){
      j = 0
      for (x_2 in X_2){
        j = j +1
        
        para.opt = array(0,dim=c(5,n.tree))
        for (k in 1:n.tree){
          tree.i = tree.seq[[k]]
          tree.level = length(tree.i)
          x.low = tree.i[[tree.level]][[1]] 
          x.high = tree.i[[tree.level]][[2]]
          
          case1 = which( (x_1 - x.low[1,]) * (x_1 - x.high[1,])<=0)
          case2 = which( (x_2 - x.low[2,]) * (x_2 - x.high[2,])<=0)
          case = Reduce(intersect, list(case1,case2))[1]
          
          para.opt[,k] = tree.i[[tree.level]][[3]][,case]
        }
        para = matrix( rowSums(para.opt), ncol=1)
        output1[i,j] = para[1]
        output2[i,j] = para[2]
        output3[i,j] = para[3]
        output4[i,j] = para[4]
        output5[i,j] = para[5]
      }
      i = i + 1
    }
    
    MAX = max(output1, output2, output3, output4, output5)
    MIN = min(output1, output2, output3, output4, output5)
    
    par(mfrow=c(2,3))
    par(las=TRUE)
    par(mar=c(4,4,2,1))
    rasterImage2(x = X_1,
                 y = X_2,
                 z = output1,xlab=expression(x[1]),ylab=expression(x[2]),
                 main="(a)",zlim=c(MIN,MAX))
    rasterImage2(x = X_1,
                 y = X_2,
                 z = output2,xlab=expression(x[1]),ylab=expression(x[2]),
                 main="(b)",zlim=c(MIN,MAX))
    rasterImage2(x = X_1,
                 y = X_2,
                 z = output3,xlab=expression(x[1]),ylab=expression(x[2]),
                 main="(c)",zlim=c(MIN,MAX))
    rasterImage2(x = X_1,
                 y = X_2,
                 z = output4,xlab=expression(x[1]),ylab=expression(x[2]),
                 main="(d)",zlim=c(MIN,MAX))
    rasterImage2(x = X_1,
                 y = X_2,
                 z = output5,xlab=expression(x[1]),ylab=expression(x[2]),
                 main="(e)",zlim=c(MIN,MAX))
  }
  
  
}








