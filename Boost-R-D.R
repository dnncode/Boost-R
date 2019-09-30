# ---------------------------------------------------------------
# ---------------------------------------------------------------
# An example of running Boost-R algorithm

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Data
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# DataA is the dataset A in the paper
# DataA is available on GitHub

# CaseStudy is the dataset used for the case study in the paper
# This dataset is too big to be uploaded to GitHub (with 25Mb size limit)
# Data are available from my Google Drive (https://drive.google.com/file/d/1z_DdNrHdX6zF844SfEfqt3ZOqClFoTF4/view?usp=drive_web)

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Boost-R-D: Boosting for Recurrence Data
# This version of code handles the static and dynamic features
# (Boost-R.R is the code that handles only the static features)
# Last modification: 09/26/2019

# ---------------------------------------------------------------
# Dependencies: 
# SNOW (Simple Network of Workstations) package in R for simple parallel computing
# doParallel: Foreach Parallel Adaptor for the 'parallel' Package
# foreach: Provides Foreach Looping Construct
# splines2: Regression Spline Functions and Classes
# splines: Regression Spline Functions and Classes
# spectral: Common Methods of Spectral Data Analysis (needed for visualization only)

# ---------------------------------------------------------------
# For questions, please contact:
# Xiao Liu, Dept of Industrial Engineering
# 4174 Bell Engineering Center
# University of Arkansas
# e-mail: xl027@uark.edu
# https://sites.google.com/site/liuxiaosite1/
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Stage 1: data loading and preparation
# Stage 2: running Boost-R
# Stage 3: output visualization
# ---------------------------------------------------------------
# ---------------------------------------------------------------


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Stage 1: data loading and preparation
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# set work directory
rm(list=ls())
wd= "C:/Users/xl027/Desktop/24. SurvivalTree/code_boost/Boost-R"
setwd(wd)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# load libraries
library(doParallel)
library(snow)
library(foreach)
library(optimParallel)
library("splines2")
library("splines")
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# creates a parallel computing cluster of one of the supported types
n.node = 8 # number of computing nodes in the cluster/number of threads in a single machine # use the SNOW (Simple Network of Workstations) package in R for simple parallel computing
cl <- makeCluster(n.node, type = "SOCK") 
registerDoParallel(cl)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Load data and process the data for analysis
data = read.csv("CaseStudy/data.csv",header=TRUE,stringsAsFactors = FALSE)
X = read.csv("CaseStudy/x.table.csv",header=TRUE,stringsAsFactors = FALSE)
load("CaseStudy/z.list.RData")
X = X[1:1000,]
z.list = z.list[1:1000]
data = data[1:1000,]
for (j in 1:length(z.list)){
  z.list[[j]] = z.list[[j]][,c(1)]
}
t.max = ceiling( max(apply(data, 1, max, na.rm=TRUE), na.rm=TRUE) )
n.data = nrow(data)
n.x = ncol(X)

# ---------------------------------------------------------------
# user-specified tuning parameters for Boost-R
K.value =10 # number of trees
lambda.value =300 # gamma_1 in the paper
gamma.value = 100 # gamma_2 in the paper
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# explore variables to all clusters
clusterExport(cl, "X", envir = .GlobalEnv)
clusterExport(cl, "data", envir = .GlobalEnv)
clusterExport(cl, "t.max", envir = .GlobalEnv)
clusterExport(cl, "n.data", envir = .GlobalEnv)
clusterExport(cl, "n.x", envir = .GlobalEnv)
clusterExport(cl, "K.value", envir = .GlobalEnv)
clusterExport(cl, "lambda.value", envir = .GlobalEnv)
clusterExport(cl, "gamma.value", envir = .GlobalEnv)
# ---------------------------------------------------------------


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Stage 2: running Boost-R
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# ---------------------------------------------------------------
source("sub_BoostR-D.R")
BoostR.out = BoostR2(data, X, Z=z.list, K.value, lambda.value, gamma.value, u.value=2, v.value=3, D.max=4)
# ---------------------------------------------------------------


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Stage 3: output visualization
# ---------------------------------------------------------------
# ---------------------------------------------------------------
source("main_visualization_fun_2.R")
Plot_Partition(i=3, BoostR.out)  # binary paritions of individual trees (dim(X)=2 only)
Plot_Leaf(BoostR.out) # number of leafs for all trees in the ensemble 
Plot_Individual(i=3, BoostR.out) # estimated cumulative intensity for selected individuals; i specifies the individual ID
Plot_Imp(BoostR.out, standardize=FALSE) # variable important; standardize=TRUE generates the standardized variable importance
Plot_Interaction(BoostR.out) # interaction between splines coefficients and static features (dim(X)=2 only)


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# THE END
# ---------------------------------------------------------------
# ---------------------------------------------------------------





