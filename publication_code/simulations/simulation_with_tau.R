# library(SSN)
# library(ggplot2)
library(dplyr)
library(reshape2)
# library(car)
# library(maptools)
# library(forecast)
# library(Hmisc)
library(fields)
library(lubridate)
library(nimble)
library(parallel)
library(scoringRules)
library(MASS)
library(mvtnorm)
try(setwd("C:/Users/philaw/Box/Research/Porcu/network/share_code/simulations"),silent = TRUE)
rm(list = ls())
# dat = read.csv("../../data/river/data_use.csv",na.strings=c("","NA"))
# dat$SampleDate = mdy(dat$SampleDate)
# dat$SampleDate[which.max(apply(dat[,-1],1,function(x){mean(!is.na(x))}))]
# date_use = seq(as.Date('2013-11-15'), as.Date('2013-11-24'), by = "1 day")
nt = 10
dist_mat_s = as.matrix(read.csv("dist_mat_sim.csv",header = FALSE)) 
dist_mat_euc_s = as.matrix(read.csv("dist_mat_euc_sim.csv",header = FALSE)) 

ns = nrow(dist_mat_s)

dist_mat = dist_mat_s %x% matrix(1,nt,nt) 
dist_mat_euc = dist_mat_euc_s %x% matrix(1,nt,nt) 

n = ns * nt
set.seed(1)
time_mat = as.matrix(rdist(runif(n, 0,1)))

gc_gc_cov <- nimbleFunction(     
  run = function(dists = double(2),time = double(2), pars = double(1),
                 sigma2 = double(0),tau2 = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        temp1 <- (1 + (time[i,j]/pars[4])^(2*pars[1]))
        temp <- sigma2/temp1^(pars[5] + 1) * (1 + (dists[i,j]/pars[3])/temp1)^(-pars[2])
        
        result[i, j] <- temp
        result[j, i] <- temp
      }
    }
    for(i in 1:(n)){
      
      temp1 <- (1 + (time[i,i]/pars[4])^(2*pars[1]))
      temp <- sigma2/temp1^(pars[5] + 1) * (1 + (dists[i,i]/pars[3])/temp1)^(-pars[2])
      
      result[i, i] <- temp + tau2
    }
    return(result)
  })

dagum_power_cov <- nimbleFunction(     
  run = function(dists = double(2),time = double(2), pars = double(1),
                 sigma2 = double(0),tau2 = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        temp1 <- (pars[5] + (time[i,j]/pars[4])^(2*pars[1]))
        temp2 <- (dists[i,j]/pars[3]) / temp1
        temp <- sigma2/temp1^2 * (1 - temp2^pars[2]/(1 + temp2^pars[2]))
        
        result[i, j] <- temp
        result[j, i] <- temp
      }
    }
    for(i in 1:(n)){
      
      temp1 <- (pars[5] + (time[i,i]/pars[4])^(2*pars[1]))
      temp2 <- (dists[i,i]/pars[3]) / temp1
      temp <- sigma2/temp1^2 * (1 - temp2^pars[2]/(1 + temp2^pars[2]))
      result[i, i] <- temp + tau2
    }
    return(result)
  })


tau2_true = 0.1
sig2_true = 0.9

# y = c(scale(mvrnorm(1,rep(0,n),cov_mat),scale = FALSE))


fn_optim_gc = function(par,y){
  cov_mat = gc_gc_cov(dists = dist_mat,time = time_mat, 
                      pars = c(1,2,exp(par[2]),exp(par[3]),1),
                      sigma2 = exp(par[1]),
                      tau2 = exp(par[4]))
  dmvnorm(y,rep(0,n),cov_mat,log = TRUE)
}

fn_optim_gc_euc = function(par,y){
  cov_mat = gc_gc_cov(dists = dist_mat_euc,time = time_mat, 
                      pars = c(1,2,exp(par[2]),exp(par[3]),1),
                      sigma2 = exp(par[1]),
                      tau2 = exp(par[4]))
  dmvnorm(y,rep(0,n),cov_mat,log = TRUE)
}

fn_optim_dag = function(par,y){
  cov_mat = dagum_power_cov(dists = dist_mat,time = time_mat, 
                      pars = c(1/4,1/2,exp(par[2]),exp(par[3]),1/2),
                      sigma2 = exp(par[1]),
                      tau2 = exp(par[4]))
  dmvnorm(y,rep(0,n),cov_mat,log = TRUE)
}

fn_optim_dag_euc = function(par,y){
  cov_mat = dagum_power_cov(dists = dist_mat_euc,time = time_mat, 
                      pars = c(1/4,1/2,exp(par[2]),exp(par[3]),1/2),
                      sigma2 = exp(par[1]),
                      tau2 = exp(par[4]))
  dmvnorm(y,rep(0,n),cov_mat,log = TRUE)
}

n_sim = 1e3
# n_sim = 5#1e3

cs = c(20,50,100,200)

sim_out = mclapply(1:length(cs),function(j){
  
  
  cov_pars_true = c(1,2,cs[j],.2,1)

  ### Fix cov_pars 1,2,5 and tau2
  
  cov_mat = gc_gc_cov(dists = dist_mat,time = time_mat, pars = cov_pars_true,
                      sigma2 = sig2_true,tau2 = tau2_true)
  
  true_log_par = c(log(0.9),log(cs[j]),log(.2),log(.1))
  
  
  sim_res = mclapply(1:n_sim,function(i){
    
    set.seed(i)
    y = c(scale(mvrnorm(1,rep(0,n),cov_mat),scale = FALSE))
    
    gc_network = optim(par=true_log_par,fn=fn_optim_gc,control=list(fnscale=-1),
                       y = y,method = "BFGS",hessian = TRUE)
    se = sqrt(diag(solve(-gc_network$hessian)))
    gc_network$lo_ci = gc_network$par - qnorm(0.975) * se
    gc_network$hi_ci = gc_network$par + qnorm(0.975) * se
    gc_network$cover = (true_log_par > gc_network$lo_ci) & (true_log_par < gc_network$hi_ci)
    
    gc_euc = optim(par=true_log_par,fn=fn_optim_gc_euc,control=list(fnscale=-1),
                   y = y,method = "BFGS")
    
    dag_network = optim(par=true_log_par,fn=fn_optim_dag,control=list(fnscale=-1),
                        y = y,method = "BFGS")
    
    dag_euc = optim(par=true_log_par,fn=fn_optim_dag_euc,control=list(fnscale=-1),
                    y = y,method = "BFGS")

    
    res = rbind(c(gc_network$value,exp(gc_network$par)),
                c(gc_euc$value,exp(gc_euc$par)),
                c(dag_network$value,exp(dag_network$par)),
                c(dag_euc$value,exp(dag_euc$par)),
                c(NA,gc_network$cover),
                c(NA,exp(gc_network$hi_ci) - exp(gc_network$lo_ci))
                )
    return(res)
  },mc.cores = 112)
  
  return(sim_res)
},mc.cores = 1)


# n_sim = 1e3
# cs = c(30,100,300)
# 
# sim_out
rm(list=setdiff(ls(), c("cs","n_sim","sim_out")))

save.image("dist_sim_with_tau.RData")
rm(list = ls())

# plot(post_samples$samples[,1],type = "l")
# plot(post_samples$samples[,10],type = "l")
