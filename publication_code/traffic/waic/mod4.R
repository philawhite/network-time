library(MASS)
library(RcppArmadillo)
library(RcppEigen)
library(Rcpp)
library(geosphere)
library(nimble)
library(fields)
try(setwd("C:/Users/philaw/Box/Research/Porcu/network/code/traffic"))
rm(list = ls())

source("models.R")

constants <- list(n = n,
                  dists = euc_dist, time = time_mat, n_cov = n_cov,
                  cov_pars_mean = c(0,0,0),
                  cov_pars_cov = diag(c(10,10,10)),
                  zeros = rep(0,n))

data <- list(y = dat$counts)
attach(data)
attach(constants)


inits <- list(beta = coef(glm_mod),
              log_cov_pars = log(c(1,10,5)),
              w = res.w,
              r = 20,
              zero_prob = 1/2)


model <- nimbleModel(gc_pois,
                     constants = constants,
                     data = data,
                     inits = inits)

mcmcConf <- configureMCMC(model)
mcmcConf$removeSamplers(c('w',"log_cov_pars"))
mcmcConf$addSampler(target = c("log_cov_pars[1]"), type = 'RW')
mcmcConf$addSampler(target = c("log_cov_pars[2]"), type = 'RW')
mcmcConf$addSampler(target = c("log_cov_pars[3]"), type = 'RW')
mcmcConf$addSampler(target = c("w"), type = 'ess')

mcmcConf$addMonitors(c('beta',"log_cov_pars","w"))
mcmcConf$enableWAIC = TRUE
mcmcConf
codeMCMC <- buildMCMC(mcmcConf)
Cmodel = compileNimble(codeMCMC,model)

n_tot = 30000
n_burn = 20000
# n_tot = 200
# n_burn = 100
n_post = n_tot - n_burn
n_thin = 10
post_samples <- runMCMC(Cmodel$codeMCMC,niter = n_tot,nburnin = n_burn,
                        thin = n_thin,WAIC = TRUE) 

WAIC_out = data.frame(WAIC = post_samples$WAIC$WAIC,elpd = post_samples$WAIC$lppd,pd = post_samples$WAIC$pWAIC)

rm(list=setdiff(ls(), c("WAIC_out")))
save.image("mod5.RData")
