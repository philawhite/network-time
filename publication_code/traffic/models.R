# library(ggmap)
# library(ggplot2)
# library(maps)
# library(mapdata)
# library(sf)
# library(lubridate)
# library(scales)
library(rgdal)
library(scoringRules)
library(fields)

dat = read.csv("waic/grid_data_215.csv")[,-1]
dat$loc_ind = dat$grid_loc

dist_mat = as.matrix(rdist(dat$milepoint))
time_mat = as.matrix(rdist(dat$time_numeric))/7

cord.dec = SpatialPoints(cbind(dat$long, dat$lat), proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:2850"))
cord.UTM@coords[,1] = (cord.UTM@coords[,1] - min(cord.UTM@coords[,1]))/1000
cord.UTM@coords[,2] = (cord.UTM@coords[,2] - min(cord.UTM@coords[,2]))/1000

euc_dist = as.matrix(rdist(cbind(cord.UTM@coords[,1],cord.UTM@coords[,2])))*0.62137119224


# sort(unique(c(which(diff_mat > 0,arr.ind = TRUE)))

dat$route = factor(dat$route)

n = nrow(dat)

gc_gc_cov_simp <- nimbleFunction(
  run = function(dists = double(2),time = double(2), pars = double(1),tau2 = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        temp1 <- (1 + (time[i,j]/pars[3])^(2))
        temp <- pars[1]/temp1^(2) * (1 + (dists[i,j]/pars[2])/temp1)^(-2)
        
        result[i, j] <- temp
        result[j, i] <- temp
      }
    }
    for(i in 1:(n)){
      
      temp1 <- (1 + (time[i,i]/pars[3])^(2))
      temp <- pars[1]/temp1^(2) * (1 + (dists[i,i]/pars[2])/temp1)^(-2)
      
      result[i, i] <- temp + tau2
    }
    return(result)
  })

dagum_power_cov_simp <- nimbleFunction(     
  run = function(dists = double(2),time = double(2), pars = double(1),tau2 = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        temp1 <- (0.5 + (time[i,j]/pars[3])^(0.5))
        temp2 <- (dists[i,j]/pars[2]) / temp1
        temp <- pars[1]/temp1^2 * (1 - temp2^0.5/(1 + temp2^0.5))
        
        result[i, j] <- temp
        result[j, i] <- temp
      }
    }
    for(i in 1:(n)){
      
      temp1 <- (0.5 + (time[i,i]/pars[3])^(0.5))
      temp2 <- (dists[i,i]/pars[2]) / temp1
      temp <- pars[1]/temp1^2 * (1 - temp2^0.5/(1 + temp2^0.5))
      result[i, i] <- temp + tau2
    }
    return(result)
  })


# chol_out = chol(gc_gc_cov_simp(dist_mat,time_mat,c(1,5,5),0))



cdagum_power_cov <- compileNimble(dagum_power_cov_simp)
cgc_gc_cov <- compileNimble(gc_gc_cov_simp)


dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(), 
                 zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dpois(x, lambda, log = TRUE) + log(1 - zeroProb))
      ## or the probability if log = FALSE
      else return((1 - zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1 - zeroProb) * dpois(0, lambda, log = FALSE)
    if (log) return(log(totalProbZero))
    return(totalProbZero)
  })


rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rpois(1, lambda))
  })


dZINB <- nimbleFunction(
  run = function(x = integer(), lambda = double(),r = double(), 
                 zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    
    
    p = r/(r+lambda) 
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dnbinom(x, r,p, log = TRUE) + log(1 - zeroProb))
      ## or the probability if log = FALSE
      else return((1 - zeroProb) * dnbinom(x, r,p, log = TRUE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1 - zeroProb) * dnbinom(0, r,p, log = TRUE)
    if (log) return(log(totalProbZero))
    return(totalProbZero)
  })


rZINB <- nimbleFunction(
  run = function(n = integer(), lambda = double(),r = double(), zeroProb = double()) {
    p = r/(r+lambda) 
    
    
    returnType(integer())
    isStructuralZero <- rnbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rnbinom(1,r,p))
  })


###### Generalized Cauchy -- Network

gc_nb <- nimbleCode({
  
  # alpha ~ dnorm(0, sd = 100)
  # rho ~ dunif(0, 10000)
  
  Sigma[1:n, 1:n] <- gc_gc_cov_simp(dists[1:n,1:n], time[1:n,1:n],
                                    cov_pars[1:n_cov],0)
  # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
  
  # for(i in 1:p){
  beta ~ dnorm(0, sd = 10)
  # }
  
  # tau2 ~ T(dnorm(0,sd = 10),0,)
  # sigma2 ~ T(dnorm(0,sd = 10),0,)
  log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
  
  # sigma_beta ~ dunif(0, 1000)
  
  # linpred[1:n] <- beta
  
  w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
  r ~ dunif(0,50)
  
  for(i in 1:n){
    # mu[i] <- alpha + beta + w[i]
    mu[i] <- exp(beta + w[i])
    prob[i] <- r/(r+mu[i]) 
    
    y[i] ~ dnegbin(prob[i],r)
  }
})


gc_pois <- nimbleCode({
  
  # alpha ~ dnorm(0, sd = 100)
  # rho ~ dunif(0, 10000)
  
  Sigma[1:n, 1:n] <- gc_gc_cov_simp(dists[1:n,1:n], time[1:n,1:n],
                                    cov_pars[1:n_cov],0)
  # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
  
  # for(i in 1:p){
  beta ~ dnorm(0, sd = 10)
  # }
  
  # tau2 ~ T(dnorm(0,sd = 10),0,)
  # sigma2 ~ T(dnorm(0,sd = 10),0,)
  log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
  
  # sigma_beta ~ dunif(0, 1000)
  
  # linpred[1:n] <- beta
  
  w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
  
  
  for(i in 1:n){
    # mu[i] <- alpha + beta + w[i]
    mu[i] <- exp(beta + w[i])
    y[i] ~ dpois(mu[i])
  }
  
})

gc_ZIP <- nimbleCode({
  
  # alpha ~ dnorm(0, sd = 100)
  # rho ~ dunif(0, 10000)
  
  Sigma[1:n, 1:n] <- gc_gc_cov_simp(dists[1:n,1:n], time[1:n,1:n],
                                    cov_pars[1:n_cov],0)
  # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
  
  # for(i in 1:p){
  beta ~ dnorm(0, sd = 10)
  # }
  
  # tau2 ~ T(dnorm(0,sd = 10),0,)
  # sigma2 ~ T(dnorm(0,sd = 10),0,)
  log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
  
  # sigma_beta ~ dunif(0, 1000)
  
  # linpred[1:n] <- beta
  
  w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
  
  zero_prob ~ dunif(0,1)
  
  for(i in 1:n){
    # mu[i] <- alpha + beta + w[i]
    mu[i] <- exp(beta + w[i])
    y[i] ~ dZIP(mu[i],zero_prob)
  }
  
})

gc_ZINB <- nimbleCode({
  
  # alpha ~ dnorm(0, sd = 100)
  # rho ~ dunif(0, 10000)
  
  Sigma[1:n, 1:n] <- gc_gc_cov_simp(dists[1:n,1:n], time[1:n,1:n],
                                    cov_pars[1:n_cov],0)
  # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
  
  # for(i in 1:p){
  beta ~ dnorm(0, sd = 10)
  # }
  
  # tau2 ~ T(dnorm(0,sd = 10),0,)
  # sigma2 ~ T(dnorm(0,sd = 10),0,)
  log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
  
  # sigma_beta ~ dunif(0, 1000)
  
  # # linpred[1:n] <- beta
  
  w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
  r ~ dunif(0,50)
  
  zero_prob ~ dunif(0,1)
  
  for(i in 1:n){
    # mu[i] <- alpha + beta + w[i]
    mu[i] <- exp(beta + w[i])
    y[i] ~ dZINB(mu[i],r,zero_prob)
  }
  
})


###### Generalized Cauchy -- Euclidean

# euc_gc_nb <- nimbleCode({
#   
#   # alpha ~ dnorm(0, sd = 100)
#   # rho ~ dunif(0, 10000)
#   
#   Sigma[1:n, 1:n] <- gc_gc_cov_simp(dists[1:n,1:n], time[1:n,1:n],
#                                     cov_pars[1:n_cov],0)
#   # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
#   
# # for(i in 1:p){
#     beta ~ dnorm(0, sd = 10)
#   }
#   
#   # tau2 ~ T(dnorm(0,sd = 10),0,)
#   # sigma2 ~ T(dnorm(0,sd = 10),0,)
#   log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
#   
#   # sigma_beta ~ dunif(0, 1000)
#   
#   # linpred[1:n] <- beta
#   
#   w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
#   r ~ dunif(0,50)
#   
#   for(i in 1:n){
#     # mu[i] <- alpha + beta + w[i]
#     mu[i] <- exp(beta + w[i])
#     prob[i] <- r/(r+mu[i]) 
#     
#     y[i] ~ dnegbin(prob[i],r)
#   }
# })
# 
# 
# euc_gc_pois <- nimbleCode({
#   
#   # alpha ~ dnorm(0, sd = 100)
#   # rho ~ dunif(0, 10000)
#   
#   Sigma[1:n, 1:n] <- gc_gc_cov_simp(dists[1:n,1:n], time[1:n,1:n],
#                                     cov_pars[1:n_cov],0)
#   # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
#   
# # for(i in 1:p){
#     beta ~ dnorm(0, sd = 10)
#   }
#   
#   # tau2 ~ T(dnorm(0,sd = 10),0,)
#   # sigma2 ~ T(dnorm(0,sd = 10),0,)
#   log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
#   
#   # sigma_beta ~ dunif(0, 1000)
#   
#   # linpred[1:n] <- beta
#   
#   w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
#   
#   
#   for(i in 1:n){
#     # mu[i] <- alpha + beta + w[i]
#     mu[i] <- exp(beta + w[i])
#     y[i] ~ dpois(mu[i])
#   }
#   
# })
# 
# euc_gc_ZIP <- nimbleCode({
#   
#   # alpha ~ dnorm(0, sd = 100)
#   # rho ~ dunif(0, 10000)
#   
#   Sigma[1:n, 1:n] <- gc_gc_cov_simp(dists[1:n,1:n], time[1:n,1:n],
#                                     cov_pars[1:n_cov],0)
#   # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
#   
# # for(i in 1:p){
#     beta ~ dnorm(0, sd = 10)
#   }
#   
#   # tau2 ~ T(dnorm(0,sd = 10),0,)
#   # sigma2 ~ T(dnorm(0,sd = 10),0,)
#   log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
#   
#   # sigma_beta ~ dunif(0, 1000)
#   
#   # linpred[1:n] <- beta
#   
#   w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
#   
#   zero_prob ~ dunif(0,1)
#   
#   for(i in 1:n){
#     # mu[i] <- alpha + beta + w[i]
#     mu[i] <- exp(beta + w[i])
#     y[i] ~ dZIP(mu[i],zero_prob)
#   }
#   
# })
# 
# euc_gc_ZINB <- nimbleCode({
#   
#   # alpha ~ dnorm(0, sd = 100)
#   # rho ~ dunif(0, 10000)
#   
#   Sigma[1:n, 1:n] <- gc_gc_cov_simp(dists[1:n,1:n], time[1:n,1:n],
#                                     cov_pars[1:n_cov],0)
#   # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
#   
# # for(i in 1:p){
#     beta ~ dnorm(0, sd = 10)
#   }
#   
#   # tau2 ~ T(dnorm(0,sd = 10),0,)
#   # sigma2 ~ T(dnorm(0,sd = 10),0,)
#   log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
#   
#   # sigma_beta ~ dunif(0, 1000)
#   
#   # linpred[1:n] <- beta
#   
#   w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
#   r ~ dunif(0,50)
#   
#   zero_prob ~ dunif(0,1)
#   
#   for(i in 1:n){
#     # mu[i] <- alpha + beta + w[i]
#     mu[i] <- exp(beta + w[i])
#     y[i] ~ dZINB(mu[i],r,zero_prob)
#   }
#   
# # })

###### dagum_power -- Network

dag_nb <- nimbleCode({
  
  # alpha ~ dnorm(0, sd = 100)
  # rho ~ dunif(0, 10000)
  
  Sigma[1:n, 1:n] <- dagum_power_cov_simp(dists[1:n,1:n], time[1:n,1:n],
                                          cov_pars[1:n_cov],0)
  # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
  
  # for(i in 1:p){
  beta ~ dnorm(0, sd = 10)
  # }
  
  # tau2 ~ T(dnorm(0,sd = 10),0,)
  # sigma2 ~ T(dnorm(0,sd = 10),0,)
  log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
  
  # sigma_beta ~ dunif(0, 1000)
  
  # linpred[1:n] <- beta
  
  w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
  r ~ dunif(0,50)
  
  for(i in 1:n){
    # mu[i] <- alpha + beta + w[i]
    mu[i] <- exp(beta + w[i])
    prob[i] <- r/(r+mu[i]) 
    
    y[i] ~ dnegbin(prob[i],r)
  }
})


dag_pois <- nimbleCode({
  
  # alpha ~ dnorm(0, sd = 100)
  # rho ~ dunif(0, 10000)
  
  Sigma[1:n, 1:n] <- dagum_power_cov_simp(dists[1:n,1:n], time[1:n,1:n],
                                          cov_pars[1:n_cov],0)
  # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
  
  # for(i in 1:p){
  beta ~ dnorm(0, sd = 10)
  # }
  
  # tau2 ~ T(dnorm(0,sd = 10),0,)
  # sigma2 ~ T(dnorm(0,sd = 10),0,)
  log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
  
  # sigma_beta ~ dunif(0, 1000)
  
  # linpred[1:n] <- beta
  
  w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
  
  
  for(i in 1:n){
    # mu[i] <- alpha + beta + w[i]
    mu[i] <- exp(beta + w[i])
    y[i] ~ dpois(mu[i])
  }
  
})

dag_ZIP <- nimbleCode({
  
  # alpha ~ dnorm(0, sd = 100)
  # rho ~ dunif(0, 10000)
  
  Sigma[1:n, 1:n] <- dagum_power_cov_simp(dists[1:n,1:n], time[1:n,1:n],
                                          cov_pars[1:n_cov],0)
  # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
  
  # for(i in 1:p){
  beta ~ dnorm(0, sd = 10)
  # }
  
  # tau2 ~ T(dnorm(0,sd = 10),0,)
  # sigma2 ~ T(dnorm(0,sd = 10),0,)
  log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
  
  # sigma_beta ~ dunif(0, 1000)
  
  # linpred[1:n] <- beta
  
  w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
  
  zero_prob ~ dunif(0,1)
  
  for(i in 1:n){
    # mu[i] <- alpha + beta + w[i]
    mu[i] <- exp(beta + w[i])
    y[i] ~ dZIP(mu[i],zero_prob)
  }
  
})

dag_ZINB <- nimbleCode({
  
  # alpha ~ dnorm(0, sd = 100)
  # rho ~ dunif(0, 10000)
  
  Sigma[1:n, 1:n] <- dagum_power_cov_simp(dists[1:n,1:n], time[1:n,1:n],
                                          cov_pars[1:n_cov],0)
  # s[1:n_loc] ~ dmnorm(zeros[1:n_loc], cov = cov[1:n_loc, 1:n_loc])
  
  # for(i in 1:p){
  beta ~ dnorm(0, sd = 10)
  # }
  
  # tau2 ~ T(dnorm(0,sd = 10),0,)
  # sigma2 ~ T(dnorm(0,sd = 10),0,)
  log(cov_pars[1:3]) ~ dmnorm(cov_pars_mean[1:3],cov_pars_cov[1:3,1:3])
  
  # sigma_beta ~ dunif(0, 1000)
  
  # linpred[1:n] <- beta
  
  w[1:n] ~ dmnorm(zeros[1:n], cov = Sigma[1:n,1:n])
  r ~ dunif(0,50)
  
  zero_prob ~ dunif(0,1)
  
  for(i in 1:n){
    # mu[i] <- alpha + beta + w[i]
    mu[i] <- exp(beta + w[i])
    y[i] ~ dZINB(mu[i],r,zero_prob)
  }
  
})

glm_mod = glm(counts ~ 1,data = dat,family = "poisson")

res.w = log(dat$counts + 1e-1) - glm_mod$linear.predictors

X = model.matrix(glm_mod)

n_cov = 3

