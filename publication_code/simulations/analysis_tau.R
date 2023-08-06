library(dplyr)
library(reshape2)
library(xtable)
library(fields)
library(lubridate)
library(nimble)
library(parallel)
library(scoringRules)
library(MASS)
library(mvtnorm)
library(xtable)
rm(list = ls())


load("dist_sim_with_tau.RData")
sigma2 = 0.9
ct = 0.2
tau2 = 0.1

##### Which was best

best = sapply(1:4,function(i){
  sapply(1:n_sim,function(j){
    which.max(sim_out[[i]][[j]][1:4,1])
  })
})

best_model_table = apply(best,2,function(x){
  table(factor(x,levels = c("1","2","3","4")))
})

xtable(t(best_model_table)[,-4])

##### Distances between Correct Estimate

Abs_dist = lapply(1:4,function(i){
  Reduce("+", lapply(1:n_sim,function(j){
    sig_dist = abs(sim_out[[i]][[j]][1:2,2] - sigma2 )
    cs_dist = abs(sim_out[[i]][[j]][1:2,3] - cs[i] )
    ct_dist = abs(sim_out[[i]][[j]][1:2,4] - ct )
    tau_dist = abs(sim_out[[i]][[j]][1:2,5] - tau2 )
    
    return(cbind(sig_dist =sig_dist,cs_dist = cs_dist,
                 ct_dist = ct_dist,tau_dist = tau_dist))
  }))/ n_sim
})

Sq_dist = lapply(1:4,function(i){
  sqrt(Reduce("+", lapply(1:n_sim,function(j){
    sig_dist = (sim_out[[i]][[j]][1:2,2] - sigma2 )^2
    cs_dist = (sim_out[[i]][[j]][1:2,3] - cs[i] )^2
    ct_dist = (sim_out[[i]][[j]][1:2,4] - ct )^2
    tau_dist = (sim_out[[i]][[j]][1:2,5] - tau2 )^2
    
    return(cbind(sig_dist =sig_dist,cs_dist = cs_dist,
                 ct_dist = ct_dist,tau_dist=tau_dist))
  }))/ n_sim)
})

xtable(cbind(do.call(rbind,Abs_dist),
             do.call(rbind,Sq_dist)),digits = 3)


cover_int = lapply(1:4,function(i){
 round(t(Reduce("+", lapply(1:n_sim,function(j){
   cover =  sim_out[[i]][[j]][5,-1] 
   # width =  c(1,2,3,4)
   
   width =  sim_out[[i]][[j]][6,-1]
   
    return(cbind(cover = cover,width = width))
  }))/ n_sim),3)
})

xtable(cbind("cs",cbind(do.call(rbind,cover_int))),digits = 3)



library(reshape2)
library(tidyverse)

real_model1 = melt(t(sapply(1:n_sim,function(i){ 
  temp = sim_out[[1]][[i]][1,-1]
  names(temp) = c("sigma2","c_s","c_t","tau2")
  temp
}))) %>% mutate(simulation = paste('Simulation',1))

real_model2 = melt(t(sapply(1:n_sim,function(i){ 
  temp = sim_out[[2]][[i]][1,-1]
  names(temp) = c("sigma2","c_s","c_t","tau2")
  temp
}))) %>% mutate(simulation = paste('Simulation',2))


real_model3 = melt(t(sapply(1:n_sim,function(i){ 
  temp = sim_out[[3]][[i]][1,-1]
  names(temp) = c("sigma2","c_s","c_t","tau2")
  temp
}))) %>% mutate(simulation = paste('Simulation',3))

real_model4 = melt(t(sapply(1:n_sim,function(i){ 
  temp = sim_out[[4]][[i]][1,-1]
  names(temp) = c("sigma2","c_s","c_t","tau2")
  temp
}))) %>% mutate(simulation = paste('Simulation',4))



sim_out2 = rbind(real_model1,real_model2,
                real_model3,real_model4)


truth = data.frame(
  Var2 = c(rep(c("sigma2","c_s","c_t","tau2"),each = 4)),
  simulation = c(rep(paste("Simulation",1:4),times = 4)),
  value = c(rep(sigma2,4),cs,rep(ct,4),rep(tau2,4))
)           

vnames <-list(
  'c_s' = bquote(c[S]),
  'c_t' = bquote(c[T]),  
  'sigma2' = bquote(sigma^2),
  'tau2' = bquote(tau^2),
  'Simulation 1' = 'Simulation 1',
  'Simulation 2' = 'Simulation 2',
  'Simulation 3' = 'Simulation 3',
  'Simulation 4' = 'Simulation 4')

vlabeller <- function(variable,value){
  return(vnames[value])
}




pdf("identify_hist_tau.pdf",width=12,height = 6)
ggplot() + 
  geom_histogram(data = sim_out2,aes(x = value), bins = 100) +
  geom_vline(data = truth, mapping = aes(xintercept = value),
             linetype = "dashed",color = "red") +
  facet_grid( simulation ~ Var2,scales = "free",
              labeller = vlabeller) + 
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text =element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 15))
dev.off()


tmp = sapply(1:4,function(i){
  Reduce("+", lapply(1:n_sim,function(j){
   
    
    return(cbind(coverage = sim_out[[i]][[j]][5,-1] ))
  }))/ n_sim
})

colnames(tmp) = paste("Simulation",1:4)
rownames(tmp) = c("sigma2","c_s","c_t","tau2")
tmp
xtable(tmp[c(1,4,2,3),],digits = 3)



load("pred_sim_with_tau.RData")

lapply(1:4,function(ii){
  round(t(Reduce("+",sim_out[[ii]])/length(sim_out[[ii]])),3)
})
