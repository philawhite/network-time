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
rm(list = ls())


load("dist_sim2.RData")
sigma2 = 0.9
ct = 0.2

##### Which was best

best = sapply(1:4,function(i){
  sapply(1:n_sim,function(j){
    which.max(sim_out[[i]][[j]][,1])
  })
})

best_model_table = apply(best,2,function(x){
  table(factor(x,levels = c("1","2","3","4")))
})

xtable(t(best_model_table))

##### Distances between Correct Estimate

Abs_dist = lapply(1:4,function(i){
  Reduce("+", lapply(1:n_sim,function(j){
    sig_dist = abs(sim_out[[i]][[j]][,2] - sigma2 )
    cs_dist = abs(sim_out[[i]][[j]][,3] - cs[i] )
    ct_dist = abs(sim_out[[i]][[j]][,4] - ct )
    return(cbind(sig_dist =sig_dist,cs_dist = cs_dist,ct_dist = ct_dist))
  }))/ n_sim
})

Sq_dist = lapply(1:4,function(i){
  sqrt(Reduce("+", lapply(1:n_sim,function(j){
    sig_dist = (sim_out[[i]][[j]][,2] - sigma2 )^2
    cs_dist = (sim_out[[i]][[j]][,3] - cs[i] )^2
    ct_dist = (sim_out[[i]][[j]][,4] - ct )^2
    return(cbind(sig_dist =sig_dist,cs_dist = cs_dist,ct_dist = ct_dist))
  }))/ n_sim)
})

xtable(cbind(do.call(rbind,Abs_dist),
do.call(rbind,Sq_dist)),digits = 3)

library(reshape2)
library(tidyverse)

real_model1 = melt(t(sapply(1:n_sim,function(i){ 
  temp = sim_out[[1]][[i]][1,-1]
  names(temp) = c("sigma2","c_s","c_t")
  temp
}))) %>% mutate(simulation = paste('Simulation',1))

real_model2 = melt(t(sapply(1:n_sim,function(i){ 
  temp = sim_out[[2]][[i]][1,-1]
  names(temp) = c("sigma2","c_s","c_t")
  temp
}))) %>% mutate(simulation = paste('Simulation',2))


real_model3 = melt(t(sapply(1:n_sim,function(i){ 
  temp = sim_out[[3]][[i]][1,-1]
  names(temp) = c("sigma2","c_s","c_t")
  temp
}))) %>% mutate(simulation = paste('Simulation',3))

real_model4 = melt(t(sapply(1:n_sim,function(i){ 
  temp = sim_out[[4]][[i]][1,-1]
  names(temp) = c("sigma2","c_s","c_t")
  temp
}))) %>% mutate(simulation = paste('Simulation',4))



sim_out = rbind(real_model1,real_model2,
                real_model3,real_model4)


truth = data.frame(
  Var2 = c(rep(c("sigma2","c_s","c_t"),each = 4)),
  simulation = c(rep(paste("Simulation",1:4),times = 3)),
  value = c(rep(sigma2,4),cs,rep(ct,4))
)           

vnames <-list(
  'c_s' = bquote(c[S]),
  'c_t' = bquote(c[T]),  
  'sigma2' = bquote(sigma^2),
  'Simulation 1' = 'Simulation 1',
  'Simulation 2' = 'Simulation 2',
  'Simulation 3' = 'Simulation 3',
  'Simulation 4' = 'Simulation 4')

vlabeller <- function(variable,value){
  return(vnames[value])
}




pdf("identify_hist.pdf",width=12,height = 6)
ggplot() + 
  geom_histogram(data = sim_out,aes(x = value), bins = 100) +
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




