library(ggmap)
library(ggplot2)
library(maps)
library(mapdata)
library(sf)
library(lubridate)
library(scales)
  library(xtable)
library(MASS)
library(RcppArmadillo)
library(RcppEigen)
library(Rcpp)
library(geosphere)
library(nimble)
library(fields)

# library(parallel)
rm(list = ls())
try(setwd("C:/Users/philaw/Box/Research/Porcu/network/share_code/traffic"),silent=TRUE)

source("models.R")
load("final_zip.RData")

dim(post_samples$samples)

idx_w = grep("w",colnames(post_samples$samples))
idx_beta = grep("beta",colnames(post_samples$samples))
idx_zero = grep("zero",colnames(post_samples$samples))
idx_cov = grep("cov_pars",colnames(post_samples$samples))

tmp = post_samples$samples[,c(idx_zero,idx_beta,idx_cov)]


post_sum = t(apply(cbind(tmp[,1],(1-tmp[,1])*exp(tmp[,2] + apply(post_samples$samples[,idx_w],1,mean) + exp(tmp[,3])/2),exp(tmp[,3:5])),2,function(x){
  
  c(mean(x),sd(x),quantile(x,c(.025,.975)))
  # c(NA,NA,NA,NA)
  
}))

rownames(post_sum) = c("omega","beta_0","sig2","c_s","c_t")

xtable(post_sum,digits = 4)

random_eff_mean = apply(post_samples$samples[,idx_w],2,mean)

post_mean = apply(sweep(
  exp(post_samples$samples[,idx_w] + 
        t(X %*% t(post_samples$samples[,idx_beta]))),
  1,1-post_samples$samples[,idx_zero],"*"),2,mean)

post_sd = apply(sweep(
  exp(post_samples$samples[,idx_w] + 
                        t(X %*% t(post_samples$samples[,idx_beta]))),
  1,1-post_samples$samples[,idx_zero],"*"),2,sd)

dat_fake = dat




dat$Date <- as.POSIXct(3600 * 24 * dat$time_numeric, origin = '2015-01-01', tz = "GMT")
dat$`Space Time Random Effect` = random_eff_mean
dat$post_mean = post_mean
dat$post_sd = post_sd


zero_prob = post_samples$samples[,idx_zero]

loc_pars =  apply(exp(post_samples$samples[,idx_w] + 
                  t(X %*% t(post_samples$samples[,idx_beta]))),2,mean)

  

set.seed(1)

dat_fake$counts = (1 -rbinom(4176,1,prob = mean(zero_prob) )) * rpois(4176,loc_pars)

dat_fake = dat %>% dplyr::select(milepoint,lat,long,time_numeric,grid_loc,loc_ind,Date,counts)
write.csv(dat_fake,"fake_data_215.csv")


dat$fake = dat_fake$counts

quilt.plot(dat$milepoint,dat$Date,dat$fake,
           ny =length(unique(dat$time_numeric)),
           nx = length(unique(dat$milepoint)))


# pdf("space_time_rand.pdf",width = 8,height= 6)
# ggplot(dat, aes(x = milepoint, y=as.Date(Date),
#                 fill= `Space Time Random Effect`, 
#                 col = `Space Time Random Effect`)) +  
#   geom_tile() +
#   labs(x="Milepost", y="Date") + 
#   theme(axis.text.y=element_text(angle=45, hjust=1)) +
#   scale_y_date(date_breaks = "12 months", date_labels = "%Y",
#                expand = c(0,0)) +
#   theme(axis.text = element_text(size = 16),
#         axis.title = element_text(size = 20),
#         legend.title = element_text(size = 18)) + 
#   scale_fill_gradient2(low="navy", mid="white", high="red", 
#                        midpoint=0, limits=range(random_eff_mean))+ 
#   scale_color_gradient2(low="navy", mid="white", high="red", 
#                        midpoint=0, limits=range(random_eff_mean)) + 
#   geom_hline(yintercept = 16451 + cumsum(c(365,366,365,365,365,366)) - 13.,
#              size = 0.2) +
# geom_hline(yintercept = 16451 + sum(c(365,366,365,365,365)) - 13. + 60,
#            size = 0.75,linetype = "dashed") + 
#   guides(fill=guide_colourbar(nrow=2, byrow=TRUE,title = "Random\n Effect"),
#          color=guide_colourbar(nrow=2, byrow=TRUE,title = "Random\n Effect"))
# dev.off()


pdf("post_mean.pdf",width = 8,height= 6)
ggplot(dat, aes(x = milepoint, y=as.Date(Date),
                fill= post_mean, 
                col = post_mean)) +  
  geom_tile() +
  labs(x="Milepost", y="Date") + 
  theme(axis.text.y=element_text(angle=45, hjust=1)) +
  scale_y_date(date_breaks = "12 months", date_labels = "%Y",
               expand = c(0,0)) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) + 
  scale_fill_viridis_c(trans = "log",labels = function(x)round(x,1),
                       breaks = c(0.1,0.3,1.0,3.0,9.0))+ 
  scale_color_viridis_c(trans = "log",labels = function(x)round(x,1),
                        breaks = c(0.1,0.3,1.0,3.,9.0))+ 
  geom_hline(yintercept = 16451 + cumsum(c(365,366,365,365,365,366)) - 13.,
             size = 0.2,color = "white") +
  geom_hline(yintercept = 16451 + sum(c(365,366,365,365,365)) - 13. + 60,
             size = 0.75,linetype = "dashed",color = "white") + 
  guides(fill=guide_colourbar(nrow=2, byrow=TRUE,title = "Posterior\nMean"),
         color=guide_colourbar(nrow=2, byrow=TRUE,title = "Posterior\nMean"))
dev.off()

pdf("post_sd.pdf",width = 8,height= 6)
ggplot(dat, aes(x = milepoint, y=as.Date(Date),
                fill= post_sd, 
                col = post_sd)) +  
  geom_tile() +
  labs(x="Milepost", y="Date") + 
  theme(axis.text.y=element_text(angle=45, hjust=1)) +
  scale_y_date(date_breaks = "12 months", date_labels = "%Y",
               expand = c(0,0)) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) + 
  scale_fill_viridis_c(trans = "log",labels = function(x)round(x,2),
                       breaks =  c(.01,.03,0.1,0.3,1,2))+ 
  scale_color_viridis_c(trans = "log",labels = function(x)round(x,2),
                        breaks =  c(.01,.03,0.1,0.3,1,2))+ 
  geom_hline(yintercept = 16451 + cumsum(c(365,366,365,365,365,366)) - 13.,
             size = 0.2,color = "white") +
  geom_hline(yintercept = 16451 + sum(c(365,366,365,365,365)) - 13. + 60,
             size = 0.75,linetype = "dashed",color = "white") + 
  guides(fill=guide_colourbar(nrow=2, byrow=TRUE,title = "Posterior\nStd Dev"),
         color=guide_colourbar(nrow=2, byrow=TRUE,title = "Posterior\nStd Dev"))
dev.off()

gc_gc_cov_simp2 <- nimbleFunction(     
  run = function(dists = double(1),time = double(1), pars = double(1),tau2 = double(0)) {
    returnType(double(2))
    n <- length(dists)
    result <- matrix(nrow = n, ncol = 1, init = FALSE)
    
    for(i in 1:n){
      temp1 <- (1 + (time[i]/pars[3])^(2))
      temp <- pars[1]/temp1^(2) * (1 + (dists[i]/pars[2])/temp1)^(-2)
        
        result[i,1] <- temp
    }
    return(result)
  })

cdagum_power_cov <- compileNimble(gc_gc_cov_simp2)

dists = seq(0,5,by  = 0.1)
time = seq(0,15,by  = .25)
dt_grid = expand.grid(dists,time)

out = apply(t(sapply(1:nrow(post_samples$samples),function(xx){
c(cdagum_power_cov(dt_grid[,1], dt_grid[,2],
                   c(1,exp(post_samples$samples[xx,idx_cov[-1]])),0))
})),2,mean)



# pdf("cor_mean.pdf",width = 8,height= 6)
# par(mar = c(5,5,1,1))
# quilt.plot(dt_grid[,1],dt_grid[,2],out/max(out),nx = 61,ny = 61,
#            zlim = c(0,1),
#            col = designer.colors( n=256, col= c("white", "lightsalmon", "firebrick4"),alpha=1.0,
#            x = c(0,.3,1)),
#            xlab = "Distance (Miles)",ylab = "Time Difference (Weeks)",
#            cex.lab = 1.8)
# dev.off()



dt_grid$correlation = out#/max(out)
colnames(dt_grid) = c("Distance (Miles)","Time Difference (Weeks)","Correlation")

attach(dt_grid)
pdf("cor_mean.pdf",width = 8,height= 6)

ggplot(data = dt_grid, aes(x = `Distance (Miles)`, y = `Time Difference (Weeks)`, 
                           z = Correlation)) +
  geom_tile(aes(fill = Correlation)) +
  stat_contour(breaks = c(.05),col = "black") +
  scale_fill_gradientn(colors = c("white","lightsalmon","orangered","red","firebrick4"),
                       breaks= c(0,.2,.4,.6,.8,1),limits=c(0, 1)) +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 10))
dev.off()





contour
out/max(out) <.05

