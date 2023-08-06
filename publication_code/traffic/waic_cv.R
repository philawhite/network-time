library(xtable)

cv_out = read.csv("cv/all_CV.csv")
waic_out =  read.csv("waic/all_waic.csv")

cv_out
xtable(cbind(waic_out[,-1],cv_out[,-1]),digits = 3)
