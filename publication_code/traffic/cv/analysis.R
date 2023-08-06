load("mod1_pred.RData")
mod1_out = c(crps_out,dss_out,MAE_out)
load("mod2_pred.RData")
mod2_out = c(crps_out,dss_out,MAE_out)
load("mod3_pred.RData")
mod3_out = c(crps_out,dss_out,MAE_out)
load("mod4_pred.RData")
mod4_out = c(crps_out,dss_out,MAE_out)
load("mod5_pred.RData")
mod5_out = c(crps_out,dss_out,MAE_out)
load("mod6_pred.RData")
mod6_out = c(crps_out,dss_out,MAE_out)


all_out = cbind(
  mod1_out,
  mod2_out,
  mod3_out,
  mod4_out,
  mod5_out,
  mod6_out
)

write.csv(all_out, "all_CV.csv")