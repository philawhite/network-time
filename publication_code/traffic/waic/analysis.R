load("mod1.RData")
mod1_out = WAIC_out
load("mod2.RData")
mod2_out = WAIC_out
load("mod3.RData")
mod3_out = WAIC_out
load("mod5.RData")
mod4_out = WAIC_out
load("mod6.RData")
mod5_out = WAIC_out
load("mod7.RData")
mod6_out = WAIC_out


all_out = cbind(
  mod1_out,
  mod2_out,
  mod3_out,
  mod4_out,
  mod5_out,
  mod6_out
)

write.csv(all_out, "all_waic.csv")