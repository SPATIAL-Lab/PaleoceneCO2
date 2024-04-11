source("code/helpers.R")

# Ages
ts = ai()
ts$ts = ts$ts[-1]

#load("bigout/tsf4e3.rda")
load("bigout/ts4e3.rda")
load("bigout/tsm4e3.rda")

plot(ts$ts, post.ts$BUGSoutput$median$pCO2, type = "l")
plot(ts$ts, post.ts$BUGSoutput$mean$pCO2, type = "l")
plot(ts$ts, post.ts$BUGSoutput$mean$MAP, type = "l")
plot(ts$ts, post.ts$BUGSoutput$mean$tempC, type = "l")
plot(ts$ts, post.ts$BUGSoutput$mean$S_z, type = "l")
x = ts$ts[ts$ts_ind[[9]]]
points(x, rep(150, length(x)))
x = ts$ts[ts$ts_ind[[10]]]
points(x, rep(200, length(x)), col = "red")

plot.jpi(ts$ts, post.ts$BUGSoutput$sims.list$pCO2, xlim = c(-65, -53))
plot.jpi(ts$ts, post.tsm$BUGSoutput$sims.list$pCO2, xlim = c(-60, -53))

plot.jpi(ts$ts, post.ts$BUGSoutput$sims.list$pCO2, xlim = c(-60, -53))
par(new = TRUE)
plot(d11BGrub$age, d11BGrub$d11B, xlim = c(-60, -53), ylim = c(14, 17),
       axes = FALSE, col = "blue")
points(d11BTsac$age, d11BTsac$d11B, col = "red")
