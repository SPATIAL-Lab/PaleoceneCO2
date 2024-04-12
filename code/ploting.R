source("code/helpers.R")

# Ages
ts = ai()
ts$ts = ts$ts[-1]

#load("bigout/tsf4e3.rda")
load("bigout/ts4e3.rda")
load("bigout/tsm4e3.rda")

plot(ages, post.ts$BUGSoutput$mean$pCO2, type = "l")
plot(ages, post.ts$BUGSoutput$mean$tempC, type = "l")
plot(ages, post.ts$BUGSoutput$mean$MAP, type = "l")
plot(ages, post.ts$BUGSoutput$mean$PPCQ, type = "l")
plot(ages, post.ts$BUGSoutput$mean$ha, type = "l")
plot(ages, post.ts$BUGSoutput$mean$MAT, type = "l")
plot(ages, post.ts$BUGSoutput$mean$TmPCQ, type = "l")
plot(ages, post.ts$BUGSoutput$mean$S_z, type = "l")
plot(ages, post.ts$BUGSoutput$mean$z_m, type = "l")
plot(ages, post.ts$BUGSoutput$mean$d13Ca, type = "l")

plot(ages, post.ts$BUGSoutput$mean$d18Of, type = "l")
points(d18Of[, 1:2])
plot(ages, post.ts$BUGSoutput$mean$d11BGrub, type = "l")
points(d11BGrub[, 1:2])
plot(ages, post.ts$BUGSoutput$mean$d11BTsac, type = "l")
points(d11BTsac[, 1:2])
plot(ages, post.ts$BUGSoutput$mean$d13Cf, type = "l")
points(d13Cf[, 1:2])
plot(ages, post.ts$BUGSoutput$mean$d13Cc, type = "l")
points(d13Cc[, 1:2])
plot(ages, post.ts$BUGSoutput$mean$d18Oc, type = "l")
points(d18Oc[, 1:2])
plot(ages, post.ts$BUGSoutput$mean$D47c, type = "l", ylim = range(D47c$D47))
points(D47c[, 1:2])


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
