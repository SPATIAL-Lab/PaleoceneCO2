library(R2jags)

# Example datum from ~58.4 Ma

data = list(d11Bf.obs = 15.95, d11Bf.se = 0.13,
            d18Of.obs = -1.19, d18Of.se = 0.05,
            mgcaf.obs = 3.98, mgcaf.se = 0.1,
            d13Cc.obs = -9.7, d13Cc.se = 0.02,
            d18Oc.obs = -3.7, d18Oc.se = 0.1,
            D47c.obs = 0.6043, D47c.se = 0.017)

parms = c("tempC", "pCO2", "MAT", "MAP",
          "TmPCQ", "PPCQ", "d18.p", "z_m", "d18O.s", "AET_PCQ", "S_z", "d13Cr",
          "pH", "d11Bsw", "sal", "d18Osw.sc", "d18Of.pr", "mgcasw")

# Test model versions ----
## Only sample a few environmental parameters 
post.s = jags(data, NULL, parms, "code/single_sample_simplest.R", n.iter = 1e5)
View(post$BUGSoutput$summary)
save(post.s, file = "out/sss.rda")

## Vary all primary environmental parameters
post.pE = jags(data, NULL, parms, "code/single_sample_primaryEnv.R", n.iter = 5e5)
View(post.pE$BUGSoutput$summary)
save(post.pE, file = "out/sspE.rda")

plot(density(post.pE$BUGSoutput$sims.list$pCO2),
     ylim = c(0, 0.015), main = "", xlab = "")
lines(density(post.s$BUGSoutput$sims.list$pCO2), col = "red")

## Vary all marine and primary environmental parameters
post.mE = jags(data, NULL, parms, "code/single_sample_marEnv.R", n.iter = 5e5)
View(post.mE$BUGSoutput$summary)
save(post.mE, file = "out/ssmE.rda")

plot(density(post.mE$BUGSoutput$sims.list$pCO2), col = "blue", lwd = 2, 
     main = "", xlab = "", ylim = c(0, 0.015))
lines(density(post.s$BUGSoutput$sims.list$pCO2), col = "dark green", lwd = 2)

### Interesting, mgca of seawater is tighter if mg and ca are allowed to vary
plot(density(post.mE$BUGSoutput$sims.list$mgcasw), col = "blue", lwd = 2, 
     main = "", xlab = "")
lines(density(post.s$BUGSoutput$sims.list$mgcasw), col = "dark green", lwd = 2)

## Vary all environmental parameters
post.aE = jags(data, NULL, parms, "code/single_sample_allEnv.R", n.iter = 5e5)
View(post.mE$BUGSoutput$summary)
save(post.aE, file = "out/ssaE.rda")

plot(density(post.mE$BUGSoutput$sims.list$pCO2), col = "blue", lwd = 2, 
     main = "", xlab = "", ylim = c(0, 0.015))
lines(density(post.aE$BUGSoutput$sims.list$pCO2), col = "dark green", lwd = 2)

# Test data versions ----
## Marine only
data = list(d11Bf.obs = 15.95, d11Bf.se = 0.13,
            d18Of.obs = -1.19, d18Of.se = 0.05,
            mgcaf.obs = 3.98, mgcaf.se = 0.1)

post.mO = jags(data, NULL, parms, "code/single_sample_mO.R", n.iter = 5e5)
View(post.mO$BUGSoutput$summary)
save(post.mO, file = "out/ssmO.rda")

plot(density(post.mO$BUGSoutput$sims.list$pCO2), col = "blue", lwd = 2, 
     main = "", xlab = "", ylim = c(0, 0.015))
lines(density(post.aE$BUGSoutput$sims.list$pCO2), col = "dark green", lwd = 2)

## Terrestrial only
data = list(d13Cc.obs = -9.7, d13Cc.se = 0.02,
            d18Oc.obs = -3.7, d18Oc.se = 0.1,
            D47c.obs = 0.6043, D47c.se = 0.017)

post.tO = jags(data, NULL, parms, "code/single_sample_tO.R", n.iter = 5e5)
View(post.tO$BUGSoutput$summary)
save(post.tO, file = "out/sstO.rda")

plot(density(post.tO$BUGSoutput$sims.list$pCO2), col = "blue", lwd = 2, 
     main = "", xlab = "", ylim = c(0, 0.015))
lines(density(post.aE$BUGSoutput$sims.list$pCO2), col = "dark green", lwd = 2)
