library(R2jags)
source("code/constructors.R")

## Read data
td = read.csv("data/BBNP_data.csv")
md = read.csv("data/Shatsky_data.csv")

## Protect against overly optimistic uncertainties
td$d18O.stdev = sqrt(0.4 ^ 2 + td$d18O.stdev^2)
td$d13C.stdev = sqrt(0.4 ^ 2 + td$d13C.stdev^2)

## Remove PETM and earliest terrestrial data, make ages negative
md = md[md$age <= 55.741 | md$age >= 55.942,]
td = td[td$Age <= 59 & td$Age >= 53, ]
md$age = -md$age
td$Age = -td$Age

## Parse data into series
### Marine
d18Of = na.exclude(md[c("age", "d18O")])
d18Of$d18O.stdev = rep(0.1)
d13Cf = na.exclude(md[c("age", "d13C")])
d13Cf$d13C.stdev = rep(0.1)
mgcaf = na.exclude(md[c("age", "MgCa")])
mgcaf$MgCa.stdev = rep(0.15)
d11BGrub = d11BTsac = na.exclude(md[c("age", "d11B", "d11Bse", "species")])
d11BGrub = d11BGrub[d11BGrub$species == "Grub", -4]
d11BTsac = d11BTsac[d11BTsac$species == "Tsac", -4]

### Terrestrial
d13Cc = na.exclude(td[c("Age", "d13C", "d13C.stdev")])
d18Oc = na.exclude(td[c("Age", "d18O", "d18O.stdev")])
D47c = na.exclude(td[c("Age", "D47", "D47.stderr")])

ages = ts(d18Of$age, d13Cf$age, mgcaf$age, d11BGrub$age, d11BTsac$age,
          d13Cc$Age, d18Oc$Age, D47c$Age,
          c(d18Of$age, d13Cf$age, mgcaf$age, d11BGrub$age, d11BTsac$age),
          c(d13Cc$Age, d18Oc$Age, D47c$Age),
          seq(-59, -53, by = 0.2))
tsi = ages$ts_ind

d = list(ai = ages$ts, 
         d18Of.obs = d18Of[, 2:3], d18Of.ai = match(tsi[[1]], tsi[[9]]),
         d13Cf.obs = d13Cf[, 2:3], d13Cf.ai = match(tsi[[2]], tsi[[9]]),
         mgcaf.obs = mgcaf[, 2:3], mgcaf.ai = match(tsi[[3]], tsi[[9]]),
         d11BGrub.obs = d11BGrub[, 2:3], d11BGrub.ai = match(tsi[[4]], tsi[[9]]),
         d11BTsac.obs = d11BTsac[, 2:3], d11BTsac.ai = match(tsi[[5]], tsi[[9]]),
         d13Cc.obs = d13Cc[, 2:3], d13Cc.ai = match(tsi[[6]], tsi[[10]]),
         d18Oc.obs = d18Oc[, 2:3], d18Oc.ai = match(tsi[[7]], tsi[[10]]),
         D47c.obs = D47c[, 2:3], D47c.ai = match(tsi[[8]], tsi[[10]]),
         mar.ai = tsi[[9]], ter.ai = tsi[[10]])

parms = c("tempC", "pCO2", "MAT", "MAP", "TmPCQ", "PPCQ", "d18.p", 
          "z_m", "d18O.s", "AET_PCQ", "S_z", "d13Cr",
          "pH", "d11Bsw", "sal", "d18Osw.sc", "d18Of.pr", "mgcasw",
          "d18Of", "d13Cf", "mgcaf", "d11BGrub", "d11BTsac",
          "d13Cc", "d18Oc", "D47c")

p = proc.time()
post.ts = jags.parallel(d, NULL, parms, "code/models/time_series.R", 
                        n.iter = 2e4, n.chains = 3)
proc.time() - p

save(post.ts, file = "bigout/ts2e4.rda")
