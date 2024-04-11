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
td = td[td$Age <= 65 & td$Age >= 53, ]
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

dt = 0.5
ages = seq(-66, -53, by = dt)
d18Of.ai = get.ind(d18Of$age, ages)
d13Cf.ai = get.ind(d13Cf$age, ages)
mgcaf.ai = get.ind(mgcaf$age, ages)
d11BGrub.ai = get.ind(d11BGrub$age, ages)
d11BTsac.ai = get.ind(d11BTsac$age, ages)
d18Oc.ai = get.ind(d18Oc$Age, ages)
d13Cc.ai = get.ind(d13Cc$Age, ages)
D47c.ai = get.ind(D47c$Age, ages)

d = list(ai = ages, dt = dt,
         d18Of.obs = d18Of[, 2:3], d18Of.ai = d18Of.ai,
         d13Cf.obs = d13Cf[, 2:3], d13Cf.ai = d13Cf.ai,
         mgcaf.obs = mgcaf[, 2:3], mgcaf.ai = mgcaf.ai,
         d11BGrub.obs = d11BGrub[, 2:3], d11BGrub.ai = d11BGrub.ai,
         d11BTsac.obs = d11BTsac[, 2:3], d11BTsac.ai = d11BTsac.ai,
         d13Cc.obs = d13Cc[, 2:3], d13Cc.ai = d13Cc.ai,
         d18Oc.obs = d18Oc[, 2:3], d18Oc.ai = d18Oc.ai,
         D47c.obs = D47c[, 2:3], D47c.ai = D47c.ai)

parms = c("tempC", "pCO2", "MAT", "MAP", "TmPCQ", "PPCQ", "d18.p", 
          "z_m", "d18O.s", "AET_PCQ", "S_z", "d13Cr", "d13Cepsilon",
          "pH", "d11Bsw", "sal", "d18Osw.sc", "d18Of.pr", "mgcasw",
          "d18Of", "d13Cf", "mgcaf", "d11BGrub", "d11BTsac",
          "d13Cc", "d18Oc", "D47c")

system.time({post.ts = jags.parallel(d, NULL, parms, "code/models/time_series_discrete.R", 
                        n.iter = 2e3, n.chains = 3)})

save(post.ts, file = "bigout/ts2e3.rda")
