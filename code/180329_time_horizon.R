library(R2OpenBUGS)
library(coda)
#library(lattice)
library(rjags)
library(R2jags)
#library(msm)
#library(ggmcmc)
#library(mcmcplots)
set.seed(1234)

setwd("C:/Users/u0133977/Dropbox/HypoMirror/Talks/2018/EGU/talk_jpi")
setwd("C:/Users/gjbowen/Dropbox/HypoMirror/Talks/2018/EGU/talk_jpi")

soil_model = " model {

dO_Carb ~ dnorm(dO_Carb_m, 1 / stdevO ^ 2)
dC_Carb ~ dnorm(dC_Carb_m, 1 / stdevC ^ 2)

dO_Carb_m <- (R_O_Carb / 0.0020672 - 1) * 1000
dC_Carb_m <- (R_Carb / 0.011237 - 1) * 1000

R_O_Carb <- R_O_soil * A_O
A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)

dO_soil <- (R_O_soil/0.0020052 - 1) * 1000
R_O_soil <- ((R_O_surface - R_O_P) * 2.71828 ^ (-z_m / z_i)) + R_O_P
R_O_surface <- ((1 - h) * DRF * R_O_P + h * R_O_atm) / (1 / A_atmP)
DRF <- 1 + 0.8 * (1 / 0.9723 - 1)

z_i <- DIFO / E_s
DIFO <- EPSmax * p * 2.3e-9
E_s <- E * (1/2.592e6) * 0.001

E ~ dnorm(E_m, 1 / 3 ^ 2)
E_m <- ETA * 0.3

R_Carb <- R_Soil / A_CO2_Carb
R_Soil <- (dC_Soil / 1000 + 1) * 0.011237
A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)
dC_Soil <- (((((R_hourly / DIFC) * deltaP_hat * (100 * z_100 - z_100 ^ 2 / 2) * (DIFC / DIFC13)) + pCO2 * deltaA_hat) / (R_hourly * (100 * z_100 - z_100 ^ 2 / 2) * (1 - DIFC * deltaP_hat / DIFC13) / DIFC + pCO2 * (1 - deltaA_hat))) / 0.011237 - 1) * 1000
deltaP_hat <- (deltaP / 1000 + 1) * 0.011237 / (1 + 0.011237 * (deltaP / 1000 + 1))
deltaA_hat <- (deltaA / 1000 + 1) * 0.011237 / (1 + 0.011237 * (deltaA / 1000 + 1))

deltaP <- deltaA - (deltaP_pCO2 - W)
deltaP_pCO2 ~ dnorm(deltaP_pCO2_m, 1 / 0.5 ^ 2)
deltaP_pCO2_m <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
W ~ dnorm(W_m, 1 / 0.5 ^ 2)
W_m <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))

DIFC13 <- DIFC * (1 / 1.004443)
DIFC <- EPS * p * 0.14
EPS <- ifelse(EPSmax * (ETA / CMP_mm) > EPSmax, EPSmax, EPSmax * (ETA / CMP_mm))

ETA ~ dnorm(ETA_m, 1 / 3 ^ 2)
ETA_m <- CMP_mm * (1 / (sqrt(1 + (1 / (ETP_M / CMP_mm)) ^ 2)))

ETP_M <- ETP_D * 30

ETP_D ~ dnorm(ETP_D_m, 1 / 0.2 ^ 2)
ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))

A <- 1/(2.71828^((-2.988e6 / CQT_K ^ 2 + 7.6663e3 / CQT_K - 2.4612) / 1000))

R_hourly <- R_month / 24 / 12.01
R_month ~ dnorm(R_month_m, 1 / 0.1 ^ 2)
R_month_m <- 1.25 * exp(0.05452 * CQT) * CMP_cm / (4.259 + CMP_cm)

z_m <- z_0 * 0.01
z_100 <- ifelse (z_0 < 100, z_0, 100) 
z_0 <- ifelse (z > 5 , z , 5)
z ~ dnorm(z_mean, 1 / 10 ^ 2)
z_mean <- (MAP - 257) / 2.775

dO_atm <- (R_O_atm / 0.0020052 - 1) * 1000
R_O_atm <- R_O_P / A_atmP
R_O_P <- (dO_P / 1000 + 1) * 0.0020052
A_atmP <- 2.71828 ^ ((5.9702e6 / CQT_K ^ 2 - 3.2801e4 / CQT_K + 52.227) / 1000)

RH <- h * 100
h <- ifelse(h_m2 < 0.95, h_m2, 0.95)
h_m2 ~ dnorm(h_m1, 1 / 0.05 ^ 2)
h_m1 <- 0.25 + 0.7 * (CQP / 900)

dO_P ~ dnorm(dO_P_m, 1)
dO_P_m <- -14.5 + 0.55 * MAT 

CQT_K <- CQT + 273
CQT <- MAT + T_seas
CMP_cm <- CQP / 30
CMP_mm <- CQP / 3
CQP <- MAP * P_seas
p <- 0.5 + ST * 0.3
EPSmax <- 0.1 + ST * 0.5


MAP ~ dunif(200, 2000)
P_seas ~ dbeta(2, 18)
MAT ~ dnorm(20, 1 / 5 ^ 2)I(10, 35)
T_seas ~ dnorm(8, 1 / 5 ^ 2)I(0, 20)
deltaA ~ dnorm(-6, 1 / 0.3 ^ 2)I(-7, -5)
pCO2 ~ dnorm(750, 1 / 300 ^ 2)I(200, 2000)
ST ~ dnorm(0.5, 1 / 0.05 ^ 2)I(0, 1)

Rs <- 20.35

}
"
#MAP ~ dnorm(1200, 1 / 500 ^ 2)I(200, 2000)
#P_seas ~ dnorm(0.10, 1 / 0.08 ^ 2)I(0, 0.30)
#MAT ~ dnorm(20, 1 / 5 ^ 2)I(10, 35)
#T_seas ~ dnorm(8, 1 / 5 ^ 2)I(0, 20)


clumped_model = " model {

dO_Carb ~ dnorm(dO_Carb_m, 1 / stdevO ^ 2)
dC_Carb ~ dnorm(dC_Carb_m, 1 / stdevC ^ 2)
CQTobs ~ dnorm(CQT, 1 / stdevCQTobs ^ 2)

dO_Carb_m <- (R_O_Carb / 0.0020672 - 1) * 1000
dC_Carb_m <- (R_Carb / 0.011237 - 1) * 1000

R_O_Carb <- R_O_soil * A_O
A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)

dO_soil <- (R_O_soil/0.0020052 - 1) * 1000
R_O_soil <- ((R_O_surface - R_O_P) * 2.71828 ^ (-z_m / z_i)) + R_O_P
R_O_surface <- ((1 - h) * DRF * R_O_P + h * R_O_atm) / (1 / A_atmP)
DRF <- 1 + 0.8 * (1 / 0.9723 - 1)

z_i <- DIFO / E_s
DIFO <- EPSmax * p * 2.3e-9
E_s <- E * (1/2.592e6) * 0.001

E ~ dnorm(E_m, 1 / 3 ^ 2)
E_m <- ETA * 0.3

R_Carb <- R_Soil / A_CO2_Carb
R_Soil <- (dC_Soil / 1000 + 1) * 0.011237
A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)
dC_Soil <- (((((R_hourly / DIFC) * deltaP_hat * (100 * z_100 - z_100 ^ 2 / 2) * (DIFC / DIFC13)) + pCO2 * deltaA_hat) / (R_hourly * (100 * z_100 - z_100 ^ 2 / 2) * (1 - DIFC * deltaP_hat / DIFC13) / DIFC + pCO2 * (1 - deltaA_hat))) / 0.011237 - 1) * 1000
deltaP_hat <- (deltaP / 1000 + 1) * 0.011237 / (1 + 0.011237 * (deltaP / 1000 + 1))
deltaA_hat <- (deltaA / 1000 + 1) * 0.011237 / (1 + 0.011237 * (deltaA / 1000 + 1))

deltaP <- deltaA - (deltaP_pCO2 - W)
deltaP_pCO2 ~ dnorm(deltaP_pCO2_m, 1 / 0.5 ^ 2)
deltaP_pCO2_m <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
W ~ dnorm(W_m, 1 / 0.5 ^ 2)
W_m <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))

DIFC13 <- DIFC * (1 / 1.004443)
DIFC <- EPS * p * 0.14
EPS <- ifelse(EPSmax * (ETA / CMP_mm) > EPSmax, EPSmax, EPSmax * (ETA / CMP_mm))

ETA ~ dnorm(ETA_m, 1 / 3 ^ 2)
ETA_m <- CMP_mm * (1 / (sqrt(1 + (1 / (ETP_M / CMP_mm)) ^ 2)))

ETP_M <- ETP_D * 30

ETP_D ~ dnorm(ETP_D_m, 1 / 0.2 ^ 2)
ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))

A <- 1/(2.71828^((-2.988e6 / CQT_K ^ 2 + 7.6663e3 / CQT_K - 2.4612) / 1000))

R_hourly <- R_month / 24 / 12.01
R_month ~ dnorm(R_month_m, 1 / 0.1 ^ 2)
R_month_m <- 1.25 * exp(0.05452 * CQT) * CMP_cm / (4.259 + CMP_cm)

z_m <- z_0 * 0.01
z_100 <- ifelse (z_0 < 100, z_0, 100) 
z_0 <- ifelse (z > 5 , z , 5)
z ~ dnorm(z_mean, 1 / 10 ^ 2)
z_mean <- (MAP - 257) / 2.775

dO_atm <- (R_O_atm / 0.0020052 - 1) * 1000
R_O_atm <- R_O_P / A_atmP
R_O_P <- (dO_P / 1000 + 1) * 0.0020052
A_atmP <- 2.71828 ^ ((5.9702e6 / CQT_K ^ 2 - 3.2801e4 / CQT_K + 52.227) / 1000)

RH <- h * 100
h <- ifelse(h_m2 < 0.95, h_m2, 0.95)
h_m2 ~ dnorm(h_m1, 1 / 0.05 ^ 2)
h_m1 <- 0.25 + 0.7 * (CQP / 900)

dO_P ~ dnorm(dO_P_m, 1)
dO_P_m <- -14.5 + 0.55 * MAT 

CQT_K <- CQT + 273
CQT <- MAT + T_seas
CMP_cm <- CQP / 30
CMP_mm <- CQP / 3
CQP <- MAP * P_seas
p <- 0.5 + ST * 0.3
EPSmax <- 0.1 + ST * 0.5


MAP ~ dunif(200, 2000)
P_seas ~ dbeta(2, 18)
MAT ~ dnorm(20, 1 / 5 ^ 2)I(10, 35)
T_seas ~ dnorm(8, 1 / 5 ^ 2)I(0, 20)
deltaA ~ dnorm(-6, 1 / 0.3 ^ 2)I(-7, -5)
pCO2 ~ dnorm(750, 1 / 300 ^ 2)I(200, 2000)
ST ~ dnorm(0.5, 1 / 0.05 ^ 2)I(0, 1)

Rs <- 20.35

}
"

BBNP_model = " model {

dO_Carb ~ dnorm(dO_Carb_m, 1 / stdevO ^ 2)
dC_Carb ~ dnorm(dC_Carb_m, 1 / stdevC ^ 2)
CQTobs ~ dnorm(CQT, 1 / stdevCQTobs ^ 2)

dO_Carb_m <- (R_O_Carb / 0.0020672 - 1) * 1000
dC_Carb_m <- (R_Carb / 0.011237 - 1) * 1000

R_O_Carb <- R_O_soil * A_O
A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)

dO_soil <- (R_O_soil/0.0020052 - 1) * 1000
R_O_soil <- ((R_O_surface - R_O_P) * 2.71828 ^ (-z_m / z_i)) + R_O_P
R_O_surface <- ((1 - h) * DRF * R_O_P + h * R_O_atm) / (1 / A_atmP)
DRF <- 1 + 0.8 * (1 / 0.9723 - 1)

z_i <- DIFO / E_s
DIFO <- EPSmax * p * 2.3e-9
E_s <- E * (1/2.592e6) * 0.001

E ~ dnorm(E_m, 1 / 3 ^ 2)
E_m <- ETA * 0.3

R_Carb <- R_Soil / A_CO2_Carb
R_Soil <- (dC_Soil / 1000 + 1) * 0.011237
A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)
dC_Soil <- (((((R_hourly / DIFC) * deltaP_hat * (100 * z_100 - z_100 ^ 2 / 2) * (DIFC / DIFC13)) + pCO2 * deltaA_hat) / (R_hourly * (100 * z_100 - z_100 ^ 2 / 2) * (1 - DIFC * deltaP_hat / DIFC13) / DIFC + pCO2 * (1 - deltaA_hat))) / 0.011237 - 1) * 1000
deltaP_hat <- (deltaP / 1000 + 1) * 0.011237 / (1 + 0.011237 * (deltaP / 1000 + 1))
deltaA_hat <- (deltaA / 1000 + 1) * 0.011237 / (1 + 0.011237 * (deltaA / 1000 + 1))

deltaP <- deltaA - (deltaP_pCO2 - W)
deltaP_pCO2 ~ dnorm(deltaP_pCO2_m, 1 / 0.5 ^ 2)
deltaP_pCO2_m <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
W ~ dnorm(W_m, 1 / 0.5 ^ 2)
W_m <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))

DIFC13 <- DIFC * (1 / 1.004443)
DIFC <- EPS * p * 0.14
EPS <- ifelse(EPSmax * (ETA / CMP_mm) > EPSmax, EPSmax, EPSmax * (ETA / CMP_mm))

ETA ~ dnorm(ETA_m, 1 / 3 ^ 2)
ETA_m <- CMP_mm * (1 / (sqrt(1 + (1 / (ETP_M / CMP_mm)) ^ 2)))

ETP_M <- ETP_D * 30

ETP_D ~ dnorm(ETP_D_m, 1 / 0.2 ^ 2)
ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))

A <- 1/(2.71828^((-2.988e6 / CQT_K ^ 2 + 7.6663e3 / CQT_K - 2.4612) / 1000))

R_hourly <- R_month / 24 / 12.01
R_month ~ dnorm(R_month_m, 1 / 0.1 ^ 2)
R_month_m <- 1.25 * exp(0.05452 * CQT) * CMP_cm / (4.259 + CMP_cm)

z_m <- z_0 * 0.01
z_100 <- ifelse (z_0 < 100, z_0, 100) 
z_0 <- ifelse (z > 5 , z , 5)
z ~ dnorm(z_mean, 1 / 10 ^ 2)
z_mean <- (MAP - 257) / 2.775

dO_atm <- (R_O_atm / 0.0020052 - 1) * 1000
R_O_atm <- R_O_P / A_atmP
R_O_P <- (dO_P / 1000 + 1) * 0.0020052
A_atmP <- 2.71828 ^ ((5.9702e6 / CQT_K ^ 2 - 3.2801e4 / CQT_K + 52.227) / 1000)

RH <- h * 100
h <- ifelse(h_m2 < 0.95, h_m2, 0.95)
h_m2 ~ dnorm(h_m1, 1 / 0.05 ^ 2)
h_m1 <- 0.25 + 0.7 * (CQP / 900)

dO_P ~ dnorm(dO_P_m, 1)
dO_P_m <- -14.5 + 0.55 * MAT 

CQT_K <- CQT + 273
CQT <- MAT + T_seas
CMP_cm <- CQP / 30
CMP_mm <- CQP / 3
CQP <- MAP * P_seas
p <- 0.5 + ST * 0.3
EPSmax <- 0.1 + ST * 0.5


MAP ~ dnorm(1200, 1 / 500 ^ 2)I(200, 2000)
P_seas ~ dbeta(2, 18)
MAT ~ dnorm(24, 1 / 5 ^ 2)I(10, 35)
T_seas ~ dnorm(8, 1 / 5 ^ 2)I(0, 20)
deltaA ~ dnorm(dA, 1 / 0.5 ^ 2)
pCO2 ~ dunif(500, 2000)
ST ~ dnorm(0.5, 1 / 0.1 ^ 2)I(0, 1)

Rs <- 20.35

}
"

# Define parameters that we want posterior distributions of 

parameters <- c("pCO2", "deltaA", "MAT", "T_seas", "MAP", "P_seas", "R_month", "ST")

##### Show prior/posterior

### Pretty good case
carbonate_data = list(dC_Carb=-10.5, dO_Carb=-5.0, stdevC = 0.5, stdevO = 0.5)

# Model Fitting

bayes_fit <- jags(model.file = textConnection(soil_model), parameters.to.save = parameters, 
                  data = carbonate_data, inits = NULL, 
                  n.chains=3, n.iter = 500000, n.burnin = 250000, n.thin = 50)

post = as.data.frame(bayes_fit$BUGSoutput$sims.list)

pd_MAP = density(post$MAP)
pd_MAT = density(post$MAT)
pd_Pseas = density(post$P_seas)
pd_Tseas = density(post$T_seas)

jpeg("pre_post1.jpg", res=600, units="in", width=8, height=5)
layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
par(mai=c(0.7,0.7,0.1,0.1))

plot(pd_MAP, main="")
dmn = 1200
dsd = 500
xvals=seq(dmn-2.5*dsd, dmn+2.5*dsd, dsd/10)
lines(xvals, dnorm(xvals, dmn, dsd), col="grey60")

plot(pd_MAT, main="")
dmn = 20
dsd = 5
xvals=seq(dmn-2.5*dsd, dmn+2.5*dsd, dsd/10)
lines(xvals, dnorm(xvals, dmn, dsd), col="grey60")

plot(pd_Pseas, main="")
dmn = 0
dmx = 0.4
xvals=seq(dmn, dmx, (dmx-dmn)/20)
lines(xvals, dbeta(xvals, 2, 18), col="grey60")

plot(pd_Tseas, main="")
dmn = 8
dsd = 5
xvals=seq(dmn-2.5*dsd, dmn+2.5*dsd, dsd/10)
lines(xvals, dnorm(xvals, dmn, dsd), col="grey60")

dev.off()


### Little different case
carbonate_data = list(dC_Carb=-3.5, dO_Carb=-13.0, stdevC = 0.5, stdevO = 0.5)

# Model Fitting

bayes_fit <- jags(model.file = textConnection(soil_model), parameters.to.save = parameters, data = carbonate_data, inits = NULL, 
                  n.chains=3, n.iter = 500000, n.burnin = 250000, n.thin = 50)

post = as.data.frame(bayes_fit$BUGSoutput$sims.list)

pd_MAP = density(post$MAP)
pd_MAT = density(post$MAT)
pd_Pseas = density(post$P_seas)
pd_Tseas = density(post$T_seas)

jpeg("pre_post2.jpg", res=600, units="in", width=8, height=5)
layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
par(mai=c(0.7,0.7,0.1,0.1))

plot(pd_MAP, main="")
dmn = 1200
dsd = 500
xvals=seq(dmn-2.5*dsd, dmn+2.5*dsd, dsd/10)
lines(xvals, dnorm(xvals, dmn, dsd), col="grey60")

plot(pd_MAT, main="")
dmn = 20
dsd = 5
xvals=seq(dmn-2.5*dsd, dmn+2.5*dsd, dsd/10)
lines(xvals, dnorm(xvals, dmn, dsd), col="grey60")

plot(pd_Pseas, main="")
dmn = 0
dmx = 0.4
xvals=seq(dmn, dmx, (dmx-dmn)/20)
lines(xvals, dbeta(xvals, 2, 18), col="grey60")

plot(pd_Tseas, main="")
dmn = 8
dsd = 5
xvals=seq(dmn-2.5*dsd, dmn+2.5*dsd, dsd/10)
lines(xvals, dnorm(xvals, dmn, dsd), col="grey60")

dev.off()

###Now run some sythetic data sets!
###Each synth set is made using 180406_forward_model.R, then inverted in the following loop

for(i in 1:nrow(dat)){
  carbonate_data <- list(dC_Carb=dat$dC[i], dO_Carb=dat$dO[i], stdevC = dat$sdC[i], stdevO = dat$sdO[i])
  
  # Model Fitting
  
  bayes_fit <- jags(model.file = textConnection(soil_model), parameters.to.save = parameters, data = carbonate_data, inits = NULL, 
                    n.chains=3, n.iter = 100000, n.burnin = 50000, n.thin = 50)
  
  post = as.data.frame(bayes_fit$BUGSoutput$sims.list)

  #gather up some stats on the posterior
  post.med =  sapply(post, median)
  post.low = sapply(post, quantile, probs = 0.25, names=FALSE)
  post.high = sapply(post, quantile, probs = 0.75, names=FALSE)
  
  #make space to store results for each level
  if(i == 1){
    post.med.all = post[1,]
    post.low.all = post[1,]
    post.high.all = post[1,]
  }
  
  #store results for this level
  for(j in 1:ncol(post)){
    post.med.all[i,j] = post.med[j]
    post.low.all[i,j] = post.low[j]
    post.high.all[i,j] = post.high[j]
    post.med.all[i, "Level"] = dat$level[i]
  }
}

#vector of variable names for convenience
nms = names(post.med.all)

#4 panel plot for the MAP case
jpeg("p_synth.jpg", res=600, units="in", width=8, height=5)
layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
par(mai=c(0.7,0.7,0.1,0.1))
for(i in c(1,2,3,6)){
  plot(post.med.all[,i], post.med.all$Level, type="l", xlim=c(min(post.low.all[,i]),max(post.high.all[,i])),
       xlab=nms[i], ylab = nms[10])
  polygon(c(post.low.all[,i], rev(post.high.all[,i])), c(post.med.all$Level, rev(post.med.all$Level)), col="grey70", lty=0)
  lines(post.med.all[,i], post.med.all$Level)
  lines(post.low.all[,i], post.med.all$Level)
  lines(post.high.all[,i], post.med.all$Level)
  if(i==1){lines(MAPs, post.med.all$Level, col="red")}
}
dev.off()

#4 panel plot for the MAT case
jpeg("t_synth.jpg", res=600, units="in", width=8, height=5)
layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
par(mai=c(0.7,0.7,0.1,0.1))
for(i in c(1,2,3,6)){
  plot(post.med.all[,i], post.med.all$Level, type="l", xlim=c(min(post.low.all[,i]),max(post.high.all[,i])),
       xlab=nms[i], ylab = nms[10])
  polygon(c(post.low.all[,i], rev(post.high.all[,i])), c(post.med.all$Level, rev(post.med.all$Level)), col="grey70", lty=0)
  lines(post.med.all[,i], post.med.all$Level)
  lines(post.low.all[,i], post.med.all$Level)
  lines(post.high.all[,i], post.med.all$Level)
  if(i==2){lines(MATs, post.med.all$Level, col="red")}
}
dev.off()

#4 panel plot for the T_seas case
jpeg("tseas_synth.jpg", res=600, units="in", width=8, height=5)
layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
par(mai=c(0.7,0.7,0.1,0.1))
for(i in c(1,2,3,6)){
  plot(post.med.all[,i], post.med.all$Level, type="l", xlim=c(min(post.low.all[,i]),max(post.high.all[,i])),
       xlab=nms[i], ylab = nms[10])
  polygon(c(post.low.all[,i], rev(post.high.all[,i])), c(post.med.all$Level, rev(post.med.all$Level)), col="grey70", lty=0)
  lines(post.med.all[,i], post.med.all$Level)
  lines(post.low.all[,i], post.med.all$Level)
  lines(post.high.all[,i], post.med.all$Level)
  if(i==6){lines(T_seass, post.med.all$Level, col="red")}
}
dev.off()

###Now try the MAP case with clumped
###Using synth set from 180406_forward_model.R

for(i in 1:nrow(dat)){
  carbonate_data <- list(dC_Carb=dat$dC[i], dO_Carb=dat$dO[i], CQTobs=dat$cqt[i], 
                         stdevC = dat$sdC[i], 
                         stdevO = dat$sdO[i], stdevCQTobs = dat$sdcqt[i])
  
  # Model Fitting
  
  bayes_fit <- jags(model.file = textConnection(clumped_model), parameters.to.save = parameters, data = carbonate_data, inits = NULL, 
                    n.chains=3, n.iter = 100000, n.burnin = 50000, n.thin = 50)
  
  post = as.data.frame(bayes_fit$BUGSoutput$sims.list)
  
  #gather up some stats on the posterior
  post.med =  sapply(post, median)
  post.low = sapply(post, quantile, probs = 0.25, names=FALSE)
  post.high = sapply(post, quantile, probs = 0.75, names=FALSE)
  
  #make space to store results for each level
  if(i == 1){
    post.med.all = post[1,]
    post.low.all = post[1,]
    post.high.all = post[1,]
  }
  
  #store results for this level
  for(j in 1:ncol(post)){
    post.med.all[i,j] = post.med[j]
    post.low.all[i,j] = post.low[j]
    post.high.all[i,j] = post.high[j]
    post.med.all[i, "Level"] = dat$level[i]
  }
}

#vector of variable names for convenience
nms = names(post.med.all)

#4 panel plot
jpeg("pclump_synth.jpg", res=600, units="in", width=8, height=5)
layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
par(mai=c(0.7,0.7,0.1,0.1))
for(i in c(1,2,3,6)){
  plot(post.med.all[,i], post.med.all$Level, type="l", xlim=c(min(post.low.all[,i]),max(post.high.all[,i])),
       xlab=nms[i], ylab = nms[10])
  polygon(c(post.low.all[,i], rev(post.high.all[,i])), c(post.med.all$Level, rev(post.med.all$Level)), col="grey70", lty=0)
  lines(post.med.all[,i], post.med.all$Level)
  lines(post.low.all[,i], post.med.all$Level)
  lines(post.high.all[,i], post.med.all$Level)
  if(i==1){lines(MAPs, post.med.all$Level, col="red")}
}
dev.off()

#single panel showing carbonate quarter results
jpeg("cqpclump_synth.jpg", res=600, units="in", width = 8/1.3, height = 5/1.3)
par(mai=c(1,1,0.1,0.1))
CQP = post.med.all$MAP * post.med.all$P_seas
CQP_low = post.low.all$MAP * post.low.all$P_seas
CQP_high = post.high.all$MAP * post.high.all$P_seas
plot(CQP, post.med.all$Level, type="l", xlim=c(min(CQP_low),max(CQP_high)),
     xlab="Carbonate quarter precipitation", ylab = "Level")
polygon(c(CQP_low, rev(CQP_high)), c(post.med.all$Level, rev(post.med.all$Level)), col="grey70", lty=0)
lines(CQP, post.med.all$Level)
lines(CQP_high, post.med.all$Level)
lines(CQP_low, post.med.all$Level)
lines(MAPs*0.1, post.med.all$Level, col="red")
dev.off()

###Now BBNP analysis

library(xlsx)
dat = read.xlsx("BBNP_data.xlsx", sheetIndex = 1)

for(i in 1:nrow(dat)){
  carbonate_data <- list(dC_Carb=dat$d13Ccarb[i], dO_Carb=dat$d18Ocarb[i], 
                         CQTobs=dat$T[i], 
                         stdevC = dat$d13C_SD[i], stdevO = dat$d18O_SD[i], 
                         stdevCQTobs = dat$T_95CI[i]/2, 
                         dA = dat$dC_A[i])
  
  # Model Fitting
  
  bayes_fit <- jags(model.file = textConnection(BBNP_model), parameters.to.save = parameters, data = carbonate_data, inits = NULL, 
                    n.chains=3, n.iter = 1000000, n.burnin = 500000, n.thin = 100)
  
  post = as.data.frame(bayes_fit$BUGSoutput$sims.list)
  
  #gather up some stats on the posterior
  post.med =  sapply(post, median)
  post.low = sapply(post, quantile, probs = 0.25, names=FALSE)
  post.high = sapply(post, quantile, probs = 0.75, names=FALSE)
  
  #make space to store results for each level
  if(i == 1){
    post.med.all = post[1,]
    post.low.all = post[1,]
    post.high.all = post[1,]
  }
  
  #store results for this level
  for(j in 1:ncol(post)){
    post.med.all[i,j] = post.med[j]
    post.low.all[i,j] = post.low[j]
    post.high.all[i,j] = post.high[j]
    post.med.all[i, "Level"] = dat$Age[i]
  }
}

#vector of variable names for convenience
nms = names(post.med.all)

#9 panel plot
jpeg("BBNP.jpg", res=600, units="in", width=9, height=7)
layout(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3, byrow=TRUE))
par(mai=c(0.6,0.6,0.05,0.05))
for(i in c(1,2,4,3,6,5)){
  plot(post.med.all[,i], post.med.all$Level, type="l", xlim=c(min(post.low.all[,i]),max(post.high.all[,i])),
       ylim=c(max(post.med.all$Level), min(post.med.all$Level)), xlab=nms[i], ylab = "Age (Ma)")
  polygon(c(post.low.all[,i], rev(post.high.all[,i])), c(post.med.all$Level, rev(post.med.all$Level)), col="grey70", lty=0)
  lines(post.med.all[,i], post.med.all$Level)
  lines(post.low.all[,i], post.med.all$Level)
  lines(post.high.all[,i], post.med.all$Level)
}

CQP = post.med.all$MAP * post.med.all$P_seas
CQP_low = post.low.all$MAP * post.low.all$P_seas
CQP_high = post.high.all$MAP * post.high.all$P_seas
plot(CQP, post.med.all$Level, type="l", xlim=c(min(CQP_low),max(CQP_high)),
     ylim=c(max(post.med.all$Level), min(post.med.all$Level)),
     xlab="Carbonate quarter precipitation", ylab = "Age (Ma)")
polygon(c(CQP_low, rev(CQP_high)), c(post.med.all$Level, rev(post.med.all$Level)), col="grey70", lty=0)
lines(CQP, post.med.all$Level)
lines(CQP_high, post.med.all$Level)
lines(CQP_low, post.med.all$Level)

CQT = post.med.all$MAT + post.med.all$T_seas
CQT_low = post.low.all$MAT + post.low.all$T_seas
CQT_high = post.high.all$MAT + post.high.all$T_seas
plot(CQT, post.med.all$Level, type="l", xlim=c(min(CQT_low),max(CQT_high)),
     ylim=c(max(post.med.all$Level), min(post.med.all$Level)),
     xlab="Carbonate quarter temperature", ylab = "Age (Ma)")
polygon(c(CQT_low, rev(CQT_high)), c(post.med.all$Level, rev(post.med.all$Level)), col="grey70", lty=0)
lines(CQT, post.med.all$Level)
lines(CQT_high, post.med.all$Level)
lines(CQT_low, post.med.all$Level)

i=9
plot(post.med.all[,i], post.med.all$Level, type="l", xlim=c(min(post.low.all[,i]),max(post.high.all[,i])),
     ylim=c(max(post.med.all$Level), min(post.med.all$Level)), xlab=nms[i], ylab = "Age (Ma)")
polygon(c(post.low.all[,i], rev(post.high.all[,i])), c(post.med.all$Level, rev(post.med.all$Level)), col="grey70", lty=0)
lines(post.med.all[,i], post.med.all$Level)
lines(post.low.all[,i], post.med.all$Level)
lines(post.high.all[,i], post.med.all$Level)

dev.off()

write.csv(post.med.all, "bbnp_post_med.csv", row.names = FALSE)
write.csv(post.low.all, "bbnp_post_low.csv", row.names = FALSE)
write.csv(post.high.all, "bbnp_post_high.csv", row.names = FALSE)
