

nsynth = 500

#place to store synthetic data
dat= data.frame(level = 0, dC = 0, dO = 0, sdC = 0, sdO = 0)

#4-panel plot showing sensitivity to variables
jpeg("sensitivity.jpg", res=600, units="in", width=8, height=5)
layout(matrix(c(1,2,3,4), 2, 2, byrow=TRUE))
par(mai=c(0.6,0.6,0.1,0.1), mgp=c(2,0.5,0))

#default MAP example
#each of these blocks can also be used to generate synth data for inversion
case = 1
MAPs = seq(800,1400, length.out = 20)
MAT = rnorm(nsynth, 20, 1)
P_seas = rnorm(nsynth, 0.10, 0.02)
T_seas = rnorm(nsynth, 8, 1)
pCO2 = rnorm(nsynth, 750, 25)
dat = sm_forward()

#plot it
plot(dat$dC, dat$dO, xlab = expression(paste("Carbonate ",delta^{13}, "C (\u2030)")),
     ylab = expression(paste("Carbonate ",delta^{18}, "O (\u2030)")))

#default MAT example
case = 2
MATs = seq(10, 30, length.out = 20)
MAP = rnorm(nsynth, 1200, 100)
P_seas = rnorm(nsynth, 0.10, 0.02)
T_seas = rnorm(nsynth, 8, 1)
pCO2 = rnorm(nsynth, 750, 25)
dat = sm_forward()

#plot it
plot(dat$dC, dat$dO, xlab = expression(paste("Carbonate ",delta^{13}, "C (\u2030)")),
     ylab = expression(paste("Carbonate ",delta^{18}, "O (\u2030)")))

#default P_seas example
case = 3
P_seass = seq(0.01, 0.21, length.out = 20)
MAP = rnorm(nsynth, 1200, 100)
MAT = rnorm(nsynth, 20, 1)
T_seas = rnorm(nsynth, 8, 1)
pCO2 = rnorm(nsynth, 750, 25)
dat = sm_forward()

#plot it
plot(dat$dC, dat$dO, xlab = expression(paste("Carbonate ",delta^{13}, "C (\u2030)")),
     ylab = expression(paste("Carbonate ",delta^{18}, "O (\u2030)")))

#default T_seas example
case = 4
T_seass = seq(5, 20, length.out = 20)
MAP = rnorm(nsynth, 1200, 100)
MAT = rnorm(nsynth, 20, 1)
P_seas = rnorm(nsynth, 0.10, 0.02)
pCO2 = rnorm(nsynth, 750, 25)
dat = sm_forward()

#plot it
plot(dat$dC, dat$dO, xlab = expression(paste("Carbonate ",delta^{13}, "C (\u2030)")),
     ylab = expression(paste("Carbonate ",delta^{18}, "O (\u2030)")))

dev.off()

#pCO2 example as an extra to explore
case = 5
MAP = rnorm(nsynth, 1200, 100)
MAT = rnorm(nsynth, 20, 1)
P_seas = rnorm(nsynth, 0.10, 0.02)
T_seas = rnorm(nsynth, 8, 1)
pCO2s = seq(200,2000, length.out = 20)
dat = sm_forward()


#now do MAP but spit out clumped 'data' as well

#need to update dat to store more variables
dat= data.frame(level = 0, dC = 0, dO = 0, cqt = 0, 
                sdC = 0, sdO = 0, sdcqt = 0)

#run it
case = 1
MAPs = seq(800,1400, length.out = 20)
MAT = rnorm(nsynth, 20, 1)
P_seas = rnorm(nsynth, 0.10, 0.02)
T_seas = rnorm(nsynth, 8, 1)
dat = sm_forward_clump()



###############

sm_forward = function(){
  for(i in 1:20){
  
    #Pick parm values
  
    if(case == 1){
      MAP = MAPs[i]
    } else if(case == 2){
      MAT = MATs[i]
    } else if(case == 3){
      P_seas = P_seass[i]
    } else if(case ==4){
      T_seas = T_seass[i]
    } else if(case ==5){
      pCO2 = pCO2s[i]
    }  
      
    deltaA = rnorm(nsynth, -6, 0.3)
    ST = rnorm(nsynth, 0.5, 0.05)
    
    Rs = 20.35
    
    #Climate stuff
    EPSmax <- 0.1 + ST * 0.5
    p <- 0.5 + ST * 0.3
    CQP <- MAP * P_seas
    CMP_mm <- CQP / 3
    CMP_cm <- CQP / 30
    CQT <- MAT + T_seas
    CQT_K <- CQT + 273
    
    h_m1 <- 0.25 + 0.7 * (CQP / 900)
    h_m2 = rnorm(nsynth, h_m1, 0.05)
    h <- ifelse(h_m2 < 0.95, h_m2, 0.95)
    RH <- h * 100
    
    dO_P_m <- -14.5 + 0.55 * MAT 
    dO_P = rnorm(nsynth, dO_P_m, 1)
    
    A_atmP <- 2.71828 ^ ((5.9702e6 / CQT_K ^ 2 - 3.2801e4 / CQT_K + 52.227) / 1000)
    R_O_P <- (dO_P / 1000 + 1) * 0.0020052
    R_O_atm <- R_O_P / A_atmP
    dO_atm <- (R_O_atm / 0.0020052 - 1) * 1000
    
    z_mean <- (MAP - 257) / 2.775
    z = rnorm(nsynth, z_mean, 10)
    z_0 <- ifelse (z > 5 , z , 5)
    z_100 <- ifelse (z_0 < 100, z_0, 100) 
    z_m <- z_0 * 0.01
    
    R_month_m <- 1.25 * exp(0.05452 * CQT) * CMP_cm / (4.259 + CMP_cm)
    R_month = rnorm(nsynth, R_month_m, 0.1)
    R_hourly <- R_month / 24 / 12.01
    
    A <- 1/(2.71828^((-2.988e6 / CQT_K ^ 2 + 7.6663e3 / CQT_K - 2.4612) / 1000))
    
    ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))
    ETP_D = rnorm(nsynth, ETP_D_m, 0.2)
    
    ETP_M <- ETP_D * 30
    
    ETA_m <- CMP_mm * (1 / (sqrt(1 + (1 / (ETP_M / CMP_mm)) ^ 2)))
    ETA = rnorm(nsynth, ETA_m, 3)
    
    EPS <- ifelse(EPSmax * (ETA / CMP_mm) > EPSmax, EPSmax, EPSmax * (ETA / CMP_mm))
    DIFC <- EPS * p * 0.14
    DIFC13 <- DIFC * (1 / 1.004443)
    
    W_m <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))
    W = rnorm(nsynth, W_m, 0.5)
    deltaP_pCO2_m <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
    deltaP_pCO2 = rnorm(nsynth, deltaP_pCO2_m, 0.5)
    deltaP <- deltaA - (deltaP_pCO2 - W)
    
    deltaA_hat <- (deltaA / 1000 + 1) * 0.011237 / (1 + 0.011237 * (deltaA / 1000 + 1))
    deltaP_hat <- (deltaP / 1000 + 1) * 0.011237 / (1 + 0.011237 * (deltaP / 1000 + 1))
    dC_Soil <- (((((R_hourly / DIFC) * deltaP_hat * (100 * z_100 - z_100 ^ 2 / 2) * (DIFC / DIFC13)) + pCO2 * deltaA_hat) / (R_hourly * (100 * z_100 - z_100 ^ 2 / 2) * (1 - DIFC * deltaP_hat / DIFC13) / DIFC + pCO2 * (1 - deltaA_hat))) / 0.011237 - 1) * 1000
    A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)
    R_Soil <- (dC_Soil / 1000 + 1) * 0.011237
    R_Carb <- R_Soil / A_CO2_Carb
    
    E_m <- ETA * 0.3
    E = max(rnorm(nsynth, E_m, 3), 1)
    
    E_s <- E * (1/2.592e6) * 0.001
    DIFO <- EPSmax * p * 2.3e-9
    z_i <- DIFO / E_s
    
    DRF <- 1 + 0.8 * (1 / 0.9723 - 1)
    R_O_surface <- ((1 - h) * DRF * R_O_P + h * R_O_atm) / (1 / A_atmP)
    R_O_soil <- ((R_O_surface - R_O_P) * 2.71828 ^ (-z_m / z_i)) + R_O_P
    dO_soil <- (R_O_soil/0.0020052 - 1) * 1000
    
    A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)
    R_O_Carb <- R_O_soil * A_O
    
    dC_Carb_m <- (R_Carb / 0.011237 - 1) * 1000
    dO_Carb_m <- (R_O_Carb / 0.0020672 - 1) * 1000
  
    dat[i,] = c(i, mean(dC_Carb_m), mean(dO_Carb_m), sd(dC_Carb_m), sd(dO_Carb_m) )
    
  }
  
  return(dat)
}

###############

sm_forward_clump = function(){
  for(i in 1:20){
    
    #Pick parm values
    
    if(case == 1){
      MAP = MAPs[i]
    } else if(case == 2){
      MAT = MATs[i]
    } else if(case == 3){
      P_seas = P_seass[i]
    } else if(case ==4){
      T_seas = T_seass[i]
    }
    
    deltaA = rnorm(nsynth, -6, 0.3)
    pCO2 = rnorm(nsynth, 750, 25)
    ST = rnorm(nsynth, 0.5, 0.05)
    
    Rs = 20.35
    
    #Climate stuff
    EPSmax <- 0.1 + ST * 0.5
    p <- 0.5 + ST * 0.3
    CQP <- MAP * P_seas
    CMP_mm <- CQP / 3
    CMP_cm <- CQP / 30
    CQT <- MAT + T_seas
    CQT_K <- CQT + 273
    
    h_m1 <- 0.25 + 0.7 * (CQP / 900)
    h_m2 = rnorm(nsynth, h_m1, 0.05)
    h <- ifelse(h_m2 < 0.95, h_m2, 0.95)
    RH <- h * 100
    
    dO_P_m <- -14.5 + 0.55 * MAT 
    dO_P = rnorm(nsynth, dO_P_m, 1)
    
    A_atmP <- 2.71828 ^ ((5.9702e6 / CQT_K ^ 2 - 3.2801e4 / CQT_K + 52.227) / 1000)
    R_O_P <- (dO_P / 1000 + 1) * 0.0020052
    R_O_atm <- R_O_P / A_atmP
    dO_atm <- (R_O_atm / 0.0020052 - 1) * 1000
    
    z_mean <- (MAP - 257) / 2.775
    z = rnorm(nsynth, z_mean, 10)
    z_0 <- ifelse (z > 5 , z , 5)
    z_100 <- ifelse (z_0 < 100, z_0, 100) 
    z_m <- z_0 * 0.01
    
    R_month_m <- 1.25 * exp(0.05452 * CQT) * CMP_cm / (4.259 + CMP_cm)
    R_month = rnorm(nsynth, R_month_m, 0.1)
    R_hourly <- R_month / 24 / 12.01
    
    A <- 1/(2.71828^((-2.988e6 / CQT_K ^ 2 + 7.6663e3 / CQT_K - 2.4612) / 1000))
    
    ETP_D_m <- ifelse (RH < 50, 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50) * (1 + ((50 - RH) / 70)), 0.0133 * (CQT / (CQT + 15)) * (23.885 * Rs + 50))
    ETP_D = rnorm(nsynth, ETP_D_m, 0.2)
    
    ETP_M <- ETP_D * 30
    
    ETA_m <- CMP_mm * (1 / (sqrt(1 + (1 / (ETP_M / CMP_mm)) ^ 2)))
    ETA = rnorm(nsynth, ETA_m, 3)
    
    EPS <- ifelse(EPSmax * (ETA / CMP_mm) > EPSmax, EPSmax, EPSmax * (ETA / CMP_mm))
    DIFC <- EPS * p * 0.14
    DIFC13 <- DIFC * (1 / 1.004443)
    
    W_m <- 22.65 - (1.2 * (MAP + 975)) / (27.2 + 0.04 * (MAP + 975))
    W = rnorm(nsynth, W_m, 0.5)
    deltaP_pCO2_m <- 28.26 * 0.35 * (pCO2 + 15) / (28.26 + 0.35 * (pCO2 + 15))
    deltaP_pCO2 = rnorm(nsynth, deltaP_pCO2_m, 0.5)
    deltaP <- deltaA - (deltaP_pCO2 - W)
    
    deltaA_hat <- (deltaA / 1000 + 1) * 0.011237 / (1 + 0.011237 * (deltaA / 1000 + 1))
    deltaP_hat <- (deltaP / 1000 + 1) * 0.011237 / (1 + 0.011237 * (deltaP / 1000 + 1))
    dC_Soil <- (((((R_hourly / DIFC) * deltaP_hat * (100 * z_100 - z_100 ^ 2 / 2) * (DIFC / DIFC13)) + pCO2 * deltaA_hat) / (R_hourly * (100 * z_100 - z_100 ^ 2 / 2) * (1 - DIFC * deltaP_hat / DIFC13) / DIFC + pCO2 * (1 - deltaA_hat))) / 0.011237 - 1) * 1000
    A_CO2_Carb <- 2.71828 ^ (-2.988e3 / CQT_K ^ 2 + 7.6663 / CQT_K - 0.0024612)
    R_Soil <- (dC_Soil / 1000 + 1) * 0.011237
    R_Carb <- R_Soil / A_CO2_Carb
    
    E_m <- ETA * 0.3
    E = max(rnorm(nsynth, E_m, 3), 1)
    
    E_s <- E * (1/2.592e6) * 0.001
    DIFO <- EPSmax * p * 2.3e-9
    z_i <- DIFO / E_s
    
    DRF <- 1 + 0.8 * (1 / 0.9723 - 1)
    R_O_surface <- ((1 - h) * DRF * R_O_P + h * R_O_atm) / (1 / A_atmP)
    R_O_soil <- ((R_O_surface - R_O_P) * 2.71828 ^ (-z_m / z_i)) + R_O_P
    dO_soil <- (R_O_soil/0.0020052 - 1) * 1000
    
    A_O <- 2.71828 ^ ((2.78e6 / CQT_K ^ 2 - 2.89) / 1000)
    R_O_Carb <- R_O_soil * A_O
    
    dC_Carb_m <- (R_Carb / 0.011237 - 1) * 1000
    dO_Carb_m <- (R_O_Carb / 0.0020672 - 1) * 1000
    
    dat[i,] = c(i, mean(dC_Carb_m), mean(dO_Carb_m), mean(CQT), 
                sd(dC_Carb_m), sd(dO_Carb_m), sd(CQT))
    
  }
  
  return(dat)
}
