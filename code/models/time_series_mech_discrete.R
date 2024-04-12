model{
  # Data model ----
  for(i in 1:length(d18Of.ai)){
    d18Of.obs[i, 1] ~ dnorm(d18Of[d18Of.ai[i]], d18Of.pre[i])
    d18Of.pre[i] = 1 / (d18Of.obs[i, 2] ^ 2 + d18Of.sd ^ 2)
  }
  d18Of.sd ~ dgamma(5, 10)
  
  
  for(i in 1:length(mgcaf.ai)){
    mgcaf.obs[i, 1] ~ dnorm(mgcaf[mgcaf.ai[i]], mgcaf.pre[i])
    mgcaf.pre[i] = 1 / (mgcaf.obs[i, 2] ^ 2 + mgcaf.sd ^ 2)
  }
  mgcaf.sd ~ dgamma(2, 10)
  
  for(i in 1:length(d11BGrub.ai)){
    d11BGrub.obs[i, 1] ~ dnorm(d11BGrub[d11BGrub.ai[i]], d11BGrub.pre[i])
    d11BGrub.pre[i] = 1 / (d11BGrub.obs[i, 2] ^ 2 + d11Bf.sd ^ 2)
  }
  
  for(i in 1:length(d11BTsac.ai)){
    d11BTsac.obs[i, 1] ~ dnorm(d11BTsac[d11BTsac.ai[i]], d11BTsac.pre[i])
    d11BTsac.pre[i] = 1 / (d11BTsac.obs[i, 2] ^ 2 + d11Bf.sd ^ 2)
  }
  d11Bf.sd ~ dgamma(3, 10)
  
  for(i in 1:length(d13Cf.ai)){
    d13Cf.obs[i, 1] ~ dnorm(d13Cf[d13Cf.ai[i]], d13Cf.pre[i])
    d13Cf.pre[i] = 1 / (d13Cf.obs[i, 2] ^ 2 + d13Cf.sd ^ 2)
  }
  d13Cf.sd ~ dgamma(5, 10)
  
  for(i in 1:length(d13Cc.ai)){
    d13Cc.obs[i, 1] ~ dnorm(d13Cc[d13Cc.ai[i]], d13Cc.pre[i])
    d13Cc.pre[i] = 1 / (d13Cc.obs[i, 2] ^ 2 + d13Cc.sd ^ 2)
  }
  d13Cc.sd ~ dgamma(5, 10)
  
  for(i in 1:length(d18Oc.ai)){
    d18Oc.obs[i, 1] ~ dnorm(d18Oc[d18Oc.ai[i]], d18Oc.pre[i])
    d18Oc.pre[i] = 1 / (d18Oc.obs[i, 2] ^ 2 + d18Oc.sd ^ 2)
  }
  d18Oc.sd ~ dgamma(5, 10)
  
  for(i in 1:length(D47c.ai)){
    D47c.obs[i, 1] ~ dnorm(D47c[D47c.ai[i]], D47c.pre[i])
    D47c.pre[i] = 1 / (D47c.obs[i, 2] ^ 2 + D47c.sd ^ 2)
  }
  D47c.sd ~ dgamma(5, 5000)
  
  for(i in 1:length(ai)){  
    # Soil carbonate ----
    ## Depth to carbonate formation based on Retallack (2005) data, meters
    z[i] = (0.093 * MAP[i] + 13.12)
    z_m[i] = z[i] / 100
    
    ## Soil temperatures at depth z
    Tsoil[i] = MAT[i] + (PCQ_to[i] * sin(2 * 3.1415 * tsc[i] - z[i] / d)) / 
      exp(z[i] / d) 
    Tsoil.K[i] = Tsoil[i] + 273.15
    
    ## Potential Evapotranspiration - Hargreaves and Samani (1982) and Turc (1961)
    Tair_PCQ[i] = MAT[i] + PCQ_to[i]
    PET_PCQ_D.1[i] = ifelse(ha[i] < 0.5, 
                            0.013 * (Tair_PCQ[i] / (Tair_PCQ[i] + 15)) * (23.885 * Rs + 50) * (1 + ((0.5 - ha[i]) / 0.7)),
                            0.013 * (Tair_PCQ[i] / (Tair_PCQ[i] + 15)) * (23.885 * Rs + 50))
    PET_PCQ_D[i] = max(PET_PCQ_D.1[i], 0.01)
    PET_PCQ[i] = PET_PCQ_D[i] * 90
    
    ## AET in mm/quarter from Budyko curve - Pike (1964)
    AET_PCQ[i] = PPCQ[i] * (1 / (sqrt(1 + (1 / ((PET_PCQ[i] / (PPCQ[i])))) ^ 2)))
    
    ## Carbon isotopes ----
    ### Free air porosity
    FAP.1[i] = min((pore - ((PPCQ[i] - AET_PCQ[i]) / (L * 10 * pore))), pore - 0.05)
    FAP[i] = max(FAP.1[i], 0.01)
    
    ### Soil respiration rate 
    R_PCQ_D_m1[i] = 1.25 * exp(0.05452 * TmPCQ[i]) * PPCQ[i] / (127.77 + PPCQ[i])
    R_PCQ_D[i] = R_PCQ_D_m1[i] * f_R[i] # (gC/m2/d)
    
    ### Convert to molC/cm3/s
    R_PCQ_D.1[i] = R_PCQ_D[i] / (12.01 * 100 ^ 2)  # from gC/m2/d to molC/cm2/d
    R_PCQ_S[i] = R_PCQ_D.1[i] / (24 * 3600)  # molC/ cm2 / s
    R_PCQ_S_0[i]= R_PCQ_S[i] / (L * pore) # Quade et al. (2007)
    
    ### CO2 diffusion
    Dair[i] = 0.1369 * (Tsoil.K[i] / 273.15) ^ 1.958
    DIFC[i] = FAP[i] * tort * Dair[i]
    
    ### S(z)
    S_z_mol[i] = k ^ 2 * R_PCQ_S_0[i] / DIFC[i] * (1 - exp(-z[i] / k)) # (mol/cm3)
    S_z[i] = S_z_mol[i] * (0.08206 * Tsoil.K[i] * 10^9) # ppmv
    
    ### d13C of soil-respired CO2
    DD13_water[i] = 25.09 - 1.2 * (MAP[i] + 975) / (27.2 + 0.04 * (MAP[i] + 975))
    D13C_plant[i] = (28.26 * 0.22 * (pCO2[i] + 23.9)) / (28.26 + 0.22 * (pCO2[i] + 23.9)) - DD13_water[i] # schubert & Jahren (2015)
    d13Cr[i] = d13Ca[i] - D13C_plant[i]
    
    ### d13C of pedogenic carbonate
    d13Cs[i] = (pCO2[i] * d13Ca[i] + S_z[i] * (1.0044 * d13Cr[i] + 4.4))/(S_z[i] + pCO2[i])
    d13Cc[i] = ((1 + (11.98 - 0.12 * Tsoil[i]) / 1000) * (d13Cs[i] + 1000)) - 1000
    
    ## Oxygen isotopes ----
    ### Rainfall isotopes
    R18.p[i] = (d18.p[i] / 1000 + 1) * R18.VSMOW
    
    ### Equilibrium fractionation (Horita and Wesolowski 1994)
    alpha18.eq[i] = 1 / exp(((1.137e6 / (Tsoil.K[i] ^ 2) - 0.4156e3/Tsoil.K[i] - 2.0667) /1000))
    
    ### Atmospheric water vapor isotopes
    R18.a[i] = R18.p[i] * alpha18.eq[i]
    
    ### Soil evaporation from AET
    E1[i] = ETR[i] * AET_PCQ[i]
    E[i] = max(E1[i], 1) # minimum of 1 mm
    E_s[i] = E[i] / (1000 * 90 * 24 * 3600) # soil evaporation rate in m/sec
    
    ### Water vapor diffusivity
    es[i] = (0.611 * exp(17.502 * Tsoil[i] / (Tsoil[i] + 240.97))) * 1000 # saturated water vapor pressure from Tetens formula
    N.sat[i] = 0.01802 * es[i] / (Rgas * Tsoil.K[i]) # saturated water vapor concentration at a given temperature
    z.bar[i] = N.sat[i] * Dv.soil / (E_s[i] * rho) # penetration depth (m)
    z.ef1[i] = (1 - ha[i]) * z.bar[i] # the thickness of the water vapor phase region (m)
    z.ef[i] = max(z.ef1[i], 1e-10)
    
    ### Liquid water diffusivity (m2/s) (Easteal 1984)
    Dl[i] = exp(1.6766 + 1.6817 * (1000 / Tsoil.K[i]) - 0.5773 * (1000 / Tsoil.K[i]) ^ 2) * 1e-9 
    Dl.soil[i] = Dl[i] * pore * tort # effective diffusivity of liquid water (m2/s)
    z.hat[i] = Dl.soil[i] / E_s[i] # the decay length (mean penetration depth)
    
    ### The evaporation front
    h.ef[i] = ha[i] + z.ef[i] / z.bar[i] # humidity at the evaporation front
    R18.ef[i] = (alpha18.diff * R18.p[i] * (z.ef[i] / z.bar[i]) + 
                   ha[i] * R18.a[i]) / (h.ef[i] * alpha18.eq[i]) # isotopic composition at the evaporation front
    
    ### Isotope composition of soil water at depth z
    hs[i] = min(ha[i] + z_m[i] / z.bar[i], 1)
    z.f[i] = (pore / a.theta) * log(z_m[i] / z.ef[i]) # the modified depth function
    R18.s[i] = ifelse(z_m[i] <= z.ef[i], 
                      (alpha18.diff * R18.p[i] * z_m[i] / z.bar[i] + ha[i] * R18.a[i]) / 
                        (hs[i] * alpha18.eq[i]),
                      (R18.ef[i] - R18.p[i]) * exp(-z.f[i] / z.hat[i]) + R18.p[i])
    d18O.s[i] = ((R18.s[i] / R18.VSMOW) - 1) * 1000
    
    ### Isotope composition of soil carbonate
    alpha18_c_w_eq[i] = exp((1.61e4 / Tsoil.K[i] - 24.6) / 1000) # Wostbrock (2020)
    R18.c[i] = R18.s[i] * alpha18_c_w_eq[i]
    d18Oc[i] = (R18.c[i] / R18.VPDB - 1) * 1000
    D47c[i] = 0.0417e6 / Tsoil.K[i] ^ 2 + 0.139
    
    # Marine carbonate system ----
    ## Equilibrium constants ----
    Ks1m_st[i] = exp(2.83655 - 2307.1266 / temp[i] - 1.5529413 * 
                       (log(temp[i])) - ((0.20760841 + 4.0484 / temp[i]) * 
                                           sqrt(sal[i])) + 
                       0.0846834 * sal[i] - 0.00654208 * 
                       (sal[i]^1.5) + 
                       log(1 - (0.001005 * sal[i])))
    Ks2m_st[i] = exp(-9.226508 - 3351.6106 / temp[i] - 0.2005743 * (log(temp[i])) - 
                       ((0.106901773 + 23.9722 / temp[i]) * sqrt(sal[i])) + 
                       0.1130822 * sal[i] - 0.00846934 * 
                       (sal[i]^1.5) + 
                       log(1 - (0.001005 * sal[i])))
    logKsspcm_st[i] = ((-171.9065 - 0.077993 * temp[i] + 2839.319 / temp[i] + 
                          71.595 * (log(temp[i]) / log(10)) + 
                          (-0.77712 + 0.0028426 * temp[i] + 178.34 / temp[i]) * 
                          (sal[i]^0.5) - 0.07711 * sal[i] + 
                          0.0041249 * (sal[i]^1.5)))
    Ksspcm_st[i] = 10 ^ (logKsspcm_st[i])
    lnKsB_st[i] = ((-8966.9 - 2890.53 * sal[i] ^ 0.5 - 77.942 * 
                      sal[i] + 1.728 * sal[i] ^ 1.5 - 0.0996 * 
                      sal[i] ^ 2) / temp[i]) + 148.0248 + 137.1942 * 
      sal[i] ^ 0.5 + 1.62142 * sal[i] - 
      (24.4344 + 25.085 * sal[i] ^ 0.5 + 0.2474 * sal[i]) * 
      (log(temp[i])) + (0.053105 * sal[i] ^ 0.5 * temp[i])
    KsB_st[i] = exp(lnKsB_st[i])
    Ksw_st[i] = exp(148.96502 - 13847.26 / temp[i] - 23.6521 * (log(temp[i])) + 
                      (118.67 / temp[i] - 5.977 + 1.0495 * (log(temp[i]))) * 
                      (sal[i] ^ 0.5) - 0.01615 * sal[i])
    K0[i] = exp(9345.17 / temp[i] - 60.2409 + 23.3585 * (log(temp[i] / 100)) + 
                  sal[i] * (0.023517 - 0.00023656 * temp[i] + 
                              0.0047036 * ((temp[i] / 100) ^ 2)))
    
    ### Adjust for the effect of pressure (Millero 1995)
    delV1[i] = (-25.50) + 0.1271 * tempC[i]
    delV2[i] = (-15.82) + (-0.0219 * tempC[i])
    delVspc[i] = (-48.76) + (0.5304 * tempC[i])
    delVB[i] = (-29.48) + 0.1622 * tempC[i] + (2.608 / 1000) * 
      tempC[i] ^ 2
    delVw[i] = (-25.60) + 0.2324 * tempC[i] + (-3.6246 / 1000) * 
      tempC[i] ^ 2
    
    delk1[i] = (-3.08 / 1000) + (0.0877 / 1000) * tempC[i]
    delk2[i] = (1.13 / 1000) + (-0.1475 / 1000) * tempC[i]
    delkspc[i] = (-11.76 / 1000) + (0.3692 / 1000) * tempC[i]
    delkB[i] = -2.84 / 1000
    delkw[i] = (-5.13 / 1000) + (0.0794 / 1000) * tempC[i]
    
    Ks1m[i] = (exp(-((delV1[i] / (R * temp[i])) * press) + 
                     ((0.5 * delk1[i]) / (R * temp[i])) * press^2)) * Ks1m_st[i]
    Ks2m[i] = (exp(-((delV2[i] / (R * temp[i])) * press) + 
                     ((0.5 * delk2[i]) / (R * temp[i])) * press^2)) * Ks2m_st[i]
    Ksspcm[i] = (exp(-((delVspc[i] / (R * temp[i])) * press) + 
                       ((0.5 * delkspc[i]) / (R * temp[i])) * press^2)) * Ksspcm_st[i]
    KsB[i] = (exp(-((delVB[i] / (R * temp[i])) * press) + 
                    ((0.5 * delkB[i]) / (R * temp[i])) * press^2)) * KsB_st[i]
    Ksw[i] = (exp(-((delVw[i] / (R * temp[i])) * press) + 
                    ((0.5 * delkw[i]) / (R * temp[i])) * press^2)) * Ksw_st[i]
    
    ### Correct for past seawater [Ca], [Mg], and [SO4] following ZT19
    Ks1[i] = Ks1m[i] * (1 + (5e-3 * (xca[i] / xcam - 1) + 
                               1.7e-2 * (xmg[i] / xmgm - 1) + 
                               0.208 * (xso4[i] / xso4m - 1)))
    Ks2[i] = Ks2m[i] * (1 + (0.157 * (xca[i] / xcam - 1) + 
                               0.42 * (xmg[i] / xmgm - 1) + 
                               0.176 * (xso4[i] / xso4m - 1)))
    Ksspc[i] = Ksspcm[i] * (1 + (0.185 * (xca[i] / xcam - 1) + 
                                   0.518 * (xmg[i] / xmgm - 1) + 
                                   0.106 * (xso4[i] / xso4m - 1)))
    
    ## Compute pH from DIC and CO2 ----
    fco2[i] = pco2[i] * 0.9968
    co2[i] = fco2[i] * K0[i]
    hyd[i] = (2 * Ks1[i] * Ks2[i]) / 
      (-Ks1[i] + sqrt(Ks1[i] ^ 2 - 4 * Ks1[i] * Ks2[i] * (1 - dic[i] / co2[i])))
    pH[i] = -(log(hyd[i]) / log(10))
    
    ## Marine Proxies ----
    ### d11Bforam
    pKsB[i] = -(log(KsB[i]) / log(10))
    t1[i] = 10 ^ (pKsB[i] - pH[i])
    d11Bb[i] = ((t1[i] * epsilon) - (t1[i] * d11Bsw[i]) - 
                  d11Bsw[i]) / (-((t1[i] * alpha) + 1))
    d11BGrub[i] = mGrub * d11Bb[i] + bGrub
    d11BTsac[i] = mTsac * d11Bb[i] + bTsac
    
    ### d18Oforam (Bemis et al., 1998; Kim and O'Neil et al., 1997; Hollis et al., 2019)
    d18Osw.sc[i] = d18Osw[i] + (sw.sens * (sal[i] - 35))
    d18Oswpdb[i] = d18Osw.sc[i] - 0.27
    d18Of.pr[i] = d18Oswpdb[i] + 
      ((4.64 - (((4.64 ^ 2) - (4 * 0.09 * (16.1 - tempC[i]))) ^ 0.5)) / 
         (2 * 0.09))
    d18Of[i] = d18Of.pr[i] * (1 - indexop) + indexop * d18Oseccal
    
    ### d13Cforam
    d13Cf[i] = d13Ca[i] + d13Cepsilon[i]
    
    ### Mg/Caforam following Hollis et al. (2019) Mg/Ca carb chem correction approach
    mgcasw[i] = (xmg[i] / xca[i])     
    Bcorr[i] = ((mgcasw[i] ^ Hp) / (mgcaswm ^ Hp)) * Bmod
    mgca_corr[i] = Bcorr[i] * (exp(A * tempC[i]))
    mgcaf[i] = mgca_corr[i] / (1 - (sal[i] - 35) * salcorrco)
  }
  
  # Time dependent variables, time series ----
  for(i in 2:length(ai)){
    ## Derived values ----
    pCO2[i] = pco2[i] * 1e6 # atmospheric CO2 mixing ratio, ppm
    MAT[i] = tempC[i] + MAT_off[i] # mean annual terrestrial site air temperature, C 
    d18.p[i] = -15 + 0.58 * (MAT[i] * (1 - PCQ_pf[i]) +
                               TmPCQ[i] * PCQ_pf[i]) # Precipitation d18O, ppt
    PPCQ[i] = MAP[i] * PCQ_pf[i] 
    TmPCQ[i] = MAT[i] + PCQ_to[i] 
    temp[i] = tempC[i] + 273.15
    
    ## Primary environmental ----
    ### Temperature coupled to pco2 through climate sensitivity
    tempC[i] = tempC[1] + sens * log(pco2[i] / pco2[1]) / log(2) + tempC.eps[i]
    tempC.eps[i] ~ dnorm(tempC.eps[i - 1] * (tempC.phi ^ dt), tempC.pc[i])
    tempC.pc[i] = tempC.tau * ((1 - tempC.phi ^ 2) / (1 - tempC.phi ^ (2 * dt)))
    
    pco2[i] = max(min(pco2.p[i], 0.003), 0.0001)
    pco2.p[i] = pco2[i - 1] + pco2.eps[i]
    pco2.eps[i] ~ dnorm(pco2.eps[i - 1] * (pco2.phi ^ dt), pco2.pc[i])
    pco2.pc[i] = pco2.tau * ((1 - pco2.phi ^ 2) / (1 - pco2.phi ^ (2 * dt)))
    
    MAT_off[i] = MAT_off[i - 1] + MAT_off.eps[i]
    MAT_off.eps[i] ~ dnorm(MAT_off.eps[i - 1] * (MAT_off.phi ^ dt), MAT_off.pc[i])
    MAT_off.pc[i] = MAT_off.tau * ((1 - MAT_off.phi ^ 2) / (1 - MAT_off.phi ^ (2 * dt)))
    
    PCQ_to[i] = PCQ_to[i - 1] + PCQ_to.eps[i]
    PCQ_to.eps[i] ~ dnorm(PCQ_to.eps[i - 1] * (PCQ_to.phi ^ dt), PCQ_to.pc[i])
    PCQ_to.pc[i] = PCQ_to.tau * ((1 - PCQ_to.phi ^ 2) / (1 - PCQ_to.phi ^ (2 * dt)))
    
    MAP[i] = MAP[i - 1] * (1 + MAP.eps[i])
    MAP.eps[i] ~ dnorm(MAP.eps[i - 1] * (MAP.phi ^ dt), MAP.pc[i])T(-1,)
    MAP.pc[i] = MAP.tau * ((1 - MAP.phi ^ 2) / (1 - MAP.phi ^ (2 * dt)))
    
    PCQ_pf[i] = max(min(PCQ_pf.p[i], 0.25), 0.01)
    PCQ_pf.p[i] = PCQ_pf[i - 1] + PCQ_pf.eps[i]
    PCQ_pf.eps[i] ~ dnorm(PCQ_pf.eps[i - 1] * (PCQ_pf.phi ^ dt), PCQ_pf.pc[i])
    PCQ_pf.pc[i] = PCQ_pf.tau * ((1 - PCQ_pf.phi ^ 2) / (1 - PCQ_pf.phi ^ (2 * dt)))
    
    ## Secondary marine ----
    sal[i] = sal[1] * (1 + sal.eps[i])
    xca[i] = xca[1] * (1 + xca.eps[i])
    xmg[i] = xmg[1] * (1 + xmg.eps[i])
    xso4[i] = xso4[1] * (1 + xso4.eps[i])
    dic[i] = dic[1] * (1 + dic.eps[i])

    sal.eps[i] ~ dnorm(wrng[i], 1 / 1e-6 ^ 2)T(-0.1, 0.1)
    xca.eps[i] ~ dnorm(wrng[i] * 5, 1 / 1e-5 ^ 2)T(-0.2, 0.2)
    xmg.eps[i] ~ dnorm(wrng[i], 1 / 1e-6 ^ 2)T(-0.1, 0.1)
    xso4.eps[i] ~ dnorm(wrng[i], 1 / 1e-6 ^ 2)T(-0.1, 0.1)
    dic.eps[i] ~ dnorm(-wrng[i] * 5, 1 / 1e-5 ^ 2)T(-0.2, 0.2)
    
    wrng[i] = (tempC[i] - tempC[1]) * wrngK

    d11Bsw[i] = d11Bsw[i - 1] + d11Bsw.eps[i]
    d11Bsw.eps[i] ~ dnorm(d11Bsw.eps[i - 1] * (d11Bsw.phi ^ dt), d11Bsw.pc[i])
    d11Bsw.pc[i] = d11Bsw.tau * ((1 - d11Bsw.phi ^ 2) / (1 - d11Bsw.phi ^ (2 * dt)))
    
    d18Osw[i] = d18Osw[i - 1] + d18Osw.eps[i]
    d18Osw.eps[i] ~ dnorm(d18Osw.eps[i - 1] * (d18Osw.phi ^ dt), d18Osw.pc[i])
    d18Osw.pc[i] = d18Osw.tau * ((1 - d18Osw.phi ^ 2) / (1 - d18Osw.phi ^ (2 * dt)))
    
    d13Cepsilon[i] = d13Cepsilon[i - 1] + d13Cepsilon.eps[i]
    d13Cepsilon.eps[i] ~ dnorm(d13Cepsilon.eps[i - 1] * (d13Cepsilon.phi ^ dt), d13Cepsilon.pc[i])
    d13Cepsilon.pc[i] = d13Cepsilon.tau * ((1 - d13Cepsilon.phi ^ 2) / (1 - d13Cepsilon.phi ^ (2 * dt)))
    
    ## Secondary soil ----
    tsc[i] = tsc[i - 1] + tsc.eps[i]
    tsc.eps[i] ~ dnorm(tsc.eps[i - 1] * (tsc.phi ^ dt), tsc.pc[i])
    tsc.pc[i] = tsc.tau * ((1 - tsc.phi ^ 2) / (1 - tsc.phi ^ (2 * dt)))
    
    ha[i] = max(min(ha.p[i], 0.6), 0.1)
    ha.p[i] = ha[i - 1] + ha.eps[i]
    ha.eps[i] ~ dnorm(ha.eps[i - 1] * (ha.phi ^ dt), ha.pc[i])
    ha.pc[i] = ha.tau * ((1 - ha.phi ^ 2) / (1 - ha.phi ^ (2 * dt)))
    
    f_R[i] = max(min(f_R.p[i], 0.3), 0.02)
    f_R.p[i] = f_R[i - 1] + f_R.eps[i]
    f_R.eps[i] ~ dnorm(f_R.eps[i - 1] * (f_R.phi ^ dt), f_R.pc[i])
    f_R.pc[i] = f_R.tau * ((1 - f_R.phi ^ 2) / (1 - f_R.phi ^ (2 * dt)))
    
    d13Ca[i] = d13Ca[i - 1] - pco2.eps[i] * orgEff + d13Ca.eps[i]
    d13Ca.eps[i] ~ dnorm(d13Ca.eps[i - 1] * (d13Ca.phi ^ dt), d13Ca.pc[i])
    d13Ca.pc[i] = d13Ca.tau * ((1 - d13Ca.phi ^ 2) / (1 - d13Ca.phi ^ (2 * dt)))
    
    ETR[i] = max(min(ETR.p[i], 0.1), 0.01)
    ETR.p[i] = ETR[i - 1] + ETR.eps[i]
    ETR.eps[i] ~ dnorm(ETR.eps[i - 1] * (ETR.phi ^ dt), ETR.pc[i])
    ETR.pc[i] = ETR.tau * ((1 - ETR.phi ^ 2) / (1 - ETR.phi ^ (2 * dt)))
  }
  
  # Time dependent variables, ts parameters ----
  tempC.tau ~ dgamma(10, 5)
  tempC.phi ~ dbeta(2, 5)
  
  pco2.tau ~ dgamma(5, 5e-7)
  pco2.phi ~ dbeta(2, 5)
  
  MAT_off.tau ~ dgamma(10, 1)
  MAT_off.phi ~ dbeta(2, 5)
  
  PCQ_to.tau ~ dgamma(10, 1)
  PCQ_to.phi ~ dbeta(2, 5)
  
  MAP.tau ~ dgamma(10, 2e-3)
  MAP.phi ~ dbeta(2, 5)
  
  PCQ_pf.tau ~ dgamma(10, 1e-3)
  PCQ_pf.phi ~ dbeta(2, 5)
  
  wrng.tau ~ dgamma(10, 5e-4)
  wrng.phi ~ dbeta(10, 2)

  d11Bsw.tau ~ dgamma(10, 1e-2)
  d11Bsw.phi ~ dbeta(10, 2)
  
  d18Osw.tau ~ dgamma(10, 1e-2)
  d18Osw.phi ~ dbeta(10, 2)
  
  d13Cepsilon.tau ~ dgamma(10, 1e-2)
  d13Cepsilon.phi ~ dbeta(5, 2)
  
  tsc.tau ~ dgamma(10, 1e-5)
  tsc.phi ~ dbeta(2, 5)
  
  ha.tau ~ dgamma(10, 1e-2)
  ha.phi ~ dbeta(2, 5)
  
  f_R.tau ~ dgamma(10, 1e-5)
  f_R.phi ~ dbeta(2, 5)
  
  d13Ca.tau ~ dgamma(10, 1e-3)
  d13Ca.phi ~ dbeta(2, 5)
  
  ETR.tau ~ dgamma(10, 1e-5)
  ETR.phi ~ dbeta(2, 5)

  # Time dependent variables, initial conditions ----
  ## Derived values ----
  pCO2[1] = pco2[1] * 1e6 # atmospheric CO2 mixing ratio, ppm
  MAT[1] = tempC[1] + MAT_off[1] # mean annual terrestrial site air temperature, C 
  d18.p[1] = -17 + 0.58 * (MAT[1] * (1 - PCQ_pf[1]) +
                             TmPCQ[1] * PCQ_pf[1]) + MAP[1] / 250 # Precipitation d18O, ppt
  PPCQ[1] = MAP[1] * PCQ_pf[1] 
  TmPCQ[1] = MAT[1] + PCQ_to[1] 
  temp[1] = tempC[1] + 273.15
  
  ## Primary environmental ----
  tempC[1] ~ dunif(30, 38) # surface water temperature, C
  tempC.eps[1] = 0
  pco2[1] ~ dunif(0.0006, 0.0012) # atmospheric CO2 mixing ratio
  pco2.eps[1] = 0
  MAT_off[1] ~ dunif(-20, -8) # offset between terrestrial and marine temperatures, C
  MAT_off.eps[1] = 0
  PCQ_to[1] ~ dunif(-2, 14) # PCQ temperature offset, C
  PCQ_to.eps[1] = 0
  MAP[1] ~ dunif(500, 1200) # mean annual terrestrial site precipitation, mm
  MAP.eps[1] = 0
  PCQ_pf[1] ~ dunif(0.05, 0.2) # PCQ precipitation fraction
  PCQ_pf.eps[1] = 0
  
  ## Secondary marine ----
  sal[1] ~ dnorm(35, 1 / 2 ^ 2)T(25, 45) # surface water salinity, ppt
  xca[1] ~ dnorm(21, 1 / 1 ^ 2)T(14, 28) # seawater [Ca], mmol/kg
  xmg[1] ~ dnorm(68, 1 / 2 ^ 2)T(40, 90) # seawater [Mg], mmol/kg 
  xso4[1] ~ dnorm(14, 1 / 0.5 ^ 2)T(10, 18) # seawater [SO4], mmol/kg
  dic[1] ~ dnorm(0.00205, 1 / 0.0001 ^ 2)T(0.0015, 0.0025) # seawater DIC
  wrng[1] ~ dnorm(0, 1 / 0.01 ^ 2)T(-0.02, 0.02)
  d11Bsw[1] ~ dnorm(38.45, 1 / 0.5 ^ 2) # seawater d11B, ppt
  d11Bsw.eps[1] = 0
  d18Osw[1] ~ dnorm(-1.2, 1 / 0.1 ^ 2) # seawater d18O, ppt
  d18Osw.eps[1] = 0
  d13Cepsilon[1] ~ dnorm(12, 1 / 0.5 ^ 2) # offset between foram calcite and d13Catm
  d13Cepsilon.eps[1] = 0

  ## Secondary soil ----
  tsc[1] ~ dbeta(0.25 * 1000 / 0.75, 1000) # seasonal offset of PCQ for thermal diffusion
  tsc.eps[1] = 0
  ha[1] ~ dunif(0.2, 0.5)
  ha.eps[1] = 0
  f_R[1] ~ dbeta(0.15 * 500 / 0.85, 500) # ratio of PCQ to mean annual respiration rate
  f_R.eps[1] = 0
  d13Ca[1] ~ dnorm(-9, 1 / 1 ^ 2) # Atmospheric d13C, ppt
  d13Ca.eps[1] = 0
  ETR[1] ~ dbeta(0.06 * 1000 / 0.94, 1000) # Soil evaporation / AET
  ETR.eps[1] = 0
  orgEff ~ dnorm(2e4, 2e-6)
  
  # Not time dependent ----
  ## Marine ----
  press = 6 # pressure at habitation depth, bar
  indexop = 0.7
  d18Oseccal = 0.85
  sw.sens = 0.558 # regression slope of GEOSECS obs. reported in Charles and Fairbanks (1990) after Duplessey et al. (1991)
  salcorrco = 0.042
  Hp = 0.74 # nonlinearity of the relationship b/w shell and Mg/Casw (Evans and Muller 2012, T. sacculifer)
  Bmod = 0.38 # modern pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016) 
  A = 0.0757 # Exponential constant in Mg/Ca-SST calibration (Evans et al., 2016)
  mGrub = 0.62
  bGrub = 5.76
  mTsac = 0.82
  bTsac = 2.04
  sens ~ dgamma(4.5 * 5, 5)
  wrngK ~ dgamma(0.01 * 10, 10)
  
  ## Terrestrial ----
  lat = 30 # terrestrial site latitude
  Ra = 42.608 - 0.3538 * abs(lat) # total radiation at the top of the atmosphere
  Rs = Ra * 0.16 * sqrt(12) # daily temperature range assumed to be 12
  L ~ dgamma(50, 1) # mean rooting depth, cm
  k = L / 2 / log(2)   # Respiration characteristic production depth (cm) - Quade (2007)
  pore ~ dbeta(0.35 * 100 / 0.65, 100)T(0.06,) # soil porosity
  tort ~ dbeta(0.7 * 100 / 0.3, 100) # soil tortuosity
  Dv.soil = Dv.air * tort * (pore - 0.05) # effective diffusivity of water vapor in soil (m2/s)
  
  ## Constants ----
  R = 83.131 # constant (cm^3 bar mol^-1 K^-1)
  alpha = 1.0272  # B isotope fractionation factor: Klochko et al. (2006)
  epsilon = (alpha - 1) * 1000  # Compute epsilon from alpha
  d = sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
  
  ### Modern seawater chemistry
  xcam = 10.2821 # modern [Ca] (mmol kg^-1)
  xmgm = 52.8171 # modern [Mg] (mmol kg^-1)
  xso4m = 28.24  # modern [SO4] (mmol kg^-1)
  mgcaswm = xmgm/xcam # modern Mg/Ca of seawater
  
  ### Isotope ratio constants
  R13.VPDB = 0.011237
  R18.VSMOW = 0.0020052
  R18.VPDB = 0.0020672
  alpha18.diff = 1.028489
  a.theta = 0.05 # rate of increase of water content with depth (m-1) (Barnes and Allison, 1983)
  Rgas = 8.314462 # gas constant
  rho = 1000 # liquid water density (kg/m3)
  Dv.air = 2.44E-05 # water vapor diffusivity in air (m2/s) (Merlivat, 1978)
  
}


