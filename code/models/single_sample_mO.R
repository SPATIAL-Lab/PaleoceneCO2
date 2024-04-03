model{
  
  # Data model ----
  d11Bf.obs ~ dnorm(d11Bf, d11Bf.pre)
  d11Bf.pre = 1 / d11Bf.se^2
  
  d18Of.obs ~ dnorm(d18Of, d18Of.pre)
  d18Of.pre = 1 / d18Of.se^2
  
  d13Cf.obs ~ dnorm(d13Cf, d13Cf.pre)
  d13Cf.pre = 1 / d13Cf.se^2
  
  mgcaf.obs ~ dnorm(mgcaf, mgcaf.pre)
  mgcaf.pre = 1 / mgcaf.se^2
  
  # Soil carbonate ----
  ## Depth to carbonate formation based on Retallack (2005) data, meters
  z = (0.093 * MAP + 13.12)
  z_m = z / 100
  
  ## Soil temperatures at depth z
  d = sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
  Tsoil = MAT + (PCQ_to * sin(2 * 3.1415 * tsc - z / d)) / exp(z / d) 
  Tsoil.K = Tsoil + 273.15
  
  ## Potential Evapotranspiration - Hargreaves and Samani (1982) and Turc (1961)
  PET_A_D.1 = ifelse(ha < 0.5, 
                     0.013 * (MAT / (MAT + 15)) * (23.885 * Rs + 50) * (1 + ((0.5 - ha) / 0.7)),
                     0.013 * (MAT / (MAT + 15)) * (23.885 * Rs + 50))
  PET_A_D = max(PET_A_D.1, 0.01)
  PET_A_A = PET_A_D * 365
  
  ## PET_PCQ
  Tair_PCQ = MAT + PCQ_to
  PET_PCQ_D.1 = ifelse(ha < 0.5, 
                       0.013 * (Tair_PCQ / (Tair_PCQ + 15)) * (23.885 * Rs + 50) * (1 + ((0.5 - ha) / 0.7)),
                       0.013 * (Tair_PCQ / (Tair_PCQ + 15)) * (23.885 * Rs + 50))
  PET_PCQ_D = max(PET_PCQ_D.1, 0.01)
  PET_PCQ = PET_PCQ_D * 90
  
  ## AET in mm/quarter from Budyko curve - Pike (1964)
  AET_PCQ = PPCQ * (1 / (sqrt(1 + (1 / ((PET_PCQ / (PPCQ)))) ^ 2)))
  
  ## Carbon isotopes ----
  ### Free air porosity
  FAP.1 = min((pore - ((PPCQ - AET_PCQ) / (L * 10 * pore))), pore - 0.05)
  FAP = max(FAP.1, 0.01)
  
  ### Respiration characteristic production depth (cm) - Quade (2007)
  k = L / 2 / log(2)
  
  ### Soil respiration rate 
  R_PCQ_D_m1 = 1.25 * exp(0.05452 * TmPCQ) * PPCQ / (127.77 + PPCQ)
  R_PCQ_D = R_PCQ_D_m1 * f_R # (gC/m2/d)
  
  ### Convert to molC/cm3/s
  R_PCQ_D.1 = R_PCQ_D / (12.01 * 100^2)  # from gC/m2/d to molC/cm2/d
  R_PCQ_S = R_PCQ_D.1 / (24 * 3600)  # molC/ cm2 / s
  R_PCQ_S_0 = R_PCQ_S / (L * pore) # Quade et al. (2007)
  
  ### CO2 diffusion
  Dair = 0.1369 * (Tsoil.K / 273.15) ^ 1.958
  DIFC = FAP * tort * Dair
  
  ### S(z)
  S_z_mol = k ^ 2 * R_PCQ_S_0 / DIFC * (1 - exp(-z / k)) # (mol/cm3)
  S_z = S_z_mol * (0.08206 * Tsoil.K * 10^9) # ppmv
  
  ### d13C of soil-respired CO2
  DD13_water = 25.09 - 1.2 * (MAP + 975) / (27.2 + 0.04 * (MAP + 975))
  D13C_plant = (28.26 * 0.22 * (pCO2 + 23.9)) / (28.26 + 0.22 * (pCO2 + 23.9)) - DD13_water # schubert & Jahren (2015)
  d13Cr = d13Ca - D13C_plant
  
  ### d13C of pedogenic carbonate
  d13Cs = (pCO2 * d13Ca + S_z * (1.0044 * d13Cr + 4.4))/(S_z + pCO2)
  d13Cc = ((1 + (11.98 - 0.12 * Tsoil) / 1000) * (d13Cs + 1000)) - 1000
  
  ## Oxygen isotopes ----
  ### Rainfall isotopes
  R18.p = (d18.p / 1000 + 1) * R18.VSMOW
  
  ### Equilibrium fractionation (Horita and Wesolowski 1994)
  alpha18.eq = 1 / exp(((1.137e6 / (Tsoil.K ^ 2) - 0.4156e3/Tsoil.K - 2.0667) /1000))
  
  ### Atmospheric water vapor isotopes
  R18.a = R18.p * alpha18.eq
  
  ### Soil evaporation from AET
  E1 = ETR * AET_PCQ
  E = max(E1, 1) # minimum of 1 mm
  E_s = E / (1000 * 90 * 24 * 3600) # soil evaporation rate in m/sec
  
  ### Water vapor diffusivity
  Dv.soil = Dv.air * tort * (pore - 0.05) # effective diffusivity of water vapor in soil (m2/s)
  es = (0.611 * exp(17.502 * Tsoil / (Tsoil + 240.97))) * 1000 # saturated water vapor pressure from Tetens formula
  N.sat = 0.01802 * es / (Rgas * Tsoil.K) # saturated water vapor concentration at a given temperature
  z.bar = N.sat * Dv.soil / (E_s * rho) # penetration depth (m)
  z.ef1 = (1 - ha) * z.bar # the thickness of the water vapor phase region (m)
  z.ef = max(z.ef1, 1e-10)
  
  ### Liquid water diffusivity (m2/s) (Easteal 1984)
  Dl = exp(1.6766 + 1.6817 * (1000 / Tsoil.K) - 0.5773 * (1000 / Tsoil.K)^2) * 10^-9 
  Dl.soil = Dl * pore * tort # effective diffusivity of liquid water (m2/s)
  z.hat = Dl.soil / E_s # the decay length (mean penetration depth)
  
  ### The evaporation front
  h.ef = ha + z.ef / z.bar # humidity at the evaporation front
  R18.ef = (alpha18.diff * R18.p * (z.ef / z.bar) + ha * R18.a) / (h.ef * alpha18.eq) # isotopic composition at the evaporation front
  
  ### Isotope composition of soil water at depth z
  hs = min(ha + z_m / z.bar, 1)
  z.f = (pore / a.theta) * log(z_m / z.ef) # the modified depth function
  R18.s = ifelse(z_m <= z.ef, (alpha18.diff * R18.p * z_m / z.bar + ha * R18.a) / (hs * alpha18.eq),
                 (R18.ef - R18.p) * exp(-z.f / z.hat) + R18.p)
  d18O.s = ((R18.s / R18.VSMOW) - 1) * 1000
  
  ### Isotope composition of soil carbonate
  alpha18_c_w_eq = exp((1.61e4 / Tsoil.K - 24.6) / 1000) # Wostbrock (2020)
  R18.c = R18.s * alpha18_c_w_eq
  d18Oc = (R18.c / R18.VPDB - 1) * 1000
  D47c = 0.0422e6 / Tsoil.K^2 + 0.215
  
  # Marine carbonate system ----
  ## Equilibrium constants ----
  Ks1m_st = exp(2.83655 - 2307.1266 / temp - 1.5529413 * 
                  (log(temp)) - ((0.20760841 + 4.0484 / temp) * 
                                   sqrt(sal)) + 
                  0.0846834 * sal - 0.00654208 * 
                  (sal^1.5) + 
                  log(1 - (0.001005 * sal)))
  Ks2m_st = exp(-9.226508 - 3351.6106 / temp - 0.2005743 * (log(temp)) - 
                  ((0.106901773 + 23.9722 / temp) * sqrt(sal)) + 
                  0.1130822 * sal - 0.00846934 * 
                  (sal^1.5) + 
                  log(1 - (0.001005 * sal)))
  logKsspcm_st = ((-171.9065 - 0.077993 * temp + 2839.319 / temp + 
                     71.595 * (log(temp) / log(10)) + 
                     (-0.77712 + 0.0028426 * temp + 178.34 / temp) * 
                     (sal^0.5) - 0.07711 * sal + 
                     0.0041249 * (sal^1.5)))
  Ksspcm_st = 10^(logKsspcm_st)
  lnKsB_st = ((-8966.9 - 2890.53 * sal^0.5 - 77.942 * 
                 sal + 1.728 * sal^1.5 - 0.0996 * 
                 sal^2)/temp) + 148.0248 + 137.1942 * 
    sal^0.5 + 1.62142 * sal - 
    (24.4344+25.085 * sal^0.5 + 0.2474 * sal) * 
    (log(temp)) + (0.053105 * sal^0.5 * temp)
  KsB_st = exp(lnKsB_st)
  Ksw_st = exp(148.96502 - 13847.26 / temp - 23.6521 * (log(temp)) + 
                 (118.67 / temp - 5.977 + 1.0495 * (log(temp))) * 
                 (sal^0.5) - 0.01615 * sal)
  K0 = exp(9345.17 / temp - 60.2409 + 23.3585 * (log(temp / 100)) + 
             sal * (0.023517 - 0.00023656 * temp + 
                      0.0047036 * ((temp / 100)^2)))
  
  ### Adjust for the effect of pressure (Millero 1995)
  delV1 = (-25.50) + 0.1271 * tempC
  delV2 = (-15.82) + (-0.0219 * tempC)
  delVspc= (-48.76) + (0.5304 * tempC)
  delVB = (-29.48) + 0.1622 * tempC + (2.608 / 1000) * 
    tempC^2
  delVw = (-25.60) + 0.2324 * tempC + (-3.6246 / 1000) * 
    tempC^2
  
  delk1 = (-3.08 / 1000) + (0.0877 / 1000) * tempC
  delk2 = (1.13 / 1000) + (-0.1475 / 1000) * tempC
  delkspc = (-11.76 / 1000) + (0.3692 / 1000) * tempC
  delkB = -2.84 / 1000
  delkw = (-5.13 / 1000) + (0.0794 / 1000) * tempC
  
  Ks1m = (exp(-((delV1 / (R * temp)) * press) + 
                ((0.5 * delk1) / (R * temp)) * press^2)) * Ks1m_st
  Ks2m = (exp(-((delV2 / (R*temp)) * press) + 
                ((0.5 * delk2) / (R*temp)) * press^2)) * Ks2m_st
  Ksspcm = (exp(-((delVspc / (R * temp)) * press) + 
                  ((0.5 * delkspc) / (R * temp)) * press^2)) * Ksspcm_st
  KsB = (exp(-((delVB / (R * temp)) * press) + 
               ((0.5 * delkB) / (R * temp)) * press^2)) * KsB_st
  Ksw = (exp(-((delVw / (R * temp)) * press) + 
               ((0.5 * delkw) / (R * temp)) * press^2)) * Ksw_st
  
  ### Correct for past seawater [Ca], [Mg], and [SO4] following ZT19
  Ks1 = Ks1m * (1 + (5e-3 * (xca / xcam - 1) + 
                       1.7e-2 * (xmg / xmgm - 1) + 
                       0.208 * (xso4 / xso4m - 1)))
  Ks2 = Ks2m* (1 + (0.157 * (xca / xcam - 1) + 
                      0.42 * (xmg / xmgm - 1) + 
                      0.176 * (xso4 / xso4m - 1)))
  Ksspc = Ksspcm * (1 + (0.185 * (xca / xcam - 1) + 
                           0.518 * (xmg / xmgm - 1) + 
                           0.106 * (xso4 / xso4m - 1)))
  
  ## Compute pH from DIC and CO2 ----
  fco2 = pco2 * 0.9968
  co2 = fco2 * K0
  hyd = (2 * Ks1 * Ks2) / (-Ks1 + sqrt(Ks1^2 - 4 * Ks1 * Ks2 * (1 - dic / co2)))
  pH = -(log(hyd) / log(10))
  
  ## Marine Proxies ----
  ### d11Bforam
  pKsB = -(log(KsB) / log(10))
  t1 = 10^(pKsB - pH)
  d11Bf = ((t1 * epsilon) - (t1 * d11Bsw) - 
             d11Bsw) / (-((t1 * alpha) + 1))
  
  ### d18Oforam (Bemis et al., 1998; Kim and O'Neil et al., 1997; Hollis et al., 2019)
  d18Osw.sc = d18Osw + (sw.sens * (sal - 35))
  d18Oswpdb = d18Osw.sc - 0.27
  d18Of.pr = d18Oswpdb + 
    ((4.64 - (((4.64^2) - (4 * 0.09 * (16.1 - tempC))) ^ 0.5)) / 
       (2 * 0.09))
  d18Of = d18Of.pr * (1 - indexop) + indexop * d18Oseccal
  
  ### d13Cforam
  d13Cf = d13Ca + d13Cepsilon
  
  ### Mg/Caforam following Hollis et al. (2019) Mg/Ca carb chem correction approach
  mgcasw = (xmg / xca)     
  Bcorr = ((mgcasw^Hp) / (mgcaswm^Hp)) * Bmod
  mgca_corr = Bcorr * (exp(A * tempC))
  mgcaf = mgca_corr / (1 - (sal - 35) * salcorrco)
  
  # Constants and derived values ----
  ## Derived values
  pCO2 = pco2 * 1e6 # atmospheric CO2 mixing ratio, ppm
  MAT = tempC + MAT_off # mean annual terrestrial site air temperature, C 
  d18.p = -15 + 0.58 * MAT # Precipitation d18O, ppt
  temp = tempC + 273.15
  PPCQ = MAP * PCQ_pf 
  TmPCQ = MAT + PCQ_to 
  
  ## Constants
  R = 83.131 # constant (cm^3 bar mol^-1 K^-1)
  alpha = 1.0272  # B isotope fractionation factor: Klochko et al. (2006)
  epsilon = (alpha - 1) * 1000  # Compute epsilon from alpha
  indexop = 0.7
  d18Oseccal = 0.85
  sw.sens = 0.558 # regression slope of GEOSECS obs. reported in Charles and Fairbanks (1990) after Duplessey et al. (1991)
  salcorrco = 0.042
  Hp = 0.74 # nonlinearity of the relationship b/w shell and Mg/Casw (Evans and Muller 2012, T. sacculifer)
  Bmod = 0.38 # modern pre-exponential constant in Mg/Ca-SST calibration (Anand et al., 2003; Evans et al., 2016) 
  A = 0.0757 # Exponential constant in Mg/Ca-SST calibration (Evans et al., 2016)
  
  ## Modern seawater chemistry
  xcam = 10.2821 # modern [Ca] (mmol kg^-1)
  xmgm = 52.8171 # modern [Mg] (mmol kg^-1)
  xso4m = 28.24  # modern [SO4] (mmol kg^-1)
  mgcaswm = xmgm/xcam # modern Mg/Ca of seawater
  
  ## Isotope ratio constants
  R13.VPDB = 0.011237
  R18.VSMOW = 0.0020052
  R18.VPDB = 0.0020672
  alpha18.diff = 1.028489
  a.theta = 0.05 # rate of increase of water content with depth (m-1) (Barnes and Allison, 1983)
  Rgas = 8.314462 # gas constant
  rho = 1000 # liquid water density (kg/m3)
  Dv.air = 2.44E-05 # water vapor diffusivity in air (m2/s) (Merlivat, 1978)
  
  # Priors ----
  ## Primary environmental ----
  tempC ~ dnorm(30, 1 / 3^2) # surface water temperature, C
  #  tempC = 30
  pco2 ~ dnorm(0.000875, 1 / 0.0001)T(0.0001, 0.002) # atmospheric CO2 mixing ratio
  #  pco2 = 0.000875
  MAT_off ~ dnorm(-18, 1 / 4^2) # offset between terrestrial and marine temperatures, C
  PCQ_to ~ dnorm(15, 1 / 2^2) # PCQ temperature offset, C
  MAP ~ dnorm(500, 1/50^2)T(100,) # mean annual terrestrial site precipitation, mm
  PCQ_pf ~ dbeta(0.5 / 0.9, 5) # PCQ precipitation fraction
  
  ## Secondary marine ----
  sal ~ dnorm(35, 1 / 2^2)T(25, 45) # surface water salinity, ppt
  press = 6 # pressure at habitation depth, bar
  xca ~ dnorm(21, 1 / 1^2)T(14, 28) # seawater [Ca], mmol/kg
  xmg ~ dnorm(68, 1 / 2^2)T(40, 90) # seawater [Mg], mmol/kg 
  xso4 ~ dnorm(14, 1 / 0.5^2)T(10, 18) # seawater [SO4], mmol/kg
  dic ~ dnorm(0.00205, 1 / 0.0001^2)T(0.0015, 0.0025) # seawater DIC, 
  d11Bsw ~ dnorm(38.45, 1 / 0.5^2) # seawater d11B, ppt
  d18Osw ~ dnorm(-1.2, 1 / 0.1^2) # seawater d18O, ppt
  d13Cepsilon ~ dnorm(10, 1 / 0.5 ^ 2) # offset between foram calcite and d13Catm
  
  ## Secondary soil ----
  tsc ~ dbeta(0.25 * 1000 / 0.75, 1000) # seasonal offset of PCQ for thermal diffusion
  ha ~ dbeta(0.35 * 500 / 0.65, 500) # PCQ atmospheric humidity
  f_R ~ dbeta(0.15 * 500 / 0.85, 500) # ratio of PCQ to mean annual respiration rate
  d13Ca ~ dnorm(-6.5, 1 / 0.5^2) # Atmospheric d13C, ppt
  ETR ~ dbeta(0.06 * 1000 / 0.94, 1000) # Soil evaporation / AET
  
  lat ~ dunif(32, 38) # terrestrial site latitude
  Ra = 42.608 - 0.3538 * abs(lat) # total radiation at the top of the atmosphere
  Rs = Ra * 0.16 * sqrt(12) # daily temperature range assumed to be 12
  L ~ dgamma(40, 1) # mean rooting depth, cm
  pore ~ dbeta(0.35 * 100 / 0.65, 100)T(0.06,) # soil porosity
  tort ~ dbeta(0.7 * 100 / 0.3, 100) # soil tortuosity
  
}
