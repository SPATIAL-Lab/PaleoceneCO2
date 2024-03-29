
## System state
tempC = 30
sal = 35
press = 6
xca = 21
xmg = 68
xso4 = 14
dic = 0.00205
pco2 = 0.000875
pCO2 = pco2 * 1e6
d11Bsw = 38.45
d18Osw = -1.2
MAT = 12
PCQ_to = 15
MAP = 500
PCQ_pf = 0.1
pore = 0.35
tort = 0.7
tsc = 0.25
lat = 35
ha = 0.35
T.diff = 6
L = 40 # mean rooting depth, cm
f_R = 0.15
d13Ca = -6.5

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

# isotope ratio constants
R13.VPDB = 0.011237
R18.VSMOW = 0.0020052
R18.VPDB = 0.0020672
alpha18.diff = 1.028489
a.theta = 0.05 # rate of increase of water content with depth (m-1) (Barnes and Allison, 1983)
Rgas = 8.314462 # gas constant
rho = 1000 # liquid water density (kg/m3)
Dv.air = 2.44E-05 # water vapor diffusivity in air (m2/s) (Merlivat, 1978)


# Carbonate system ---- 

## Calculate equil. constants using salinity and temp:
temp = tempC + 273.15
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

## Adjust equil. constants for the effect of pressure (Millero 1995):
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

## K*1, K*2, and K*spc are corrected for past seawater [Ca], [Mg], and [SO4] following ZT19
Ks1 = Ks1m * (1 + (5e-3 * (xca / xcam - 1) + 
                           1.7e-2 * (xmg / xmgm - 1) + 
                           0.208 * (xso4 / xso4m - 1)))
Ks2 = Ks2m* (1 + (0.157 * (xca / xcam - 1) + 
                          0.42 * (xmg / xmgm - 1) + 
                          0.176 * (xso4 / xso4m - 1)))
Ksspc = Ksspcm * (1 + (0.185 * (xca / xcam - 1) + 
                               0.518 * (xmg / xmgm - 1) + 
                               0.106 * (xso4 / xso4m - 1)))

## Compute pH from DIC and CO2
fco2 = pco2 * 0.9968
co2 = fco2 * K0
hyd = (2 * Ks1 * Ks2) / (-Ks1 + sqrt(Ks1^2 - 4 * Ks1 * Ks2 * (1 - dic / co2)))
pH = -log(hyd, 10)

# Marine Proxies ----
## d11Bforam
pKsB = -(log(KsB) / log(10))
t1 = 10^(pKsB - pH)
d11Bb = ((t1 * epsilon) - (t1 * d11Bsw) - 
              d11Bsw) / (-((t1 * alpha) + 1))

## d18Oforam (Bemis et al., 1998; Kim and O'Neil et al., 1997; Hollis et al., 2019)
d18Osw.sc = d18Osw + (sw.sens * (sal - 35))
d18Oswpdb = d18Osw.sc - 0.27
d18Of.pr = d18Oswpdb + 
  ((4.64 - (((4.64^2) - (4 * 0.09 * (16.1 - tempC))) ^ 0.5)) / 
     (2 * 0.09))
d18Of = d18Of.pr * (1 - indexop) + indexop * d18Oseccal

## Mg/Caforam following Hollis et al. (2019) Mg/Ca carb chem correction approach
mgcasw = (xmg / xca)     
Bcorr = ((mgcasw^Hp) / (mgcaswm^Hp)) * Bmod
mgca_corr = Bcorr * (exp(A * tempC))
mgcaf = mgca_corr / (1 - (sal - 35) * salcorrco)

# Soil carbonate ----
## PCQ precipitation and temperature
PPCQ = MAP * PCQ_pf 
TmPCQ = MAT + PCQ_to 

## Depth to carbonate formation based on Retallack (2005) data, meters
z = (0.093 * MAP + 13.12)
z_m = z / 100

## Soil temperatures at depth z
d = sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
Tsoil = MAT + (PCQ_to * sin(2 * 3.1415 * tsc - z / d)) / exp(z / d) 
Tsoil.K = Tsoil + 273.15

## Potential Evapotranspiration - Hargreaves and Samani (1982) and Turc (1961)
Ra = 42.608 - 0.3538 * abs(lat) # total radiation at the top of the atmosphere
Rs = Ra * 0.16 * sqrt(12) # daily temperature range assumed to be 12
PET_A_D = ifelse(ha < 0.5, 0.013 * (MAT / (MAT + 15)) * (23.885 * Rs + 50) * (1 + ((0.5 - ha) / 0.7)),
                      0.013 * (MAT / (MAT + 15)) * (23.885 * Rs + 50))
PET_A_D = max(PET_A_D, 0.01)
PET_A_A = PET_A_D * 365

## PET_PCQ
Tair_PCQ = MAT + PCQ_to
PET_PCQ_D = ifelse(ha < 0.5, 0.013 * (Tair_PCQ / (Tair_PCQ + 15)) * (23.885 * Rs + 50) * (1 + ((0.5 - ha) / 0.7)),
                        0.013 * (Tair_PCQ / (Tair_PCQ + 15)) * (23.885 * Rs + 50))
PET_PCQ_D = max(PET_PCQ_D, 0.01)
PET_PCQ = PET_PCQ_D * 90

## AET - actual evapotranspiration
# AET in mm/quarter from Budyko curve - Pike (1964)
AET_PCQ = PPCQ * (1 / (sqrt(1 + (1 / ((PET_PCQ / (PPCQ)))) ^ 2)))

### water vapor pressure and relative humidity - dew point temperature - Qiu et al. (2021)
#AI = min(PPCQ / PET_PCQ, 1)
#a.t = -2.81 * log(AI) + 1.28
#Tmin = Tair_PCQ - T.diff
#Tdew = Tmin - a.t
#e_a = 0.6108 * exp(17.27 * Tdew / (Tdew + 237.3)) * 10^3 # actual water vapor pressure (kPa)
#e_s = 0.6112 * exp(17.67 * Tair_PCQ / (Tair_PCQ + 243.5)) * 10^3 # saturated water vapor pressure (kPa)
#ha = e_a / e_s

### ------------------------------------- CARBON ISOTOPE SYSTEM ------------------------------------------------
## Free air porosity
# Scales volumetrically w/ excess precipitation relative to pore space, assume a minimum of 5% volumetric water content
FAP1 = min((pore - ((PPCQ - AET_PCQ) / (L * 10 * pore))), pore - 0.05)
# At least 1% free air porosity
FAP = max(FAP1, 0.01)

### Respiration
# characteristic production depth (cm) - Quade (2007)
k = L / 2 / log(2)

## calculate soil respiration rate 
R_PCQ_D_m1 = 1.25 * exp(0.05452 * TmPCQ) * PPCQ / (127.77 + PPCQ)
# the proportion of respired CO2 during PCQ - Fischer-Femal & Bowen (2020)
R_PCQ_D = R_PCQ_D_m1 * f_R # (gC/m2/d)

# convert to molC/cm3/s
R_PCQ_D1 = R_PCQ_D / (12.01 * 100^2)  # from gC/m2/d to molC/cm2/d
R_PCQ_S = R_PCQ_D1 / (24 * 3600)  # molC/ cm2 / s
R_PCQ_S_0 = R_PCQ_S / (L * pore) # Quade et al. (2007)

# soil CO2 diffusion rate
Dair = 0.1369 * (Tsoil.K / 273.15) ^ 1.958
DIFC = FAP * tort * Dair

# S(z) 
S_z_mol = k ^ 2 * R_PCQ_S_0 / DIFC * (1 - exp(-z / k)) # (mol/cm3)
S_z = S_z_mol * (0.08206 * Tsoil.K * 10^9) # ppmv

## estimate the d13Cr of soil-respired CO2
DD13_water = 25.09 - 1.2 * (MAP + 975) / (27.2 + 0.04 * (MAP + 975))
D13C_plant = (28.26 * 0.22 * (pCO2 + 23.9)) / (28.26 + 0.22 * (pCO2 + 23.9)) - DD13_water # schubert & Jahren (2015)
d13Cr = d13Ca - D13C_plant

### d13C of pedogenic carbonate
d13Cs = (pCO2 * d13Ca + S_z * (1.0044 * d13Cr + 4.4))/(S_z + pCO2)
# T-dependent fractionation equation - Romanek et al., 1992
d13Cc = ((1 + (11.98 - 0.12 * Tsoil) / 1000) * (d13Cs + 1000)) - 1000

### ------------------------------------- OXYGEN ISOTOPE SYSTEM ------------------------------------------------

#### rainfall as reservoir water
d18.p_m = -0.0082 * MAP - 3.639  # CLP observations - Da (2023)
# d18p_m = -7 # assign a fixed value
d18.p ~ dnorm(d18.p_m, 1 / 3 ^ 2)
R18.p = (d18.p / 1000 + 1) * R18.VSMOW
dp18.p = log(d18.p/1000 + 1) * 1000
dp17.p = 0.5267 * dp18.p + 0.013 # Global meteoric water - Aron (2021)
R17.p = exp(dp17.p / 1000) * R17.VSMOW

# equilibrium fractionation factor for vapor to liquid - (Horita and Wesolowski 1994)
alpha18.eq = 1 / exp(((1.137e6 / (Tsoil.K ^ 2) - 0.4156e3/Tsoil.K - 2.0667) /1000))
alpha17.eq = alpha18.eq ^ theta.eq

# water vapor in air assumed to be in equilibrium with reservoir water
R18.a = R18.p * alpha18.eq
R17.a = R17.p * alpha17.eq

### soil water evaporation
# evaporation is 6+/-4% of AET - Good (2015)
e_mean = 0.06
e_var = 0.04^2
e_size = e_mean*(1-e_mean)/e_var - 1
e_alpha = e_mean * e_size
e_beta = (1-e_mean) * e_size
ETR ~ dbeta(e_alpha, e_beta)
E1 = ETR * AET_PCQ
E = max(E1, 1) # minimum of 1 mm
E_s = E / (1000 * 90 * 24 * 3600) # soil evaporation rate in m/sec

### water vapor diffusivity
Dv.soil = Dv.air * tort * (pore - 0.05) # effective diffusivity of water vapor in soil (m2/s)
es = (0.611 * exp(17.502 * Tsoil / (Tsoil + 240.97))) * 1000 # saturated water vapor pressure from Tetens formula
N.sat = 0.01802 * es / (Rgas * Tsoil.K) # saturated water vapor concentration at a given temperature
z.bar = N.sat * Dv.soil / (E_s * rho) # penetration depth (m)
z.ef1 = (1 - ha) * z.bar # the thickness of the water vapor phase region (m)

# liquid water diffusivity (m2/s) (Easteal 1984)
Dl = exp(1.6766 + 1.6817 * (1000 / Tsoil.K) - 0.5773 * (1000 / Tsoil.K)^2) * 10^-9 
Dl.soil = Dl * pore * tort # effective diffusivity of liquid water (m2/s)
z.hat = Dl.soil / E_s # the decay length (mean penetration depth)

# the evaporation front
h.ef = ha + z.ef / z.bar # humidity at the evaporation front
R18.ef = (alpha18.diff * R18.p * (z.ef / z.bar) + ha * R18.a) / (h.ef * alpha18.eq) # isotopic composition at the evaporation front
R17.ef = (alpha17.diff * R17.p * (z.ef / z.bar) + ha * R17.a) / (h.ef * alpha17.eq)
