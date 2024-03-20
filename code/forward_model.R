# System state
tempC = 28
sal = 35
press = 6
R = 83.131 # constant (cm^3 bar mol^-1 K^-1)
xca = 21
xmg = 68
xso4 = 14
dic = 0.00205
pco2 = 0.000875

# Set modern concentrations for Mg, Ca, and SO4
xcam = 10.2821 # modern [Ca] (mmol kg^-1)
xmgm = 52.8171 # modern [Mg] (mmol kg^-1)
xso4m = 28.24  # modern [SO4] (mmol kg^-1)
mgcaswm = xmgm/xcam # modern Mg/Ca of seawater

# Carbonate system ---- 

# Calculate equil. constants using salinity and temp:
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

# Adjust equil. constants for the effect of pressure (Millero 1995):
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

# K*1, K*2, and K*spc are corrected for past seawater [Ca], [Mg], and [SO4] following ZT19
Ks1 = Ks1m * (1 + (5e-3 * (xca / xcam - 1) + 
                           1.7e-2 * (xmg / xmgm - 1) + 
                           0.208 * (xso4 / xso4m - 1)))
Ks2 = Ks2m* (1 + (0.157 * (xca / xcam - 1) + 
                          0.42 * (xmg / xmgm - 1) + 
                          0.176 * (xso4 / xso4m - 1)))
Ksspc = Ksspcm * (1 + (0.185 * (xca / xcam - 1) + 
                               0.518 * (xmg / xmgm - 1) + 
                               0.106 * (xso4 / xso4m - 1)))

# Compute pH from DIC and CO2
fco2 = pco2 * 0.9968
co2 = fco2 * K0
hyd = (2 * Ks1 * Ks2) / (-Ks1 + sqrt(Ks1^2 - 4 * Ks1 * Ks2 * (1 - dic / co2)))
pH = -log(hyd, 10)

# Next thing ----

# Compute d11Bforam from pH and d11Bsw
pKsB = -(log(KsB) / log(10))
t1 = 10^(pKsB - pH)
d11Bb = ((t1 * epsilon) - (t1 * d11Bsw) - 
              d11Bsw) / (-((t1 * alpha) + 1))