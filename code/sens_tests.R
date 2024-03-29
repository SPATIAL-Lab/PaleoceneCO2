

vars = ctrl()
baseCase = fm(vars)

# Sensitivity to global temperature
vars = ctrl()
vars$tempC = (20:35)
sens.t = fm(vars)
sens.t$tempC = vars$tempC

# Sensitivity to CO2
vars = ctrl()
vars$pco2 = seq(250e-6, 2000e-6, by = 50e-6)
sens.co2 = fm(vars)
sens.co2$pco2 = vars$pco2

# Sensitivity to local (terrestrial) temperature
vars = ctrl()
vars$MAT_off = (-25:-5)
sens.MAT = fm(vars)
sens.MAT$MAT_off = vars$MAT_off

plot(sens.t)
plot(sens.co2)
plot(sens.MAT)

# Sensitivity with climate sensitivity
vars = ctrl()
vars$pco2 = seq(250e-6, 2000e-6, by = 50e-6)
vars$tempC = 30 + 4 * log(vars$pco2 / 875e-6) / log(2) # ESS = 4
vars$MAT_off = -18 + 2 * log(vars$pco2 / 875e-6) / log(2) # continental amplification
sens.sens = fm(vars)
sens.sens$pco2 = vars$pco2

plot(sens.sens)
