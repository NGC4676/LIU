from matplotlib import pyplot as plt
import xlrd  
import numpy as np
from astropy.modeling import models, fitting
#read and plot spectrum
data = xlrd.open_workbook(r"C:\Users\mac\Desktop\NAOC\data\Mrk-231.xlsx")   
table=data.sheets()[0]
wave=np.array(table.col_values(0)[:8163])
spec=np.array(table.col_values(1)[:8163])

#fit Gaussion continuum
gg_init = models.Gaussian1D(0.5e-14, 4600,800)+models.Gaussian1D(0.8e-14, 7500, 3000)
fitter = fitting.SLSQPLSQFitter()
gg_fit = fitter(gg_init, wave, spec)

# Plot the data with the best-fit model
plt.figure(figsize=(10,3))
plt.plot(wave, spec, wave, gg_fit(wave))
plt.xlabel('Wavelength',fontsize=15)
plt.ylabel('Flux',fontsize=15)
plt.yscale('log')

#Subtract continuum
plt.figure(figsize=(10,3))
line=spec-gg_fit(wave)
plt.plot(wave, line)
plt.xlabel('Wavelength',fontsize=15)
plt.ylabel('Flux',fontsize=15)

#Fit Ha emssion
Mu=6558
Gamma=45
B=2.5e-14
line2=B*Gamma**2/(Gamma**2 + (wave[6800:7400] - Mu)**2)
plt.plot(wave[6800:7400], line2,'r--')