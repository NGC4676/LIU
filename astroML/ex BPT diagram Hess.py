import numpy as np
from matplotlib import pyplot as plt
from astroML.datasets import fetch_sdss_specgals

data = fetch_sdss_specgals()[:]
print data.shape
print data.dtype.names
y = np.log(data['oiii_5007_flux']/data['h_beta_flux'])
x = np.log(data['nii_6584_flux']/data['h_alpha_flux'])
H, xbins, ybins = np.histogram2d(x,y,bins=(np.linspace(-5, 5, 1000),
															np.linspace(-5, 5, 1000)))																															
cmap = 'jet'
fig, ax = plt.subplots(figsize=(10, 7.5))	

#plt.plot(x,y,linestyle='none',marker=',')
																																		
ax.imshow(np.log10(H).T, origin='lower',
          extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]],
          cmap=cmap, interpolation='nearest',
          aspect='auto')

plt.ylim(-5,5)
plt.xlim(-5,5)
plt.title('Galaxies',fontsize=15)
plt.xlabel('lg nii_6584_flux/h_alpha_flux',fontsize=15)
plt.ylabel('lg oiii_5007_flux/h_beta_flux',fontsize=15)

plt.show()
