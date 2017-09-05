import numpy as np
from matplotlib import pyplot as plt
from astroML.datasets import fetch_sdss_specgals

data = fetch_sdss_specgals()[:]
print data.shape
print data.dtype.names
lgm = data['lgm_tot_p50']
uv_r = data['modelMag_u']-data['modelMag_r']

H, xbins, ybins = np.histogram2d(lgm,uv_r,
                                 bins=(np.linspace(8, 12, 200),
                                       np.linspace(1, 5, 200)))																															
cmap = 'jet'
fig, ax = plt.subplots(figsize=(10, 7.5))																																			
ax.imshow(np.log10(H).T, origin='lower',
          extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]],
          cmap=cmap, interpolation='nearest',
          aspect='auto')

plt.ylim(0, 6)
plt.xlim(6, 15)
plt.title('Galaxies',fontsize=15)
plt.xlabel('lg M',fontsize=15)
plt.ylabel('UV-r',fontsize=15)
plt.show()
