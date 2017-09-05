import numpy
from matplotlib import pyplot as plt
from astroML.datasets import fetch_sdss_specgals

data = fetch_sdss_specgals()[:50000]
print data.dtype.names,data.shape

def plot_galaxies(galaxies):
	plt.figure(figsize=(10, 7.5))
	plt.plot(galaxies['lgm_tot_p50'],
				galaxies['modelMag_u']-galaxies['modelMag_r'],
				linestyle='none',marker=',')
	plt.ylim(0, 6)
	plt.xlim(5, 15)
	plt.title('Galaxies',fontsize=15)
	plt.xlabel('lg M',fontsize=15)
	plt.ylabel('UV-r',fontsize=15)
				
plot_galaxies(data)


plt.show()
