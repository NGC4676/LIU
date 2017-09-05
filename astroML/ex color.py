import numpy as np
from matplotlib import pyplot as plt
from astroML.datasets import fetch_imaging_sample

data = fetch_imaging_sample()
print data.dtype.names
def get_stars_and_galaxies(Nstars=33000, Ngals=33000):
    """Get the subset of star/galaxy data to plot"""
    data = fetch_imaging_sample()

    objtype = data['type']

    stars = data[objtype == 6][:Nstars]
    galaxies = data[objtype == 3][:Ngals]

    return stars, galaxies
def plot_stars_and_galaxies(stars, galaxies):
	plt.figure(figsize=(10, 7.5))
	plt.plot(galaxies['uRaw']+galaxies['gRaw']+galaxies['rRaw'],
				galaxies['gRaw']-galaxies['rRaw'],
				linestyle='none',marker=',')
stars,galaxies = get_stars_and_galaxies()
plot_stars_and_galaxies(stars, galaxies)
plt.title('Galaxies',fontsize=15)
plt.xlabel('u+g+r',fontsize=15)
plt.ylabel('g-r',fontsize=15)
plt.show()

