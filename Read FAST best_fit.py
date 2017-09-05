# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 11:17:39 2017

@author: Qing Liu
"""

import numpy as np
import asciitable as asc
import seaborn as sns
import matplotlib.pyplot as plt

Ages = 10**(np.linspace(7.0,10.0,13)-9)
SFH = 10**(np.array([7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])-9)

for i in range(Ages.size):
    for j in range(SFH.size):
        plt.figure(figsize=(8,4))
        table  = asc.read('BEST_FITS/grid_exp_%d.fit'
                          %(i*SFH.size+j+1))
        wl, fl = table['col1'], table['col2']
        table_i  = asc.read('BEST_FITS/grid_exp_%d.input_res.fit'
                            %(i*SFH.size+j+1))
        wave, flux = table_i['col1'], table_i['col2']

        plt.plot(wl,fl,'k',alpha=0.7)
        plt.plot(wave,flux/1e29,'rs',alpha=0.9)
        plt.title(r'T = %.2f Gyr    $\tau$ = %.2f Gyr'
                  %(Ages[i],SFH[j]),fontsize=15)
        plt.xlabel('wavelength ($\AA$)',fontsize=12)
        plt.ylabel('flux ($10^{-19} erg s^{-1} cm^{-2} \AA^{-1}$)',fontsize=12)
        plt.xlim(0,25000)
        plt.savefig('SED-fitting/exp_tau%dT%d'%(j+1,i+1))
        plt.show()
        print i,j