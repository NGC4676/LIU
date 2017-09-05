# -*- coding: utf-8 -*-
"""
Created on Sun Jul 09 04:11:03 2017

@author: Qing Liu
"""
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.visualization import quantity_support
from scipy.integrate import quad

from smpy.ssp import BC
from smpy.smpy import CSP, LoadEAZYFilters, FilterSet, Observe
from smpy.sfh import exponential, delayed
from smpy.dust import Calzetti

#==============================================================================
# Composite Model
#==============================================================================
bc03 = BC('data/ssp/bc03/chab/lr/')

SFH_law = exponential
#==============================================================================
# Fang
#==============================================================================
Ages = np.logspace(-1., 1.3, 20)* u.Gyr
#SFH = np.logspace(-1, 1.3, 10) * u.Gyr
               
#==============================================================================
# Test 
#==============================================================================
#Ages = np.linspace(0.,13.5, 26)* u.Gyr
SFH = np.array([0.5, 1., 2., 3., 5., 10., 20.]) * u.Gyr

Dust = np.array([1.])   

models = CSP(bc03, age = Ages, sfh = SFH, 
             dust = Dust, metal_ind = 1.0, f_esc = 1.,
             sfh_law = SFH_law, dust_model = Calzetti)
#==============================================================================
# Renormalization
#==============================================================================
#use_bc = False
#
#if use_bc:
#    M_s = 10**lgMs.reshape((SFH.size, Ages.size)).T
#    SSFR = 10**lgSSFR.reshape((SFH.size, Ages.size)).T
#else:
#    M_s = np.array([quad(SFH_law, 0, T.value, args=(tau.value,))[0] \
#                    for T in Ages for tau in SFH]).reshape(Ages.size,SFH.size)
#
#    SSFR = np.squeeze(models.SFR.value)

#==============================================================================
# Mock Observation
#==============================================================================
eazy_library = LoadEAZYFilters('FILTER.RES.CANDELS')

print eazy_library.filternames

filters = FilterSet()
filters.addEAZYFilter(eazy_library, range(len(eazy_library.filternames)))  
    # FUV, NUV, U, B, V, R, I, J, H, K

synphot = Observe(models, filters, redshift=0.001)

mags = synphot.AB[0]
FUV_V = mags[0]-mags[4]
NUV_V = mags[1]-mags[4]
U_V = mags[2]-mags[4]
V_J = mags[4]-mags[7]
V_K = mags[4]-mags[9]
NUV_R = mags[1]-mags[5]
R_K = mags[5]-mags[9]

#==============================================================================
# UVJ
#==============================================================================
X, Y = V_J, U_V
plt.figure(figsize=(8,7))
with quantity_support():
    for k, Av in enumerate(Dust):
        for j, tau in enumerate(SFH):
            plt.plot(X[0,:,j,k,0], Y[0,:,j,k,0], label=r'$\tau = ${0:.1f}'.format(tau), 
                         lw=2, alpha=0.7)
            s = plt.scatter(X[0,:,j,k,0], Y[0,:,j,k,0], c=Ages/u.Gyr, 
                            s=30, cmap=plt.cm.viridis_r)
    col = plt.colorbar(s)
    col.set_label('Age(Gyr)')
    plt.plot([-1., 1.0, 1.6, 1.6], [1.3, 1.3, 2.01, 2.5], color='k', alpha=0.7)
    plt.legend(loc='best',fontsize=10,frameon=True,facecolor='w')
    plt.xlabel('V - J')
    plt.ylabel('U - V')
    plt.xlim([-0.75, 2.0])
    plt.ylim([-0.25, 2.25])
plt.show()

#plt.figure(figsize=(8,7))
#with quantity_support():
##    for k, Av in enumerate(Dust):
##            for j, tau in enumerate(SFH):
##                plt.plot(X[0,:,j,k,0], Y[0,:,j,k,0], label=r'$\tau = ${0:.2f}'.format(tau), 
##                         lw=2, alpha=0.7)
#    for i, Age in enumerate(Ages):
#        plt.plot(X[0,i,:,k,0], Y[0,i,:,k,0], label='$T = ${0:.2f}'.format(Age), 
#                 lw=2, ls='--',alpha=0.9)
#    plt.legend(loc=4,ncol=2,fontsize=10,frameon=True,facecolor='w')
#    plt.axhline(0.5,ls='--',color='k')
#    plt.xlabel('V - J')
#    plt.ylabel('U - V')
#    plt.xlim([-0.3, 1.0])
#    plt.ylim([-0.3, 1.0])
#plt.show()

#X, Y = V_J, U_V
#plt.figure(figsize=(8,7))
#with quantity_support():
#    for k, Av in enumerate(Dust):
#            for j, tau in enumerate(SFH):
#                plt.plot(X[0,:,j,k,0], Y[0,:,j,k,0], label=r'$\tau = ${0:.1f}'.format(tau), 
#                         c='grey',lw=2, alpha=0.7)
#    plt.plot([-1., 1.0, 1.6, 1.6], [1.3, 1.3, 2.01, 2.5], color='k', alpha=0.7)
#    plt.legend(loc='best',fontsize=10,frameon=True,facecolor='w')
#    plt.xlabel('V - J')
#    plt.ylabel('U - V')
#    plt.xlim([-0.75, 2.0])
#    plt.ylim([-0.25, 2.25])
#plt.show()
#==============================================================================
# Csed
#==============================================================================
theta = 34.8*u.deg
Ssed = pd.DataFrame(np.squeeze(np.sin(theta)*U_V + np.cos(theta)*V_J))
Csed = pd.DataFrame(np.squeeze(np.cos(theta)*U_V - np.sin(theta)*V_J))
UV = pd.DataFrame(np.squeeze(U_V))

plt.figure(figsize=(7,7))
for i, tau in enumerate(SFH):
    plt.plot(Ages,UV[i],label=r'$\tau = ${0:.1f}'.format(tau))
plt.legend(loc='best',fontsize=12,frameon=True,facecolor='w')
plt.xlabel('T (Gyr)',fontsize=15)
plt.ylabel('U - V',fontsize=15)
plt.axhline(0.5,color='k',ls='--')
plt.axvline(2.,color='b',ls='--',lw=1,alpha=0.3)
plt.axvline(6.,color='b',ls='--',lw=1,alpha=0.3)
plt.xlim(-1.,13.)
plt.ylim(-0.1,1.2)
plt.show()
