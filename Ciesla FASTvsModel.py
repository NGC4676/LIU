# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 11:28:09 2017

@author: Qing Liu
"""

import numpy as np
import asciitable as asc
import seaborn as sns
import matplotlib.pyplot as plt

#==============================================================================
# Read model result 
#==============================================================================
lgAge = np.log10(np.linspace(0.5,13.,26)* 1e9)
M_seed = np.linspace(6.0, 10.0, 5)

BC03_out = asc.read('Ciesla/CSP result.txt')
lgMs = BC03_out['lgM*'].reshape((M_seed.size, lgAge.size)).T.ravel()
lgSSFR = BC03_out['lgsSFR'].reshape((M_seed.size, lgAge.size)).T.ravel()
                 
lgAges = np.array([lgAge for m in M_seed]).T.ravel()
M_class =  np.array([M_seed for T in range(lgAge.size)]).T.ravel()
#==============================================================================
# Read FAST SED-fitting result
#==============================================================================
table = asc.read('Ciesla/Ciesla_del_dust.fout')
Ages_FAST = table.lage
SFH_FAST = table.ltau
M_FAST = table.lmass
SSFR_FAST = table.lssfr
Av_FAST = table.Av

#==============================================================================
# Plot FAST vs Model
#==============================================================================
#fig,axes = plt.subplots(figsize=(15,6), nrows=1, ncols=3)
#xx=np.linspace(-20,20,20)
#legends = ['log(T/Gyr)', r'log$(\rm{M_*}/\rm{M_{\odot}})$', r'log(sSFR/yr$^{-1}$)']
#for i, (l, model, FAST) in enumerate(zip(legends,[lgAges,lgMs,lgSSFR],
#                                                 [Ages_FAST,M_FAST,SSFR_FAST])):
#    ax = plt.subplot(1,3,i+1) 
#    plt.plot(xx,xx,'k--',alpha=0.5)
#    s = plt.scatter(model, FAST, c=lgAges,
#                    cmap='jet', label=l, s=20, alpha=0.7)
#    plt.xlabel('Model',fontsize=12)
#    plt.ylabel('FAST',fontsize=12)
#    if i==2:
#        plt.xlim(-14,-6)
#        plt.ylim(-14,-6)
#    elif i==1:
#        plt.xlim(0.,14.)
#        plt.ylim(0.,14.)
#    else:
#        plt.xlim(7.,11.)
#        plt.ylim(7.,11.)
#        plt.xlabel('Time since Big Bang',fontsize=12)
#        plt.ylabel('Age return by FAST',fontsize=12)
#    plt.legend(loc=4, fontsize=12, frameon=True, facecolor='w')
#fig.subplots_adjust(bottom=0.25)
#cbar_ax = fig.add_axes([0.2, 0.12, 0.6, 0.03])
#colorbar = fig.colorbar(s, orientation='horizontal',cax=cbar_ax)
#colorbar.set_label('log (Age/yr)')
#plt.suptitle('Aldo SFH Model vs FAST Fitting using Tau Model',fontsize=15,y=0.95)
#plt.show() 

#==============================================================================
# Dust
#==============================================================================
fig,axes = plt.subplots(figsize=(10,10), nrows=2, ncols=2)
xx=np.linspace(-20,20,20)
legends = ['log(T/Gyr)', r'log$(\rm{M_*}/\rm{M_{\odot}})$', r'log(sSFR/yr$^{-1}$)']
for i, (l, model, FAST) in enumerate(zip(legends,[lgAges,lgMs,lgSSFR],
                                                 [Ages_FAST,M_FAST,SSFR_FAST])):
    ax = plt.subplot(2,2,i+1) 
    plt.plot(xx,xx,'k--',alpha=0.5)
    s = plt.scatter(model, FAST, c=lgAges,
                    cmap='jet', label=l, s=20, alpha=0.7)
    plt.xlabel('Model',fontsize=12)
    plt.ylabel('FAST',fontsize=12)
    if i==2:
        plt.xlim(-14,-6)
        plt.ylim(-14,-6)
    elif i==1:
        plt.xlim(0.,14.)
        plt.ylim(0.,14.)
    else:
        plt.xlim(7.,11.)
        plt.ylim(7.,11.)
        plt.xlabel('Age of Model',fontsize=12)
        plt.ylabel('Age return by FAST',fontsize=12)
    plt.legend(loc=4, fontsize=12, frameon=True, facecolor='w')
fig.subplots_adjust(bottom=0.2)
cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.03])
colorbar = fig.colorbar(s, orientation='horizontal',cax=cbar_ax)
colorbar.set_label('log (Age/yr)')
ax = plt.subplot(2,2,4)
plt.hist(Av_FAST,alpha=0.7,label='Av')
plt.axvline(1.0,color='k',ls='--')
plt.xlabel('FAST',fontsize=12)
plt.legend(loc='best', fontsize=12, frameon=True, facecolor='w')
plt.suptitle('MS-SFH Models vs FAST Fitting with Delayed Templates',fontsize=15,y=0.92)
plt.show()