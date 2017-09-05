# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 18:33:10 2017

@author: Qing Liu
"""

import numpy as np
import seaborn as sns
import asciitable as asc
import pandas as pd
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.table import Table, Column
from astropy.visualization import quantity_support
from astropy.cosmology import FlatLambdaCDM, z_at_value

from smpy.ssp import BC
from smpy.smpy import CSP, LoadEAZYFilters, FilterSet, Observe
from smpy.sfh import gauss_lorentz_hermite, glh_prtb, glh_absprtb, glh_ST
from smpy.dust import Calzetti

# Go to smpy.smpy to set cosmo, the current is below:
cosmo = FlatLambdaCDM(H0=67.8, Om0=0.307) 

#==============================================================================
# Composite Model
#==============================================================================
bc03 = BC('data/ssp/bc03/chab/hr/')
Ages = np.linspace(1.,8.5,76)* u.Gyr
M_today = np.array([9.0,9.5,10.0,10.5,11.0])
#==============================================================================
# Test 
#==============================================================================
Ages = np.linspace(2.,8.,61)* u.Gyr
M_today = np.linspace(9.5,11.5,5)
                  
perturb = True
if perturb:
    A, P = 0.3, 1.
    SFH_law = glh_ST
    SFH_fit = asc.read('Aldo/Aldo SFH prtb.txt')
    SFH_pars = zip(SFH_fit.c1, SFH_fit.mu, SFH_fit.sigma,\
               SFH_fit.h13, SFH_fit.h14,\
               SFH_fit.c2, SFH_fit.x0, SFH_fit.gama,\
               SFH_fit.h23, SFH_fit.h24,\
               A*np.ones_like(M_today), P*np.ones_like(M_today))        
else:
    SFH_law = gauss_lorentz_hermite
    SFH_fit = asc.read('Aldo/Aldo SFH params.txt')
    SFH_pars = zip(SFH_fit.c1, SFH_fit.mu, SFH_fit.sigma, SFH_fit.h13, SFH_fit.h14,\
               SFH_fit.c2, SFH_fit.x0, SFH_fit.gama, SFH_fit.h23, SFH_fit.h24)
    

Dust = np.array([0.0])   

models = CSP(bc03, age = Ages, sfh = SFH_pars, 
             dust = Dust, metal_ind = 1., f_esc = 1.,
             sfh_law = SFH_law, dust_model = Calzetti,
             neb_cont=True, neb_met=False)

SED = np.squeeze(models.SED.value)

#iT = 25 
#plt.figure(figsize=(8,4))
#for j in range(6):
#    plt.semilogy(models.wave,SED[iT,j,:],lw=1)
#plt.fill_between(np.linspace(3e3,4.2e3,10), 1, 1e-7,facecolor='grey', alpha=0.3)
#plt.fill_between(np.linspace(4.8e3,6.9e3,10), 1, 1e-7,facecolor='grey', alpha=0.3)
#plt.fill_between(np.linspace(1.1e4,1.34e4,10), 1, 1e-7,facecolor='grey', alpha=0.3)
#plt.xlabel('Wavelength',fontsize=15)
#plt.ylabel('Flux',fontsize=15)
#plt.text(1.2e4,1e-3,'T = %.1f Gyr'%Ages[iT].value,fontsize=15)
#plt.xlim(1e3,1.5e4)
#plt.ylim(3e-6,5e-3)

#for k in range(220):
#    print (models.ta[k+1]-models.ta[k]).to('Gyr'),models.ta[k].to('Gyr')

#==============================================================================
# Renormalization
#==============================================================================
use_bc = True
if use_bc:
    if perturb:
        BC03_out = asc.read('Aldo/ST/m%.1f CSP result.txt'%A)
    else:
        BC03_out = asc.read('Aldo/CSP result.txt')
    lgMs = BC03_out['lgM*'].reshape((M_today.size, Ages.size)).T
    lgSSFR = BC03_out['lgsSFR'].reshape((M_today.size, Ages.size)).T            
    Ms = 10**lgMs
    SSFR = 10**lgSSFR

#==============================================================================
# Mock Observation
#==============================================================================
eazy_library = LoadEAZYFilters('FILTER.RES.CANDELS')

print eazy_library.filternames

filters = FilterSet()
filters.addEAZYFilter(eazy_library, range(len(eazy_library.filternames)))  
Names = ['FUV', 'NUV', 'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
Redshifts = np.arange(0.5,2.5,0.5)
synphot = Observe(models, filters, redshift=0.001)

mags = synphot.AB[0]
FUV, NUV, U, B, V, R, I, J, H, K = mags
FUV_V = FUV - V
NUV_V = NUV - V
U_V = U - V
V_J = V - J
V_K = V - K

#==============================================================================
# UVJ Tracks
#==============================================================================
X, Y = V_J, U_V
plt.figure(figsize=(8,7))
with quantity_support():
    for k, Av in enumerate(Dust):
            for j,m in enumerate(M_today):
                plt.plot(X[0,:,j,k,0], Y[0,:,j,k,0], label='log M$_{today}$ = %.1f'%m,
                         lw=2, alpha=0.7)
                s = plt.scatter(X[0,:,j,k,0], Y[0,:,j,k,0], c=Ages/u.Gyr, 
                                s=15, cmap=plt.cm.viridis_r)
    col = plt.colorbar(s)
    col.set_label('Age(Gyr)')
    plt.plot([-1., 1.0, 1.6, 1.6], [1.3, 1.3, 2.01, 2.5], color='k', alpha=0.7)
    plt.legend(loc='best',fontsize=10,frameon=True,facecolor='w')
    plt.xlabel('V - J',fontsize=12)
    plt.ylabel('U - V',fontsize=12)
    plt.xlim([-0.75, 2.0])
    plt.ylim([-0.25, 2.25])
plt.show()

#==============================================================================
# Single Perturbed UVJ color code in Av
#==============================================================================
#X, Y = V_J, U_V
#plt.figure(figsize=(9,6))
#with quantity_support():
#    for k, Av in enumerate(Dust):
#            for j,m in enumerate(M_today):
#                plt.subplot(2,3,j+1)
#                plt.plot(X[0,:,j,k,0], Y[0,:,j,k,0], label='log M$_{today}$ = %.1f'%m,
#                         lw=2, alpha=0.7)
##                s = plt.scatter(X[0,:,j,k,0], Y[0,:,j,k,0], c=Av_fit[j], 
##                                 s=30, cmap='rainbow_r')
#                
##                col = plt.colorbar(s)
##                col.set_label('$Av_{FAST}$')
#                j+=1
#                plt.plot([-1., 1.0, 1.6, 1.6], [1.3, 1.3, 2.01, 2.5], color='k', alpha=0.7)
#                plt.legend(loc='best',fontsize=10,frameon=True,facecolor='w')
#                plt.xlabel(' ',fontsize=12)
#                plt.ylabel(' ',fontsize=12)
#                plt.xlim([-0.75, 2.0])
#                plt.ylim([-0.25, 2.25])
#                
#plt.tight_layout()
#plt.show()

#==============================================================================
# Csed/UV vs Time
#==============================================================================
theta = 34.8*u.deg
Ssed = pd.DataFrame(np.squeeze(np.sin(theta)*U_V + np.cos(theta)*V_J))
Csed = pd.DataFrame(np.squeeze(np.cos(theta)*U_V - np.sin(theta)*V_J))
UV = pd.DataFrame(np.squeeze(U_V))

plt.figure(figsize=(7,7))
for i, m in enumerate(M_today):
    plt.plot(Ages,UV[i],label='log M$_{today}$ = %.1f'%m)
plt.legend(loc='best',fontsize=12,frameon=True,facecolor='w')
plt.xlabel('T (Gyr)',fontsize=15)
plt.ylabel('U - V',fontsize=15)
plt.axhline(0.5,color='k',ls='--')
plt.axvline(3.,color='b',ls='--',lw=1,alpha=0.3)
plt.axvline(7.,color='b',ls='--',lw=1,alpha=0.3)
plt.xlim(0.,10.)
plt.ylim(-0.1,2.)
plt.show()

#plt.figure(figsize=(7,7))
#for i, m in enumerate(M_today):
#    plt.plot(Ages,Csed[i],label='log M$_{today}$ = %.1f'%m)
#plt.legend(loc='best',fontsize=12,frameon=True,facecolor='w')
#plt.xlabel('T (Gyr)',fontsize=15)
#plt.ylabel('C$_{SED}$',fontsize=15)
#plt.axhline(0.25,color='k',ls='--')
#plt.axvline(3.,color='b',ls='--',lw=1,alpha=0.3)
#plt.axvline(7.,color='b',ls='--',lw=1,alpha=0.3)
#plt.xlim(-1.,13.)
#plt.ylim(-0.1,1.2)
#plt.show()

#==============================================================================
# Make Table
#==============================================================================
make_table = False
multiple = True
real = True

if real:
    table_uds = Table.read('gds_merged_v1.1.fits')
    data = table_uds.to_pandas()
    
    z_cat = np.arange(0.5,3.0,0.5)
    M_cat = np.arange(9.0,11.5,0.25)
    NS_H = np.zeros((z_cat.size, M_cat.size))
    NS_mock = np.zeros_like(Ms)
    
    for i, z in enumerate(z_cat):
        for j, m in enumerate(M_cat):
            Data = data[(abs(data.zbest-z)<0.1) & (abs(data.M_med-m)<0.25)]
            NS_H[i,j] = (Data.FLUXERR_PETRO_F160W/ Data.FLUX_PETRO_F160W).median()
            
    zs = np.array([z_at_value(cosmo.age, t+0.5*u.Gyr) for t in Ages])
    
    def Interp_NS(z,M,NS):
        x = np.argmin(abs(z-z_cat))
        y = np.argmin(abs(M-M_cat))
        return NS[x,y]

    for i in range(Ages.size):
        for j in range(M_today.size):
            NS_mock[i,j] = Interp_NS(zs[i],np.log10(Ms[i,j]),NS_H)    

if make_table:
    table = Table()
    fluxes = np.squeeze(synphot.fluxes).value                   
    for i, n in enumerate(Names):

        if (real & multiple):  
            if multiple:
                N_sample = 10
            else:
                N_sample = 1                
            flux = (fluxes[i] * Ms[:,:]).ravel()
            error = NS_mock.ravel() * flux
            mflux = np.concatenate([flux for k in range(N_sample)])
            merror = np.concatenate([error for k in range(N_sample)])
            mnoise = np.array([np.random.normal(0, err) for err in merror])
            table.add_columns([Column(mflux+mnoise,'F%s'%(i+1)), Column(merror,'E%s'%(i+1))])

        else:
            flux = (fluxes[i] * Ms[:,:]).ravel()
            error = 0.01 * flux
            noise = [np.random.normal(0, err) for err in error]
            table.add_columns([Column(flux+noise,'F%s'%(i+1)), Column(error,'E%s'%(i+1))])

    id = Column(name='id', data=np.arange(1,len(table)+1)) 
    zspec = Column(name='zspec', data= -np.ones(len(table)))  
    table.add_column(id, 0) 
    table.add_column(zspec)  

    df = table.to_pandas()
    if perturb:
        if (real & multiple):    
            np.savetxt('Aldo/ST/mAldo_SFH_ST%.1f.cat'%A, table, header=' '.join(table.colnames),
                       fmt=['%d']+['%.5e' for i in range(20)]+['%.2f'])
        else:
            np.savetxt('Aldo/ST/Aldo_SFH_ST%.1f.cat'%A, table, header=' '.join(table.colnames),
                       fmt=['%d']+['%.5e' for i in range(20)]+['%.2f'])
    else:
        np.savetxt('Aldo/Aldo_SFH_elnm.cat', table, header=' '.join(table.colnames),
                   fmt=['%d']+['%.5e' for i in range(20)]+['%.2f'])