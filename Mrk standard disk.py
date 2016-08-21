from matplotlib import pyplot as plt
import xlrd  
import numpy as np
import scipy.constants as const
from scipy import integrate

#read and plot spectrum
data = xlrd.open_workbook(r"C:\Users\mac\Desktop\NAOC\data\Mrk-231.xlsx")   
table=data.sheets()[0]
wave=np.array(table.col_values(0)[:8163])
spec=np.array(table.col_values(1)[:8163])
spec[7000:7200]=None
spec[6460:6800]=None
spec[5560:5645]=None
spec[:2000]=None

fig = plt.subplots(figsize=(10, 3))	
plt.xlabel('Wavelength',fontsize=15)
plt.ylabel('log Flux',fontsize=15)
plt.plot(wave,spec)


#define constant and parameter

#constants
G=const.gravitational_constant
c=const.c
h=const.h
sigma_SB=const.Stefan_Boltzmann
pi=const.pi
k_B=const.Boltzmann

#parameter
E_BV=0.1      #extinction
M_t=2.8e38     #total mass kg
M_sun=2e30   #solar mass kg
a_BBH=0.4e14   #semimajor axis at t
cosi=0.8   #circumbinary inclination
q=0.02    #mass ratio
M_s=q*M_t/(1+q)  #total mass
cosj=0.8   #mini inclination
f_s_Edd=0.6  #mini Eddington factor
f_c_Edd=0.5  #circumbinary Eddington factor
L_s_Edd=1.26e31*(M_s/M_sun)   #mini Eddington luminosity  kg/s
L_c_Edd=1.26e31*(M_t/M_sun)   #circum Eddington luminosity
M_s_acc=f_s_Edd*L_s_Edd/(0.1*c**2)  #mini accretion rate
M_c_acc=f_c_Edd*L_c_Edd/(0.1*c**2)  #mini accretion rate
print M_s_acc
print M_c_acc

r_in=a_BBH/(1+q)+0.69*a_BBH*q**(1.0/3)
r_out=1e5*G*M_t/c**2

r_s_in=3.5*G*M_s/c**2  #mini inner radius
R_RL=0.49*a_BBH*q**(2.0/3)/(0.6*q**(2.0/3)+np.log(1+q**0.5))  #mean roche radius
f_rs=0.11   #mini disk outer ratio to roche
r_s_out=f_rs*R_RL    #mini inner radius

print r_s_in,r_s_out
print r_in,r_out

def Extinction(lamda):    #extinction=F0/F'
	R=4.0
	lamda_mu=lamda*1e6
	def k(lamda_mu):
		if lamda_mu>0.63:
			return 2.659*(-1.857+1.040/lamda_mu)+R
		else:
			return 2.659*(-2.156+1.509/lamda_mu-0.198/lamda_mu**2+0.011/lamda_mu**3)+R
	return 10**(0.4*0.44*E_BV*k(lamda_mu))

#mini-disk model
def t_eff(r):     #effective temperature
#	return 2e5*((f_s_Edd/0.3)*(2e38/M_s)*(1-np.sqrt(r_s_in/r)))**0.25*(r_s_in/r)**0.75
	return ((3*G*M_s*M_s_acc*(1-np.sqrt(r_s_in/r)))/(8*pi*sigma_SB*r**3))**0.25

def f_lamda(lamda):   #lamda flux

	def g(r):
		return 2*pi*r*(2*h*c**2*cosj/lamda**5)/(np.exp(h*c/(lamda*k_B*t_eff(r)))-1)/Extinction(lamda)

	f_lamda=integrate.quad(g, r_s_in, r_s_out)

	return f_lamda

def T_eff(r):     #effective temperature
#	return 2e5*((f_c_Edd/0.3)*(2e38/M_t)*(1-np.sqrt(r_in/r)))**0.25*(r_in/r)**0.75
	return ((3*G*M_t*M_c_acc*(1-np.sqrt(r_in/r)))/(8*pi*sigma_SB*r**3))**0.25

def F_lamda(lamda):   #lamda flux
	def f(r):
		return 2*pi*r*(2*h*c**2*cosi/lamda**5)/(np.exp(h*c/(lamda*k_B*T_eff(r)))-1)/Extinction(lamda)	

	F_lamda=integrate.quad(f, r_in, r_out)

	return F_lamda	

#plot spectrum
x=np.linspace(1500,8000,200)
flux1=[]
for i in x:
	flux1.append(0.5e-56*F_lamda(1e-10*i)[0])

flux2=[]
for i in x:
	flux2.append(0.5e-56*f_lamda(1e-10*i)[0])	

flux3=[]
for i in range(200):
	flux3.append(flux1[i]+flux2[i])	

plt.plot(x,flux1,x,flux2,x,flux3,'.')
plt.yscale('log')
plt.xlim(1000,8000)
plt.ylim(1e-16,1e-13)
