from matplotlib import pyplot as plt
import xlrd  
import numpy as np
import scipy.constants as const
from scipy import integrate

#read and plot spectrum
data = xlrd.open_workbook(r"C:\Users\mac\Desktop\Mrk-231.xlsx")   
table=data.sheets()[0]
wavelen=table.col_values(0)[:8163]
spec=np.array(table.col_values(1)[:8163])

fig = plt.subplots(figsize=(10, 3))	
plt.xlabel('Wavelength',fontsize=15)
plt.ylabel('log Flux',fontsize=15)
plt.plot(wavelen,spec)
plt.yscale('log')


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
Nu=1      #viscosity
#Sigma=0.008   #surface density
Sigma=0.5e22   #surface density
M_t=3e38      #total mass
a_BBH=0.885e14   #semimajor axis at t
cosi=0.5   #circumbinary inclination
q=0.025     #mass ratio
M_s=q*M_t/(1+q)  #total mass
cosj=0.8   #mini inclination
f_s_Edd=0.5  #mini Eddington factor
L_s_Edd=1.26e31*(M_s/(2e30))   #Eddington luminosity
M_acc=f_s_Edd*L_s_Edd/(0.1*c**2)  #mini accretion rate

r_in=a_BBH/(1+q)+0.69*a_BBH*q**(1.0/3)
r_out=1e5*G*M_t/c**2
r_s_in=3.5*G*M_s/c**2  #mini inner radius
R_RL=0.49*a_BBH*q**(2.0/3)/(0.6*q**(2.0/3)+np.log(1+q**0.5))  #mean roche radius
f_rs=0.11    #mini disk outer ratio to roche
r_s_out=f_rs*R_RL    #mini inner radius

print r_s_in,r_s_out,r_in,r_out

def Extinction(lamda):    #extinction=Fi/F0
	R=4.0
	lamda_mu=lamda*1e6
	def k(lamda_mu):
		if lamda_mu>0.63:
			return 2.659*(-1.857+1.040/lamda_mu)+R
		else:
			return 2.659*(-2.156+1.509/lamda_mu-0.198/lamda_mu**2+0.011/lamda_mu**3)+R
	return 10**(0.4*0.44*E_BV*k(lamda_mu))
print  Extinction(3e-7)

#mini-disk model
def t_eff(r):     #effective temperature
	return ((3*G*M_s*M_acc*(1-np.sqrt(r_s_in/r)))/(8*pi*sigma_SB*r**3))**0.25

def f_lamda(lamda):   #lamda flux

	def g(x):
		return 2*pi*x*(2*h*c**2*cosj/lamda**5)/(np.exp(h*c/(lamda*k_B*t_eff(x)))-1)/Extinction(lamda)

	f_lamda=integrate.quad(g, r_s_in, r_s_out)

	return f_lamda
	
#Circumbinary disk model
def Omega(r):    #angular velocity
	return np.sqrt(G*M_t/r**3) 

def F_J(r):      #angular momentum flux 
	return 3*pi*Nu*Sigma*Omega(r)*r**2   

def T_eff(r):     #effective temperature
	return ((3.0/(8*pi))*((F_J(r)*Omega(r))/(sigma_SB*r**2)))**0.25

def F_lamda(lamda):   #lamda flux

#	def r_in(a_BBH,q):     #inner disk 
#		def R_H(a_BBH,q):   #Hill radius
#			return 0.69*a_BBH*q**(1.0/3)
#		return a_BBH/(1+q)+R_H
#	def r_out(M_t):  #outer disk
#		return 1e5*G*M_t/c**2

	def f(x):
		return 2*pi*x*(2*h*c**2*cosi/lamda**5)/(np.exp(h*c/(lamda*k_B*T_eff(x)))-1)/Extinction(lamda)	

	F_lamda=integrate.quad(f, r_in, r_out)

	return F_lamda	
	
#plot spectrum
x=np.linspace(1000,8000,200)

flux1=[]
for i in x:
	flux1.append(5e-57*F_lamda(1e-10*i)[0])

flux2=[]
for i in x:
	flux2.append(5e-57*f_lamda(1e-10*i)[0])	

flux3=[]
for i in range(200):
	flux3.append(flux1[i]+flux2[i])	
plt.plot(x,flux1,x,flux2,x,flux3,'.')

plt.yscale('log')
plt.xlim(1000,8000)
plt.ylim(1e-16,1e-13)

fig = plt.subplots(figsize=(7, 5))
r=np.linspace(r_s_in,r_s_out,200)
plt.plot(r,t_eff(r))

fig = plt.subplots(figsize=(7, 5))
r=np.linspace(r_in,r_out,200)
plt.plot(r,T_eff(r))
