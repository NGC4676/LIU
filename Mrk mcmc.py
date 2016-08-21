from matplotlib import pyplot as plt
import xlrd  
import numpy as np
import scipy.constants as const
from scipy import integrate
import pymc

#read and plot spectrum

data = xlrd.open_workbook(r"C:\Users\mac\Desktop\NAOC\data\Mrk-231 mask.xlsx")   
table = data.sheets()[0]
wave = np.array(table.col_values(0)[3000:8163:50])
lamda = 1e-10*wave
spec = np.array(table.col_values(1)[3000:8163:50])

#define constant and parameter

#constants
G = const.gravitational_constant
c = const.c
h = const.h
sigma_SB = const.Stefan_Boltzmann
pi = const.pi
k_B = const.Boltzmann
M_sun = 2e30   #solar mass kg

#parameter

cosi = 0.8   #circumbinary inclination
cosj = 0.8   #mini inclination


#M_t = 2.8e38
M_t=pymc.Uniform('M_t', 1e36, 1e39)    #total mass kg
#q = 0.02                                    
q = pymc.Uniform('q', 0.001, 1)               #mass ratio
#a_BBH=pymc.Normal('a_BBH', 0.4e14, 1e-24)    
a_BBH = pymc.Uniform('a_BBH', 0.1e13, 1e15)    #semimajor axis at t
#A=pymc.Normal('A',0.5e-56,1e114)  
A = pymc.Uniform('A', 0.1e-56, 1e-56)  #scale factor
#f_s_Edd=0.6    #mini Eddington factor
f_s_Edd = pymc.Uniform('f_s_Edd', 0.1, 1)   
#f_c_Edd=0.5    #circumbinary Eddington factor
f_c_Edd = pymc.Uniform('f_c_Edd', 0.1, 1)   
#f_RS = 0.11    #ratio of outer mini disk
f_RS = pymc.Uniform('f_RS', 0.01, 0.9)
E_BV=0.1      #extinction
#E_BV = pymc.Uniform('E_BV', 0.01, 0.9)
@pymc.deterministic
def M_s(q = q, M_t = M_t):
	return q*M_t/(1+q)           #small BH mass

@pymc.deterministic
def M_s_acc(f_s_Edd = f_s_Edd, M_s = M_s):  
	L_s_Edd = 1.26e31*(M_s/M_sun)   #mini Eddington luminosity  kg/s
	return f_s_Edd*L_s_Edd/(0.1*c**2)  #mini accretion rate
@pymc.deterministic
def M_c_acc(f_c_Edd = f_c_Edd, M_t = M_t): 
	L_c_Edd = 1.26e31*(M_t/M_sun)   #circum Eddington luminosity
	return f_c_Edd*L_c_Edd/(0.1*c**2)  #mini accretion rate

@pymc.deterministic
def r_in(a_BBH = a_BBH, q=q):
	return a_BBH/(1+q) + 0.69*a_BBH*q**(1.0/3)

@pymc.deterministic
def r_out(M_t = M_t):
	return 1e5*G*M_t/c**2

@pymc.deterministic
def r_s_in(M_s = M_s):
	return 3.5*G*M_s/c**2

@pymc.deterministic
def r_s_out(a_BBH = a_BBH, q = q, f_RS = f_RS):
	R_RL = 0.49*a_BBH*q**(2.0/3)/(0.6*q**(2.0/3) + np.log(1+q**0.5))  #mean roche radius
	return f_RS*R_RL    #mini inner radius

#Extinction

def Extinction(lamda):    #extinction=F0/F'
	R = 4.0      #
	lamda_mu = lamda*1e6
	out=[]
	for i in range(len(lamda_mu)): 	
		def k(i):
			if i > 0.63:
				return 2.659*(-1.857 + 1.040/i) + R
			else:
				return 2.659*(-2.156 + 1.509/i - 0.198/i**2 + 0.011/i**3) + R
		out.append(10**(0.4*0.44*E_BV*k(lamda_mu[i])))
	return out

#Circumbinary-disk model
	
@pymc.deterministic
def F_lamda(r_in = r_in, r_out = r_out, M_t = M_t, M_c_acc = M_c_acc):   #flux
	global lamda
	r = np.linspace(r_in, r_out, 1000).reshape(-1,1)

	def T_eff(r):     #effective temperature
		return ( (3*G*M_t*M_c_acc*(1-np.sqrt(r_in/r)))/(8*pi*sigma_SB*r**3) )**0.25

	def f(r):
		return 2*pi*r*(2*h*c**2*cosi/lamda**5)/(np.exp(h*c/(lamda*k_B*T_eff(r)))-1)/Extinction(lamda)
	
	F_lamda = integrate.simps(f(r), r, axis=0)   #integrate with Simpson method

	return F_lamda

#mini-disk model

@pymc.deterministic
def f_lamda(r_s_out = r_s_out, r_s_in = r_s_in, M_s = M_s, M_s_acc = M_s_acc):   #flux
	global lamda
	r = np.linspace(r_s_in, r_s_out, 1000).reshape(-1,1)

	def t_eff(r):     #effective temperature
		return ((3*G*M_s*M_s_acc*(1-np.sqrt(r_s_in/r)))/(8*pi*sigma_SB*r**3))**0.25
		
	def g(r):
		return 2*pi*r*(2*h*c**2*cosj/lamda**5)/(np.exp(h*c/(lamda*k_B*t_eff(r)))-1)/Extinction(lamda)
	
	f_lamda = integrate.simps(g(r), r, axis=0)
	
	return f_lamda

@pymc.deterministic
def S(f_lamda = f_lamda, F_lamda = F_lamda, A = A):
	return A*(f_lamda+F_lamda)/Extinction(lamda)

T = pymc.Normal('T', S , 1e33, value=spec, observed=True)

model = dict(a_BBH = a_BBH, q = q, M_t = M_t, M_s=M_s,
													r_in = r_in, r_out = r_out,
													r_s_out = r_s_out, r_s_in = r_s_in,
													f_s_Edd = f_s_Edd, f_c_Edd = f_c_Edd,
													M_c_acc = M_c_acc, M_s_acc = M_s_acc,
													f_lamda = f_lamda, F_lamda = F_lamda,
													f_RS=f_RS, A = A, S = S, T = T)

M = pymc.MCMC(model)

M.sample(iter = 10000, burn = 5000, thin = 10) 

print
#Trace result
trace1 = M.trace('a_BBH')
trace2 = M.trace('q')
trace3 = M.trace('f_s_Edd')
trace4 = M.trace('f_c_Edd')
trace5 = M.trace('f_RS')
trace6 = M.trace('M_t')
trace7 = M.trace('A')

#Assign params
a_bbh = np.median(trace1[:])
p = np.median(trace2[:])
f_s_edd = np.median(trace3[:])
f_c_edd = np.median(trace4[:])
f_rs = np.median(trace5[:])
m_t = np.median(trace6[:])
a = np.median(trace7[:])

print 'a_bbh = ',a_bbh , 'q = ',p 
print 'f_s_edd = ',f_s_edd , 'f_c_edd = ',f_c_edd
print 'f_rs = ',f_rs , 'm_t = ',m_t
print 'a = ',a

plt.figure(figsize = (10, 12))
ax1 = plt.subplot(421)
ax2 = plt.subplot(422)
ax3 = plt.subplot(423)
ax4 = plt.subplot(424)
ax5 = plt.subplot(425)
ax6 = plt.subplot(426)
ax7 = plt.subplot(427)

ax1.hist(trace1[:], 25, histtype='step')
ax1.set_title('a_BBH')
ax2.hist(trace2[:], 25, histtype='step')
ax2.set_title('q')
ax3.hist(trace3[:], 25, histtype='step')
ax3.set_title('f_s_Edd')
ax4.hist(trace4[:], 25, histtype='step')
ax4.set_title('f_c_Edd')
ax5.hist(trace5[:], 25, histtype='step')
ax5.set_title('f_rs')
ax6.hist(trace6[:], 25, histtype='step')
ax6.set_title('M_t')
ax7.hist(trace7[:], 25, histtype='step')
ax7.set_title('A')

plt.show()
#==============================================================================
# replot fitting result
#==============================================================================
fig = plt.subplots(figsize = (10, 3))
#secondary parameter
m_s=p*m_t/(1+p)
r_in = a_bbh/(1+p) + 0.69*a_bbh*p**(1.0/3)
r_out = 1e5*G*m_t/c**2
R_RL = 0.49*a_bbh*p**(2.0/3) / (0.6*p**(2.0/3) + np.log(1+p**0.5))
r_s_out = f_rs*R_RL    
r_s_in=3.5*G*m_s/c**2  
L_s_edd = 1.26e31*(m_s/M_sun)   
L_c_edd = 1.26e31*(m_t/M_sun)   #circum Eddington luminosity
m_s_acc = f_s_edd*L_s_edd / (0.1*c**2)  #mini accretion rate
m_c_acc = f_c_edd*L_c_edd / (0.1*c**2)  #mini accretion rate

#Extinction
def extinction(lamda):    #extinction=F0/F'
	R = 4.0
	lamda_mu = lamda*1e6
	def k(lamda_mu):
		if lamda_mu > 0.63:
			return 2.659*(-1.857 + 1.040/lamda_mu) + R
		else:
			return 2.659*(-2.156 + 1.509/lamda_mu-0.198/lamda_mu**2 + 0.011/lamda_mu**3) + R
	return 10**(0.4*0.44*E_BV*k(lamda_mu))
	
#mini-disk model
def t_eff(r):     #effective temperature
	return ( (3*G*m_s*m_s_acc*(1-np.sqrt(r_s_in/r)))/(8*pi*sigma_SB*r**3) )**0.25

def flamda(lamda):   #lamda flux

	def g(r):
		return 2*pi*r*(2*h*c**2*cosj/lamda**5)/(np.exp(h*c/(lamda*k_B*t_eff(r)))-1)/extinction(lamda)

	flamda = integrate.quad(g, r_s_in, r_s_out)

	return flamda
#Circumbinary disk model
def T_eff(r):     #effective temperature
	return ( (3*G*m_t*m_c_acc*(1-np.sqrt(r_in/r)))/(8*pi*sigma_SB*r**3) )**0.25

def Flamda(lamda):   #lamda flux
	
	def f(r):
		return 2*pi*r*(2*h*c**2*cosi/lamda**5)/(np.exp(h*c/(lamda*k_B*T_eff(r)))-1)/extinction(lamda)	

	Flamda = integrate.quad(f, r_in, r_out)

	return Flamda	

#plot spectrum
x=np.linspace(1500, 8000, 200)
flux1 = []
for i in x:
	flux1.append(a*Flamda(1e-10*i)[0])

flux2 = []
for i in x:
	flux2.append(a*flamda(1e-10*i)[0])	

flux3 = []
for i in range(200):
	flux3.append(flux1[i] + flux2[i])	

plt.plot(x,flux1, x,flux2, x,flux3,'--')
plt.yscale('log')
plt.xlim(1000,8000)
plt.ylim(1e-16,1e-13)
plt.plot(wave,spec)