import numpy as np
import matplotlib.pyplot as plt
#==============================================================================
#toolbox 	
#==============================================================================

#plot func
def plot(xdata,ydata,xname,yname):
	plt.figure()
	plt.plot(xdata,ydata,label=yname)
	plt.xlabel(xname)
	plt.ylabel(yname)
	
#ODE euler method
def euler(state,t,dt,ODE):
	state_next=state+ODE(state,t)*dt  #给出下一时刻状态向量
	return state_next

#ODE rk2 method
def rk2(state,t,dt,ODE):
	k0=dt*ODE(state,time)
	k1=dt*ODE(state+k0,time+dt)
	state_next=state+0.5*(k0+k1)
	return state_next
	