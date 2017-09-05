import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sci

N=1000                    #you can change system parameters here
y=np.zeros([4])
l0=1.0
l=1.0
v0=0.0
theta0=0.3
omega0=0.0

y[0]=l                    #you can change initial state here
y[1]=v0
y[2]=theta0
y[3]=omega0

time=np.linspace(0,50,N)  #you can change simulation time here
k=3
m=0.2
g=9.8

def spring_pendulum(y,time):
	z0=y[1]
	z1=(l0+y[0])*y[3]**2-k/m*y[0]+g*np.cos(y[2])
	z2=y[3]
	z3=-(g*np.sin(y[2])+2*y[1]*y[3])/(l0+y[0])
	return np.array([z0,z1,z2,z3])

answer=sci.odeint(spring_pendulum,y,time) # ODE,state,time

l=[answer[i,0] for i in range(N)]
v=[answer[i,1] for i in range(N)]

x_pos=(l0+answer[:,0])*np.sin(answer[:,2])
y_pos=-(l0+answer[:,0])*np.cos(answer[:,2])

plt.figure()
plt.plot(time,l,label='position')
plt.plot(time,v,label='velocity')
plt.xlabel('time')
plt.ylabel('position / velocity')
plt.legend(loc='upper left')

plt.figure()
plt.plot(x_pos,y_pos)
plt.xlabel('x')
plt.ylabel('y')

plt.show()

