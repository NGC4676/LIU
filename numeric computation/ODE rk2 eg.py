import numpy as np
import matplotlib.pyplot as plt
#==============================================================================
# A ODE solution sample of Spring F=kx-mg with Runge-Kunta 2th-order method
#==============================================================================
#define ODE going to solve
def ODE(state,time):
	z0=state[1]             #dx/dt=v
	z1=-k/m*state[0]-g      #dv/dt=-k/m*x-g
	return np.array([z0,z1]) #返回一个状态，即行向量
	
#define rk2 method
def rk2(state,t,dt,ODE):
	k0=dt*ODE(state,time)
	k1=dt*ODE(state+k0,time+dt)
	state_next=state+0.5*(k0+k1)
	return state_next

N=1000 #number of steps
[x0,v0]=[0.0,0.0] #initial state
T=10.0 #total simulation time
dt=T/float(N-1) #time step
k=3.5 #Spring constant
m=0.2 #mass
g=9.8 #gravity
time=np.linspace(0,T,N) #time
y=np.zeros([N,2]) #定义一个N行两列的状态空间
y[0]=[x0,v0] #指定第一行

for i in range(N-1):    #do calculation
	y[i+1]=rk2(y[i],time[i],dt,ODE)  #第i+1行由i行给出
x=[y[i,0] for i in range(N)]   #一维数组x，v
v=[y[i,1] for i in range(N)]

#plot
plt.figure()
plt.plot(time,x,label='position')
plt.plot(time,v,label='velocity')
plt.xlabel('time')
plt.ylabel('position / velocity')
plt.legend(loc='upper left')
plt.ylim(-6,6)
plt.show()

