import numpy as np
import matplotlib.pyplot as plt
#==============================================================================
# A ODE solution sample of Spring F=kx-mg with Euler method(easy to diverge)
#==============================================================================
#define ODE going to solve
def ODE(state,t):  #state为方程左边列向量，z为右边列向量
	z0=state[1]             #dx/dt=v
	z1=-k/m*state[0]-g      #dv/dt=-k/m*x-g
	return np.array([z0,z1]) #返回一个状态，即行向量
	
#define euler method
def euler(state,t,dt,ODE):
	state_next=state+ODE(state,t)*dt  #给出下一时刻状态向量
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
	y[i+1]=euler(y[i],time[i],dt,ODE)  #第i+1行由i行给出
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
