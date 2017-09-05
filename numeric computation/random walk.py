from random import choice
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
N=200             #steps
n=200             #particles
x=np.zeros([N])
t=range(N)
x_sum=np.zeros([N])
x2_sum=np.zeros([N])

for j in range(n):         #循环n遍模拟n个粒子
	for i in range(1,N):    #由于序列总是从x[0]=0出发，所以各个循环（粒子）间不相关
		if choice(['forward','back']) == 'back':
			x[i]=x[i-1]-1
		else:
			x[i]=x[i-1]+1
	for k in range(1,N):	
		x_sum[k]+=x[k]
		x2_sum[k]+=x[k]**2
		
x_average=[float(xi)/n for xi in x_sum]
RMS=[np.sqrt(float(x2i)/n) for x2i in x2_sum]	

plt.plot(t,x_average)
plt.plot(t,RMS,'g')
plt.xlabel('time')
plt.ylabel('Average / RMS position')

def Power(x,a,b):
	return a*x**b

opt,cov=curve_fit(Power,t,RMS)

print 'Power fit: y=At^B'
print 'A=%0.6f +/- %0.6f' %(opt[0],np.sqrt(cov[0,0]))
print 'B=%0.6f +/- %0.6f' %(opt[1],np.sqrt(cov[1,1]))

plt.plot(t,Power(t,opt[0],opt[1]),'r')
plt.show()