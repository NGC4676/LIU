import  matplotlib.pyplot as plt
import numpy as np
from random import randint

plt.ion()
plt.figure(figsize=(10,10))

atoms=np.ones([400,2])*100
line, =plt.plot(atoms[:,0],atoms[:,1],'ro')
plt.xlim(0,200)
plt.ylim(0,200)
plt.draw()
N=100
for i in range(N):
	for j in range(N):
		atoms[j,0]+=randint(-5,5)
		atoms[j,1]+=randint(-5,5)	
		x,y=(atoms[j,0],atoms[j,1])
		if x == 200 :
			atoms[j,0]=198
		elif x == 0 :
			atoms[j,0]=2
		if y == 200 :
			atoms[j,1]=198
		elif y== 0 :
			atoms[j,1]=2
line.set_xdata(atoms[:,0])
line.set_ydata(atoms[:,1])
plt.draw()	
#1wait=raw_input('continue\n')
		
