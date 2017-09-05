from random import random
import math 
#==============================================================================
#use Monte Carlo method  to find Volume of a hemisphere r=1
#==============================================================================
N=1000000  #number of random points 
count=0    #number of points in the sphere
for j in range(N):
	point=(random(),random(),random())
	if math.sqrt(point[0]**2+point[1]**2+point[2]**2)<1:
		count+=1
answer=4*float(count)/float(N) 
print 'the volume of a hemispehre is %0.4f+/- %0.4f' % (answer,4.0/math.sqrt(N))