from math import *
#==============================================================================
# 2-point method
#==============================================================================
def f(x):
	return x**2+1 
def derivative_2point(f,x,d):       
	y=(f(x+d)-f(x))/d
	return "Derivative is y'= %.6f" %y

print derivative_2point(f,2,0.1)
#==============================================================================
# 3-point method
#==============================================================================
def g(x):
	return x**2+1 
def derivative_3point(g,x,d):       
	y=(g(x+d)-g(x-d))/(2*d)
	return "Derivative is y'= %.6f" %y

print derivative_3point(g,2,0.1)
#==============================================================================
# 5-point method
#==============================================================================
def h(x):
	return x**2+1 
def derivative_5point(h,x,d):       
	y=(h(x-2*d)-8*h(x-d)+8*h(x+d)-h(x+2*d))/(12*d)
	return "Derivative is y'= %.6f" %y

print derivative_5point(h,2,0.1)