from math import *
import numpy as np
#==============================================================================
# Simple Method
#==============================================================================
def f(x):
	return x**2+1  #input f(x)
def integral_simple(f,a,b,N):        #N为长方体个数（slice）
	d=(b-a)/float(N)
	S=0
	for i in range(N)	:
		k=a+i*d     #k=xn
		S+=f(k)*d
	return "Integral is S= %.6f" %S

print integral_simple(f,0,1,100)
#==============================================================================
# Trapezoid Method
#==============================================================================
def g(x):
	return x**2+1
def integral_trapezoid(g,a,b,N):
	d=(b-a)/float(N)
	S=(g(b)+g(a))*d/2.0
	for i in range(N-1):
		k=a+(i+1)*d      #k=xn
		S+=g(k)*d        #integral=[(f(x0)+f(xn))/2+f(x1)+...+f(xn-1)]*d
	return "Integral is S= %.6f" %S

print integral_trapezoid(g,0,1,100)
#==============================================================================
# Simpson Method
#==============================================================================
def h(x):
	return x**2+1
def integral_simpson(h,a,b,N):
	d=(b-a)/float(N)
	S=0
	for i in np.arange(0,N-1,2):
		k=a+(i+1)*d
		S+=(h(k-d)+4*h(k)+h(k+d))*d/3.0   #k=xn, n取1 - N-1奇数
	if N % 2 == 1:     #若是奇数个slice还须加上最后一块
		S+=(5*h(k+2*d)+8*h(k+d)-h(k))*d/12.0	 #A=[(5f(xn)+8f(xn-1)-f(n-2)]*d/12
	return S

print "Integral is S= %.6f" % integral_simpson(h,0,1,100)