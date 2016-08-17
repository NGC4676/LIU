# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 22:15:01 2016

@author: mac
"""
from math  import *
import numpy as np
#==============================================================================
# Bisection Method f(x)=0
#==============================================================================
def f(x):     #如何把某个式子f作为参数传递给bisection？
	return exp(x)*log(x)-x**2    #input f(x)
def root_bisection(f,start_guess,end_guess):
	tolerance=1.0e-6
	a,b=start_guess,end_guess
	l=abs(b-a)
	while l>tolerance:
		x=(a+b)/2.0
		if f(a)*f(x)<0:
			b=x
		else:
			a=x
		l=abs(b-a)
	return (x,tolerance)	

print 'f(x)=0 at x = %.6f +/-%.6f' %root_bisection(f,0.1,3)     
#==============================================================================
# Newtonian Method g(x)=0
#==============================================================================
def g(x):     
	return exp(x)*log(x)-x**2
def dg(x):
	return exp(x)*(1.0/x)+exp(x)*log(x)-2*x
def root_newton(g,dg,start_guess):
	tolerance=1.0e-6
	x=start_guess
	l=2*tolerance
	while l>tolerance:
		x1=x-g(x)/dg(x)
		l=abs(x-x1)
		x=x1
	return (x,tolerance)	

print 'g(x)=0 at x = %.6f +/-%.6f' % root_newton(g,dg,0.1)      
#==============================================================================
# Secant Method h(x)=0
#==============================================================================
def h(x):     
	return exp(x)*log(x)-x**2
def root_secant(h,start_guess,end_guess):
	tolerance=1.0e-6
	x=start_guess
	x2=end_guess
	l=2*tolerance
	while l>tolerance:
		x1=x2
		k=(h(x1)-h(x))/(x1-x)
		x2=x-h(x)/k
		l=abs(x2-x1)
		x=x1
	return (x,tolerance)	

print 'h(x)=0 at x = %.6f +/-%.6f' % root_secant(h,0.1,3)      