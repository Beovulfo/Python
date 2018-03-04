#!/usr/bin/python                                                                                                                  
## module run_kut4
from numpy import array,zeros
def integrate(F,x,y,xStop,h):
	def run_kut4(F,x,y,h):
# Computes increment of y from Eqs. (7.10)
		K0 = h*F(x,y)
		K1 = h*F(x + h/2.0, y + K0/2.0)
		K2 = h*F(x + h/2.0, y + K1/2.0)
		K3 = h*F(x + h, y + K2)
		return (K0 + 2.0*K1 + 2.0*K2 + K3)/6.0
	X = []
	Y = []
	X.append(x)
	Y.append(y)
	while x < xStop:
		h = min(h,xStop - x)
		y = y + run_kut4(F,x,y,h)
		x = x + h
		X.append(x)
		Y.append(y)
	return array(X),array(Y)
