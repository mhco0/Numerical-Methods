from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import numpy as np
import matplotlib.pyplot as plt
import sys
import math


def T(t):
	return (230.0*exp(-0.19018161948*t)+70.0)

def A(t):
	return (600.0 -550.0*exp(-t/100.0))
	
def P(t):
	return (16620.0*exp(0.11*t))
	
def i1(t):
	return (-60.0*t*exp(-100.0*t)-1.2*exp(-100.0*t)+1.2)
	
def i2(t):
	return (-120.0*t*exp(-100.0*t)-1.2*exp(-100.0*t)+1.2)

def main():

	op = input()
	arrt = []
	arry = []	
	
	if(op == '1'):
		for i in np.arange(0.0,20.0,0.1):
			arrt.append(i)
			arry.append(P(i))
		
		plt.ylabel('y(t) (Habitantes)')
	
	elif(op == '2'):
		for i in np.arange(0.0,20.0,0.1):
			arrt.append(i)
			arry.append(A(i))
			
		plt.ylabel('y(t) (Lb)')
	
	elif(op == '3'):
		for i in np.arange(0.0,20.0,0.1):
			arrt.append(i)
			arry.append(T(i))
		
		plt.ylabel('y(t) (ºF)')
		
	
	else:
		for i in np.arange(0.0,20.0,0.1):
			arrt.append(i)
			arry.append(i1(i))
		
		plt.ylabel('y(t) (A)')
		plt.xlabel('t (m)')
		plt.title('Grafico da Solucao Analitica da Corrente 1')
		plt.plot(arrt,arry,color='silver')
		plt.show()
		
		arrt.clear()
		arry.clear()
		
		for i in np.arange(0.0,20.0,0.1):
			arrt.append(i)
			arry.append(i2(i))
		
		plt.ylabel('y(t) (A)')
		
		
	plt.xlabel('t (m)')
	plt.plot(arrt,arry,color='black')
	if(op != '4'):
		plt.title('Grafico da Solucao Analitica')
	else:
		plt.title('Grafico da Solucao Analitica da Corrente 2')
	plt.show()
	return

if __name__ == '__main__':
	main()
