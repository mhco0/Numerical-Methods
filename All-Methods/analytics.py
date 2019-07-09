from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import numpy as np
import matplotlib.pyplot as plt
import sys
import math


def T(t):
	return (230.0*exp(-0.19018161948*t)+70.0)
	

def main():
	arrt = []
	arry = []	
	
	for i in np.arange(0.0,20.0,0.1):
		arrt.append(i)
		arry.append(T(i))
		
	plt.xlabel('t (m)')
	plt.ylabel('y(t) (ºF)')
	plt.title('Grafico da Solucao Analitica')
	plt.plot(arrt,arry,color='black')
	plt.show()
	
	return

if __name__ == '__main__':
	main()
