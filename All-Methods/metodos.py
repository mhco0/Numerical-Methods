from sympy import *
from sympy.parsing.sympy_parser import parse_expr
#################### Auxi ####################
def EulerPd(y0,t0,h,funct):
	t,y = symbols('t y')
	yn = y0
	return yn + h*funct.subs({t:t0,y:yn})


#################### Methods ####################
def Euler(y0,t0,h,qp,funct):
	print('Metodo de Euler')
	print('y({}) = {}'.format(t0,y0))
	print('h = {}'.format(h))

	t,y = symbols('t y')
	yn = y0
	for i in range(qp+1):
		print(i,' ',yn)
		yn = yn + h*funct.subs({t:t0,y:yn})
		t0 = t0 + h
	return 

def Euler_Inverso(y0,t0,h,qp,funct):
	print('Metodo de Euler Inverso')
	print('y({}) = {}'.format(t0,y0))
	print('h = {}'.format(h))

	t,y = symbols('t y')
	yn = y0
	for i in range(qp+1):
		print(i,' ',yn)
		yaux = EulerPd(yn,t0,h,funct)
		taux = t0 + h
		yn = yn + h*funct.subs({t:taux,y:yaux})
		t0 = t0 + h 
	return

def Euler_Aprimorado(y0,t0,h,qp,funct):
	print('Metodo de Euler Aprimorado')
	print('y({}) = {}'.format(t0,y0))
	print('h = {}'.format(h))

	t,y = symbols('t y')
	yn = y0
	for i in range(qp+1):
		print(i,' ',yn)

		k1 = funct.subs({t:t0,y:yn})
		k2 = funct.subs({t:(t0+h),y:yn+h*k1})
		yn = yn + float((h*(k1+k2))/2)
		t0 = t0 + h
	return

def Runge_Kutta(y0,t0,h,qp,funct):
	print('Metodo de Runge-Kutta')
	print('y({}) = {}'.format(t0,y0))
	print('h = {}'.format(h))

	t,y = symbols('t y')
	yn = y0
	for i in range(qp+1):
		print(i,' ',yn)

		k1 = funct.subs({t:t0,y:yn})
		k2 = funct.subs({t:(t0+h/2),y:yn+(h/2)*k1})
		k3 = funct.subs({t:(t0+h/2),y:yn+(h/2)*k2})
		k4 = funct.subs({t:(t0+h),y:yn+h*k3})
		yn = yn + float((h*(k1+2*k2+2*k3+k4))/6)
		t0 = t0 + h
	return

def Adam_Bashforth(y0,t0,h,qp,funct,order):
	print('Metodo de Adam-Bashforth')
	print('y({}) = {}'.format(t0,y0))
	print('h = {}'.format(h))

	coefM = [
			[1],
			[3.0/2,-1.0/2],
			[23.0/12,-4.0/3,5.0/12],
			[55.0/24,-59.0/24,37.0/24,-3.0/8],
			[1901.0/720,-1387.0/360,109.0/30,-637.0/360,251.0/720],
			[4277.0/1440,-2641.0/480,4991.0/720,-3649.0/720,959.0/480,-95.0/288],
			[198721.0/60480,-18637.0/2520,235183.0/20160,-10754.0/945,135713.0/20160,-5603.0/2520,19087.0/60480],
			[16083.0/4480,-1152169.0/120960,242653.0/13440,-296053.0/13440,2102243.0/120960,-115747.0/13440,32863.0/13440,-5257.0/17280]
			]


	return 
def Adam_Moulton():
	return

def Formula_Inversa():
	return

def main():
	entry = str(input()).split()

	method = entry[0]

	if(method == 'euler'):
		y0 = float(entry[1])
		t0 = float(entry[2])
		h = float(entry[3])
		qp = int(entry[4])
		funct = parse_expr(entry[5])

		Euler(y0,t0,h,qp,funct)
	elif(method == 'euler_inverso'):
		y0 = float(entry[1])
		t0 = float(entry[2])
		h = float(entry[3])
		qp = int(entry[4])
		funct = parse_expr(entry[5])

		Euler_Inverso(y0,t0,h,qp,funct)
	elif(method == 'euler_aprimorado'):
		y0 = float(entry[1])
		t0 = float(entry[2])
		h = float(entry[3])
		qp = int(entry[4])
		funct = parse_expr(entry[5])

		Euler_Aprimorado(y0,t0,h,qp,funct)
	elif(method == 'runge_kutta'):
		y0 = float(entry[1])
		t0 = float(entry[2])
		h = float(entry[3])
		qp = int(entry[4])
		funct = parse_expr(entry[5])

		Runge_Kutta(y0,t0,h,qp,funct)
	elif(method == 'adam_bashforth'):
		yl = []
		size = len(entry)

		for i in range(int(entry[size-1])):
			yl.append(entry[i+1])

		y0 = float(entry[size-6])
		t0 = float(entry[size-5])
		h = float(entry[size-4])
		qp = int(entry[size-3])
		funct = parse_expr(entry[size-2])
		order = int(entry[size-1])
		
		Adam_Bashforth(y0,t0,h,qp,funct,order)

	# elif(method == 'adam_bashforth_by_euler'):

	# elif(method == 'adam_bashforth_by_euler_inverso'):

	# elif(method == 'adam_bashforth_by_euler_aprimorado'):

	# elif(method == 'adam_bashforth_by_runge_kutta'):

	# elif(method == 'adam_multon'):

	# elif(method == 'adam_multon_by_euler'):

	# elif(method == 'adam_multon_by_euler_inverso'):

	# elif(method == 'adam_multon_by_euler_aprimorado'):

	# elif(method == 'adam_multon_by_runge_kutta'):

	# elif(method == 'formula_inversa'):

	# elif(method == 'formula_inversa_by_euler'):

	# elif(method == 'formula_inversa_by_euler_inverso'):

	# elif(method == 'formula_inversa_by_euler_aprimorado'):

	# elif(method == 'formula_inversa_by_runge_kutta'):

	return 0

if __name__ == "__main__":
	main()