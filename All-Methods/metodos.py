from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import sys
#################### Auxi ####################
def EulerPd(y0,t0,h,funct):
	t,y = symbols('t y')
	return float(y0 + h*funct.subs({t:t0,y:y0}))

def Adam_BashforthPd(yl,t0,h,funct,order):

	coefM = [
			[1.0],
			[3.0/2.0,-1.0/2.0],
			[23.0/12.0,-4.0/3.0,5.0/12.0],
			[55.0/24.0,-59.0/24.0,37.0/24.0,-3.0/8.0],
			[1901.0/720.0,-1387.0/360.0,109.0/30.0,-637.0/360.0,251.0/720.0],
			[4277.0/1440.0,-2641.0/480.0,4991.0/720.0,-3649.0/720.0,959.0/480.0,-95.0/288.0],
			[198721.0/60480.0,-18637.0/2520.0,235183.0/20160.0,-10754.0/945.0,135713.0/20160.0,-5603.0/2520.0,19087.0/60480.0],
			[16083.0/4480.0,-1152169.0/120960.0,242653.0/13440.0,-296053.0/13440.0,2102243.0/120960.0,-115747.0/13440.0,32863.0/13440.0,-5257.0/17280.0]
			]

	t,y = symbols('t y')

	aux = 0.0

	yn = yl[len(yl)-1]

	for j in range(len(coefM[order-1])):
		aux += h*coefM[order-1][j]*funct.subs({t:(t0-(h*j)),y:float(yl[len(yl)-j-1])})

	yn = yn + aux

	return yn

def EulerList(y0,t0,h,funct,order):
	yl = []

	yl.append(float(y0))

	t,y = symbols('t y')

	yn = float(y0)
	for i in range(order-1):
		yn = yn + h*funct.subs({t:t0,y:yn})
		t0 = t0 + h
		yl.append(yn)

	return yl

def EulerInvList(y0,t0,h,funct,order):
	yl = []

	yl.append(float(y0))

	t,y = symbols('t y')

	yn = float(y0)
	for i in range(order-1):
		yaux = EulerPd(yn,t0,h,funct)
		taux = t0 + h
		yn = yn + h*funct.subs({t:taux,y:yaux})
		t0 = t0 + h 
		yl.append(yn)
	return yl

def EulerAprimList(y0,t0,h,funct,order):
	yl = []

	yl.append(float(y0))

	t,y = symbols('t y')

	yn = float(y0)
	for i in range(order-1):
		k1 = funct.subs({t:t0,y:yn})
		k2 = funct.subs({t:(t0+h),y:yn+h*k1})
		yn = yn + float((h*(k1+k2))/2)
		t0 = t0 + h
		yl.append(yn)
	return yl

def RungeKuttaList(y0,t0,h,funct,order):
	yl = []

	yl.append(float(y0))

	t,y = symbols('t y')

	yn = float(y0)
	for i in range(order-1):
		k1 = funct.subs({t:t0,y:yn})
		k2 = funct.subs({t:(t0+h/2),y:yn+(h/2)*k1})
		k3 = funct.subs({t:(t0+h/2),y:yn+(h/2)*k2})
		k4 = funct.subs({t:(t0+h),y:yn+h*k3})
		yn = yn + float((h*(k1+2*k2+2*k3+k4))/6)
		t0 = t0 + h
		yl.append(yn)
	return yl


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

def Adam_Bashforth(yl,t0,h,qp,funct,order,straux=''):
	print('Metodo de Adam-Bashforth'+straux)
	print('y({}) = {}'.format(t0,yl[0]))
	print('h = {}'.format(h))


	coefM = [
			[1.0],
			[3.0/2.0,-1.0/2.0],
			[23.0/12.0,-4.0/3.0,5.0/12.0],
			[55.0/24.0,-59.0/24.0,37.0/24.0,-3.0/8.0],
			[1901.0/720.0,-1387.0/360.0,109.0/30.0,-637.0/360.0,251.0/720.0],
			[4277.0/1440.0,-2641.0/480.0,4991.0/720.0,-3649.0/720.0,959.0/480.0,-95.0/288.0],
			[198721.0/60480.0,-18637.0/2520.0,235183.0/20160.0,-10754.0/945.0,135713.0/20160.0,-5603.0/2520.0,19087.0/60480.0],
			[16083.0/4480.0,-1152169.0/120960.0,242653.0/13440.0,-296053.0/13440.0,2102243.0/120960.0,-115747.0/13440.0,32863.0/13440.0,-5257.0/17280.0]
			]

	t,y = symbols('t y')

	yn = float(yl[len(yl)-1])

	for i in range(order):
		print(i,' %.17f'%yl[i])
		t0 = t0 + h

	#t0 = t0 - h
	for i in range(order,qp+1,1):
		aux = 0.0

		for j in range(len(coefM[order-1])):
			aux += h*coefM[order-1][j]*funct.subs({t:(t0-(h*j)),y:float(yl[len(yl)-j-1])})

		yn = yn + aux
		t0 = t0 + h

		yl.append(yn)

		print(i,' %.17f'%yn)
	return 

def Adam_Multon(yl,t0,h,qp,funct,order,straux=''):
	print('Metodo de Adam-Multon'+straux)
	print('y({}) = {}'.format(t0,yl[0]))
	print('h = {}'.format(h))

	coefM = [
			[1.0],
			[1.0/2.0,1.0/2.0],
			[5.0/12.0,2.0/3.0,-1.0/12.0],
			[3.0/8.0,19.0/24.0,-5.0/24.0,1.0/24.0],
			[251.0/720.0,323.0/360.0,-11.0/30.0,53.0/360.0,-19.0/720.0],
			[95.0/288.0,1427.0/1440.0,-133.0/240.0,241.0/720.0,-173.0/1440.0,3.0/160.0],
			[19087.0/60480.0,2713.0/2520.0,-15487.0/20160.0,586.0/945.0,-6737.0/20160.0,263.0/2520.0,-863.0/60480.0],
			[5257.0/17280.0,139849.0/120960.0,-4511.0/4480.0,123133.0/120960.0,-88547.0/120960.0,1537.0/4480.0,-11351.0/120960.0,275.0/24192.0]
			]

	t,y = symbols('t y')


	yn = float(yl[len(yl)-1])

	for i in range(order):
		print(i,' %.17f'%yl[i])
		t0 = t0 + h

	#t0 = t0 - h
	for i in range(order,qp+1,1):
		aux = 0.0

		implict = h*coefM[order-1][0]*funct.subs({t:(t0+h),y:Adam_BashforthPd(yl,t0,h,funct,order)})

		aux += implict

		for j in range(1,len(coefM[order-1]),1):
			aux += h*coefM[order-1][j]*funct.subs({t:(t0-(h*j)),y:float(yl[len(yl)-j-1])})

		yn = yn + aux
		t0 = t0 + h

		yl.append(yn)

		print(i,' %.17f'%yn)
	return

def Formula_Inversa():
	return

def main():

	for inpt in sys.stdin:
		entry = str(inpt).split()

		method = entry[0]

		if(method == 'euler'):
			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])

			Euler(y0,t0,h,qp,funct)
			print()

		elif(method == 'euler_inverso'):
			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])

			Euler_Inverso(y0,t0,h,qp,funct)
			print()

		elif(method == 'euler_aprimorado'):
			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])

			Euler_Aprimorado(y0,t0,h,qp,funct)
			print()

		elif(method == 'runge_kutta'):
			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])

			Runge_Kutta(y0,t0,h,qp,funct)
			print()

		elif(method == 'adam_bashforth'):
			yl = []
			size = len(entry)

			for i in range(int(entry[size-1])):
				yl.append(entry[i+1])
				yl[i] = float(yl[i])


			t0 = float(entry[size-5])
			h = float(entry[size-4])
			qp = int(entry[size-3])
			funct = parse_expr(entry[size-2])
			order = int(entry[size-1])
			
			Adam_Bashforth(yl,t0,h,qp,funct,order)
			print()

		elif(method == 'adam_bashforth_by_euler'):

			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])
			order = int(entry[6])

			yl = EulerList(y0,t0,h,funct,order)
			Adam_Bashforth(yl,t0,h,qp,funct,order,' por Euler')
			print()

		elif(method == 'adam_bashforth_by_euler_inverso'):

			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])
			order = int(entry[6])

			yl = EulerInvList(y0,t0,h,funct,order)
			Adam_Bashforth(yl,t0,h,qp,funct,order,' por Euler Inverso')
			print()

		elif(method == 'adam_bashforth_by_euler_aprimorado'):

			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])
			order = int(entry[6])

			yl = EulerAprimList(y0,t0,h,funct,order)
			Adam_Bashforth(yl,t0,h,qp,funct,order,' por Euler Aprimorado')
			print()

		elif(method == 'adam_bashforth_by_runge_kutta'):

			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])
			order = int(entry[6])

			yl = RungeKuttaList(y0,t0,h,funct,order)
			Adam_Bashforth(yl,t0,h,qp,funct,order,' por Runge-Kutta ( ordem = {} )'.format(order))
			print()

		elif(method == 'adam_multon'):

			yl = []
			size = len(entry)

			for i in range(int(entry[size-1])):
				yl.append(entry[i+1])
				yl[i] = float(yl[i])


			t0 = float(entry[size-5])
			h = float(entry[size-4])
			qp = int(entry[size-3])
			funct = parse_expr(entry[size-2])
			order = int(entry[size-1])
			
			Adam_Multon(yl,t0,h,qp,funct,order)
			print()

		elif(method == 'adam_multon_by_euler'):

			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])
			order = int(entry[6])

			yl = EulerList(y0,t0,h,funct,order)
			Adam_Multon(yl,t0,h,qp,funct,order,' por Euler')
			print()

		elif(method == 'adam_multon_by_euler_inverso'):

			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])
			order = int(entry[6])

			yl = EulerInvList(y0,t0,h,funct,order)
			Adam_Multon(yl,t0,h,qp,funct,order,' por Euler Inverso')
			print()
		elif(method == 'adam_multon_by_euler_aprimorado'):

			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])
			order = int(entry[6])

			yl = EulerAprimList(y0,t0,h,funct,order)
			Adam_Multon(yl,t0,h,qp,funct,order,' por Euler Aprimorado')
			print()
		elif(method == 'adam_multon_by_runge_kutta'):

			y0 = float(entry[1])
			t0 = float(entry[2])
			h = float(entry[3])
			qp = int(entry[4])
			funct = parse_expr(entry[5])
			order = int(entry[6])

			yl = RungeKuttaList(y0,t0,h,funct,order)
			print(yl)
			Adam_Multon(yl,t0,h,qp,funct,order,' por Runge-Kutta ( ordem = {} )'.format(order))
			print()

		# elif(method == 'formula_inversa'):

		# elif(method == 'formula_inversa_by_euler'):

		# elif(method == 'formula_inversa_by_euler_inverso'):

		# elif(method == 'formula_inversa_by_euler_aprimorado'):

		# elif(method == 'formula_inversa_by_runge_kutta'):

	return 0

if __name__ == "__main__":
	main()