import math
def f(x,y):
	return(2*math.pi**2*math.sin(math.pi*x)*math.sin(math.pi*y))

	
def Jacobi (N,f):
	h = 1/(N+1)
	#	起始迭代
	z = [[0 for i in range(N+2)]for j in range(N+2)]
	loops = 0
	while True:
		loops += 1
		#print(loops)
		z1=[]
		for i in range(len(z)):
			z1.append(z[i][:])
		for i in range(1,N+1):
			for j in range(1,N+1):
				z[i][j]=(z1[i][j+1]+z1[i][j-1]+z1[i-1][j]+z1[i+1][j])/4 + f(i*h,j*h)*h**2/4
		flag = True
		for i in range(1,N+1):
			for j in range(1,N+1):
				resid = abs(((z[i][j+1]+z[i][j-1]+z[i-1][j]+z[i+1][j])-4*z[i][j])/h**2+f(i*h,j*h))
				if resid > 0.0001:
					flag = False
					#print(resid)
					break
			if flag ==False:
				break
		if flag == True:
			break
	print (loops)
	return z

a = Jacobi(10,f)
a = Jacobi(20,f)
a = Jacobi(50,f)

		
				
	