import math

def f(x,y):
	return(2*math.pi**2*math.sin(math.pi*x)*math.sin(math.pi*y))
	
def SOR(N,f):
	h = 1/(N+1)
	#	弛豫参数
	rH = math.cos(math.pi/(N+1))**2
	omega = [1,1/(1-rH/2)]
	#	计数君
	loops = 0; rmax = 0
	z = [[0 for i in range(N+2)]for j in range(N+2)]
	b = [[f(i*h,j*h) for j in range(N+2)] for i in range(N+2)]
	while True:
		rmax = 0; loops += 1
		c = [[b[i][j]*h**2+z[i][j+1]+z[i+1][j]+4*(1/omega[0]-1)*z[i][j] for j in range(N+1)] for i in range(N+1)]
		for i in range(1,N+1):
			for j in range(1,N+1):
				z[i][j] = (c[i][j]+z[i][j-1]+z[i-1][j])/4*omega[0]
		
		omega=[1/(1-rH*omega[i]/4) for i in range(2)]
		flag = True
		for i in range(1,N+1):
			for j in range(1,N+1):
				resid = abs((z[i][j+1]+z[i][j-1]+z[i-1][j]+z[i+1][j]-4*z[i][j])/h**2+f(i*h,j*h))
				rmax = max(rmax,resid)
				if resid > 0.0001:
					flag = False
					#print(resid)
					break
			if flag ==False:
				break
		if flag == True:
			break
	print ([N,loops,rmax])
	return z

for i in range(6):
	a = SOR(10*(i+1),f)