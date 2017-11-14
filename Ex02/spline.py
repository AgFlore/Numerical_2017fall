def Runge(x):
	return (1/(1+25*x**2))

def Spline(xx,fun,boundary="align",xmin,xmax,steps,lderivative=0,rderivative=0):
	x=[]
	N=steps
	h=(xmax-xmin)/steps
	# 写好采样点
	xtmp=xmin
	for n in range(N+1):
		x.append(xtmp)
		xtmp+=h
	y=[]
	for xi in x:
		y.append(fun(xi))
	
	if boundary=="align":
		#LU分解
		a=[h/3]
		b=[0]
		for i in range(1,N):
			b.append(h/6/a[-1])
			a.append(2h/3-b[-1]*h/6)
		b.append(h/6/a[-1])
		a.append(h/3-b[-1]*h/6)
		#反代求解
		z=[(rderivative-(y[-1]-y[-2])/h)/a.pop()]
		for i in range(1,N):
			j=N-i
			Z[:0]=[((y[j+1]+y[j-1]-2*y[j])/h-Z[0]*h/6)/a.pop()]
		Z[:0]=[((y[1]-y[0])/h-lderivative-Z[0]*h/6)/a.pop()]
		M[0]=Z[0]
		for i in range(1,N+1):
			M.append(Z[i]-b[i]*M[-1])
		#至此各个节点的矩都算出来了，可以计算函数值了。
		j=math.floor((x-xmin)/h)
		return (M[j]*(x[j+1]-xx)**3/(6*h)+M[j+1]*(xx-x[j])**3/(6*h)+((y[j+1]-y[j])/h-(M[j+1]-M[j])*h/6)*(xx-x[j])+y[j]-M[j]*h**2/6)
	
	elif boundary=="periodical":
		#周期边界条件的话只需要解n-1阶的方程
		x.pop()
		y.pop()
		