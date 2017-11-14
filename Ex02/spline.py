import math

def Runge(x):
	return (1/(1+25*x**2))

def Spline(xx,fun,xmin,xmax,steps,boundary="align",lderivative=0,rderivative=0):
	x=[]
	N=steps
	h=(xmax-xmin)/steps
	# 写好采样点
	xtmp=xmin
	for n in range(N+1):
		x.append(xtmp)
		xtmp+=h
	y=[fun(xi) for xi in x]
	
	def Thomas(aa,bb,cc,dd):
		
		#LU分解
		a=[aa[0]]
		b=[0]
		for i in range(1,len(aa)):
			b.append(bb[i]/a[-1])
			a.append(aa[i]-b[-1]*cc[i-1])
		
		#反代求解
		Z=[dd.pop()/a.pop()]
		for i in range(1,len(aa)):
			j=len(aa)-i
			Z[:0]=[(dd.pop()-cc.pop()*Z[0])/a.pop()]
			
		X=[Z[0]]
		for i in range(1,len(aa)):
			X.append(Z[i]-b[i]*X[-1])
		return X
		
	if boundary=="align":
		aa=[2*h/3 for i in range(N+1)]
		aa[0]=h/3
		aa[-1]=h/3
		bb=[h/6 for i in range(N+1)]
		cc=bb
		dd=[(y[j+1]+y[j-1]-2*y[j])/h for j in range(1,N)]
		dd[:0]=[(y[1]-y[0])/h-lderivative]
		dd.append((rderivative-(y[-1]-y[-2])/h))
		M=Thomas(aa,bb,cc,dd)
		#至此各个节点的矩都算出来了，可以计算函数值了。
		j=math.floor((xx-xmin)/h)
		return (M[j]*(x[j+1]-xx)**3/(6*h)+M[j+1]*(xx-x[j])**3/(6*h)+((y[j+1]-y[j])/h-(M[j+1]-M[j])*h/6)*(xx-x[j])+y[j]-M[j]*h**2/6)
	
	elif boundary=="periodical":
		#周期边界条件的话只需要解n-1阶的方程
		x.pop()
		
		aa=[2*h/3 for i in range(N)]
		bb=[h/6 for i in range(N)]
		cc=bb
		dd=[(y[j+1]+y[j-1]-2*y[j])/h for j in range(1,N)]
		y.pop()
		dd[:0]=[(y[1]+y[-1]-2*y[0])/h]
		
		aa[-1]+=bb[0]*cc[-1]/aa[0]
		aa[0]+=aa[0]
		M1=Thomas(aa,bb,cc,dd)
		
		bb=[h/6 for i in range(N)]
		cc=bb
		ee=[0 for i in range(N)]
		ee[0]=-2*h/3
		ee[-1]=cc[-1]
		
		q=Thomas(aa,bb,cc,ee)
		M=[M1[i]-q[i]*(M1[0]-M1[-1]/4)/(1+q[0]-q[-1]/4) for i in range(N)]
		j=math.floor((xx-xmin)/h)
		x.append(xmax)
		y.append(y[0])
		M.append(M[0])
		return (M[j]*(x[j+1]-xx)**3/(6*h)+M[j+1]*(xx-x[j])**3/(6*h)+((y[j+1]-y[j])/h-(M[j+1]-M[j])*h/6)*(xx-x[j])+y[j]-M[j]*h**2/6)
	
f=open("spline.txt",'w')
wrt="Runge starts...\n"

Xlist=[-1+i/20 for i in range(40)]
for x in Xlist:
	wrt+=str(x)+"	&"+str(Runge(x))+"	&"+str(Spline(x,Runge,-1,1,20,"align",25/338,-25/338))+"\n\\\\\n" 
f.write(wrt)

wrt="\n\nCardioid starts...\n"

def Cardx(t):
	return(math.cos(t)*(1-math.cos(t)))
def Cardy(t):
	return(math.sin(t)*(1-math.cos(t)))
	
Xlist=[i*math.pi/20 for i in range(40)]
for x in Xlist:
	wrt+=str(x)+"	&"+str(Cardx(x))+"	&"+str(Spline(x,Cardx,0,2*math.pi,8,"periodical"))+"	&"+str(Cardy(x))+"	&"+str(Spline(x,Cardy,0,2*math.pi,8,"periodical"))+"\n\\\\\n" 
f.write(wrt)



f.close()