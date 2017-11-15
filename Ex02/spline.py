import math

def Runge(x):
	return (1/(1+25*x**2))

def Spline(xx,fun,xmin,xmax,steps,boundary="align",lderivative=0,rderivative=0):	# xx为样条函数的（哑）自变量，fun为被内插的函数。xmin和xmax为内插区间的端点，steps为分隔区间个数。boundary目前开发了两种选择："align"代表边界条件取为样条函数在端点的一阶导数值等于原函数，"periodical"代表边界条件取为周期边界条件。lderivative和rderivative两个参数适用于boundary='align'时，为原函数在区间左右端点的一阶导数值，作为样条函数的边界条件。
	x=[]	#节点的x存在这里
	N=steps
	h=(xmax-xmin)/steps	#每个区间长度
	xtmp=xmin
	for n in range(N+1):
		x.append(xtmp)
		xtmp+=h
	y=[fun(xi) for xi in x]	#各节点的y存放于此
	
	def Thomas(aa,bb,cc,dd):	#三对角矩阵Ax=dd求解。aa为主对角线，bb[1:]与cc分别为下与上副对角线。注意cc的长度比aa和bb短1。
		
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
		#周期边界条件的话只需要解n-1阶的方程组。但是这时的系数矩阵不是三对角矩阵，需要额外处理一下。
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
	wrt+="\\num{"+str(x)+"}	&\\num{"+str(Runge(x))+"}	&\\num{"+str(Spline(x,Runge,-1,1,20,"align",25/338,-25/338))+"}	&\\num{"+str(abs(Spline(x,Runge,-1,1,20,"align",25/338,-25/338)-Runge(x)))+"}\n\\\\\n" 
f.write(wrt)

wrt="\n\nCardioid starts...\n"

def Cardx(t):
	return(math.cos(t)*(1-math.cos(t)))
def Cardy(t):
	return(math.sin(t)*(1-math.cos(t)))
	
Xlist=[i*math.pi/20 for i in range(40)]
for x in Xlist:
	wrt+="\\num{"+str(x)+"}	&\\num{"+str(Cardx(x))+"}	&\\num{"+str(Spline(x,Cardx,0,2*math.pi,8,"periodical"))+"}	&\\num{"+str(Cardy(x))+"}	&\\num{"+str(Spline(x,Cardy,0,2*math.pi,8,"periodical"))+"}\n\\\\\n" 
f.write(wrt)

f.close()


#开始画图
#所以接下来我import numpy就不算犯规了哟～
#助教姐姐如果没有装matplotlib的话把下面的代码直接注释掉就好啦～

import numpy
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt

xpoints=numpy.arange(-1,0.999,0.001)
xlist=xpoints.tolist()

rlist=[Runge(x) for x in xlist]
plist=[Spline(x,Runge,-1,1,20,"align",25/338,-25/338) for x in xlist]
rpoints=numpy.asarray(rlist)
ppoints=numpy.asarray(plist)

fig, ax = plt.subplots()
ax.plot(xpoints,ppoints,linewidth=0.2)
ax.plot(xpoints,rpoints,linewidth=0.2)


ax.set(xlim=(-1,1),ylim=(-0.2,1.2))
ax.grid(linewidth=0.2)

fig.savefig("spline.eps",dpi=300)

#画第二题的图
tpoints=numpy.arange(0,6.283,0.001)
tlist=tpoints.tolist()
#标准的心脏线
xlist=[Cardx(t) for t in tlist]
ylist=[Cardy(t) for t in tlist]

#样条(sp)的心脏线
xlisp=[Spline(t,Cardx,0,2*math.pi,8,"periodical") for t in tlist]
ylisp=[Spline(t,Cardy,0,2*math.pi,8,"periodical") for t in tlist]

fig, ax = plt.subplots()
ax.plot(xlist,ylist,linewidth=0.2,color="blue")
ax.plot(xlisp,ylisp,linewidth=0.2,color="red")

tlist=[i*math.pi/4 for i in range(8)]
xlist=[Cardx(t) for t in tlist]
ylist=[Cardy(t) for t in tlist]
ax.scatter(xlist,ylist,s=12,marker="*",color="purple")

ax.set(xlim=(-2.5,0.5),ylim=(-1.5,1.5))
ax.grid(linewidth=0.2)

fig.savefig("cardioid.eps",dpi=300)