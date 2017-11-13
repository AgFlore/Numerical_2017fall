import math

#先封装一个Zeta函数。qs即为q^2
def Zeta (qs,N=6,steps=1000):
	#第一项
	t1=0
	for x in range(-N,N+1):
		for y in range(-N,N+1):
			for z in range (-N,N+1):
				t1+=math.exp(qs-x**2-y**2-z**2)/(x**2+y**2+z**2-qs)
	t1=t1/math.sqrt(4*math.pi)
	
	#第二项
	t2=-math.pi
	
	#第三项，直接泰勒级数算
	t3=0
	for n in range(1,30):
		t3+=qs**n/math.factorial(n)/(n-0.5)
	t3=t3*math.pi/2
	
	#第四项
	def sum4 (t):
		s=0
		for x in range(-N,N+1):
			for y in range(-N,N+1):
				for z in range (-N,N+1):
					s+=math.exp(-math.pi**2*(x**2+y**2+z**2)/t)
		s-=1
		return s
	t4 = 0
	#Cotes积分公式。尚待优化成Romberg外推emmmm
	for i in range(steps):
		x1=(i+0.5)/steps
		x2=(i+1)/steps
		t4+=(2*sum4(x1)*math.exp(x1*qs)/x1**(1.5)/(3*steps)+sum4(x2)*math.exp(x2*qs)/x2**(1.5)/(3*steps)) #t为零时被积函数为零。每个小Z区间的1+4+1合并为4+2计算，最后再扣掉多加的t=1
	t4-=sum4(1)*math.exp(qs)/(6*steps)
	t4=t4*math.sqrt(math.pi)/2
	
	return (t1+t2+t3+t4)	

f=open("Zeta.txt",'w')
f.write('Output starts...')


#做表格可视化
def Tableit (xlist,N):
	wrt="\n$q^2$"
	for n in range(1,N+1):
		wrt+=("&$N="+str(n)+"$")
	wrt+="\\\\\n"
	for x in xlist:
		wrt+=(str(x))
		for n in range(1,N+1):
			wrt+=("&"+str(Zeta(qs=x,N=n)))
		wrt+="\\\\\n"
	f.write(wrt)
	return 0

#做第一小问
xlist=[.0001,.001,.01,.05,.2,.5,.8,.99,1.01,1.3,1.7,1.99,2.01,2.3,2.7,2.9]
Tableit(xlist,7)


#做第二小问。对分法解方程。要求fun尽量单调性好，并且在xmin和xmax上异号
def Bisect (fun,xmin,xmax,threshold=1e-10,maxstep=80):
	x1=xmin
	x2=xmax
	def flg (y,thres=threshold):
		if y>thres:
			return 1
		elif y<-thres:
			return -1
		else:
			return 0
	f1=flg(fun(x1))
	f2=flg(fun(x2))
	for n in range(maxstep):
		print('Round '+str(n)+' of bisecting has started. The search range is from '+str(x1)+' to '+str(x2)+'...\n')
		if f1==0:
			return x1
		elif f2==0:
			return x2
		else:
			f3=flg(fun((x1+x2)/2))
			if(f3==0):
				return (x1+x2)/2
			elif(f3*f1>0):
				x1=(x1+x2)/2
				f1=flg(fun(x1))
			else:
				x2=(x1+x2)/2
				f2=flg(fun(x2))
	return float('nan')

#待解的函数
def targetfunc(qs):
	return (math.pi**1.5*(1+qs/4)-Zeta(qs,N=5))

	
qssolu=Bisect(targetfunc,.01,.99,1e-15)

if math.isnan(qssolu):
	f.write("Root not found within steps given!")
else:
	f.write(str(qssolu))
	
f.close()