import math
def Chebeshev (x,Fun,N):
	import math
	
	#标准节点
	samplist=[]
	for k in range(N):
		samplist.append(math.cos(math.pi*(k+.5)/N)) 
	
	#计算系数
	c=[]
	for m in range(N):
		s=0
		for k,q in enumerate(samplist):
			s+=Fun(q)*math.cos(m*math.pi*(k+.5)/N)
		c.append(s*2/N)
	
	#秦九韶算法
	b=[0,0]
	for k in reversed(range(1,N)):
		b[:0]=[c[k]+2*x*b[-2]-b[-1]]
		b.pop()
	return (c[0]/2+x*b[0]-b[1])

def Runge (x):
	return (1/(1+25*x**2))

N=20


#按老师要求需要计算的点包括内插的节点及相邻节点的中点
samplist=[]
for k in range(N):
	samplist.append(math.cos(math.pi*(k+.5)/N))
xlist=[samplist[0]]
for k in range(N-1):
	xlist.extend([(samplist[k]+samplist[k+1])/2,samplist[k+1]])

#输出
writ=""
for x in xlist:
	writ=writ+"\\num{"+str(x)+"}&	\\num{"+str(Runge(x))+"}&	\\num{"+str(Chebeshev(x,Runge,N))+"}&	\\num{"+str(abs(Chebeshev(x,Runge,N)-Runge(x)))+"}\\\\\n"

f=open("Chebeshev.txt",'w')
f.write(writ)
f.close()	

#开始画图
#所以接下来我import numpy就不算犯规了哟～
#助教姐姐如果没有装matplotlib的话把下面的代码直接注释掉就好啦～

import numpy
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt

xpoints=numpy.arange(-1.5,1.5,0.001)
xlist=xpoints.tolist()

rlist=[Runge(x) for x in xlist]
plist=[Chebeshev(x,Runge,N) for x in xlist]


fig, ax = plt.subplots()
ax.plot(xlist,plist,color="red",linewidth=.2)
ax.plot(xlist,plist,color="blue",linewidth=.2)


ax.set(xlim=(-1.05,1.05),ylim=(-0.2,1.2))
ax.grid(linewidth=0.2)

fig.savefig("Chebeshev.eps")