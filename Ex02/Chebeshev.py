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

writ=""

samplist=[]
for k in range(N):
	samplist.append(math.cos(math.pi*(k+.5)/N))
xlist=[samplist[0]]
for k in range(N-1):
	xlist.extend([(samplist[k]+samplist[k+1])/2,samplist[k+1]])

for x in xlist:
	writ=writ+str(x)+"&	"+str(Runge(x))+"&	"+str(Chebeshev(x,Runge,N))+"\\\\\n"

f=open("Chebeshev.txt",'w')
f.write(writ)
f.close()	

