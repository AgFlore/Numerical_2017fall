#Neville算法

#输入的Runge函数。
def inputfunc(x):
	y=1/(1+25*x**2)
	return y

#插值所得函数
def Neville(x,fun,sampxlist):
	nevarray=[]
	for x0 in sampxlist:
		nevarray[:0]=[fun(x0)]
		for j in range (1,len(nevarray)):
			nevarray[j]=((x-sampxlist[len(nevarray)-j-1])*nevarray[j-1]-(x-x0)*nevarray[j])/(x0-sampxlist[len(nevarray)-j-1])
	return nevarray.pop()
	
			
#计算插值节点	
Min=-1
Max=1
Step=0.1
sampxlist=[Min]
while (sampxlist[-1]<Max):
	sampxlist.append(sampxlist[-1]+Step)
#应老师要求，每两个节点中插入一个点来比较插值多项式与原函数的值。
Density=2
xlist=[Min]
while (xlist[-1]<Max):
	xlist.append(xlist[-1]+Step/Density)

f=open("Neville.txt",'w')

writ="用多项式在[-1,1]上内插Runge函数1/(1+25x^2)。\n\n$x$&$R(x)=(1+25x^2)^{-1}$	&多项式$P_{20}(x)$	&$\\lvert P_{20}(x)-R(x)\rvert$\\\\\n\\hline"
for x in xlist:
	writ+="\\num{"+str(x)+"}&	\\num{"+str(inputfunc(x))+"}&	\\num{"+str(Neville(x,inputfunc,sampxlist))+"}&	\\num{"+str(abs(Neville(x,inputfunc,sampxlist)-inputfunc(x)))+"}\\\\\n"

f.write(writ)
f.close()	


#开始画图
#所以接下来我import numpy就不算犯规了哟～
#助教姐姐如果没有装某些包导致跑不了的话，把下面的代码直接注释掉就好啦～

import numpy
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt

xpoints=numpy.arange(-1.1,1.1,0.001)
xlist=xpoints.tolist()

rlist=[inputfunc(x) for x in xlist]
plist=[Neville(x,inputfunc,sampxlist) for x in xlist]

fig, ax = plt.subplots()
ax.plot(xlist,plist,color="red",linewidth=.2)
ax.plot(xlist,rlist,color="blue",linewidth=.2)


ax.set(xlim=(-1.05,1.05),ylim=(-1,2))
ax.grid(linewidth=0.2)


fig.savefig("Neville.eps")	#保存为Neville.eps