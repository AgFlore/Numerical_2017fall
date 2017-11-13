#Neville算法

#输入的函数。函数形式可以随时修改。
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
	
			
	
Min=-1
Max=1
Step=0.1
Density=2

sampxlist=[Min]
while (sampxlist[-1]<Max):
	sampxlist.append(sampxlist[-1]+Step)

xlist=[Min]
while (xlist[-1]<Max):
	xlist.append(xlist[-1]+Step/Density)

print (xlist)	


writ="writestart...\n"  #提示语待修改

for x in xlist:
	writ=writ+str(x)+"&	"+str(inputfunc(x))+"&	"+str(Neville(x,inputfunc,sampxlist))+"\\\\\n"

f=open("Neville.txt",'w')
f.write(writ)
f.close()	
