def xrate(x,y):
	return(2/3*x-4/3*x*y)
def yrate(x,y):
	return(x*y-y)

def DSolve(dxdt,dydt,xi,yi,err=1e-11,tmax=10):
	def RK(h,f,g,x,y):
		[k1,l1]=[h*f(x,y),h*g(x,y)]
		[k2,l2]=[h*f(x+k1/2,y+l1/2),h*g(x+k1/2,y+l1/2)]
		[k3,l3]=[h*f(x+k2/2,y+l2/2),h*g(x+k2/2,y+l2/2)]
		[k4,l4]=[h*f(x+k3,y+l3),h*g(x+k3,y+l3)]
		return([x+k1/6+k2/3+k3/3+k4/6,y+l1/6+l2/3+l3/3+l4/6])
	sol = [[0,xi,yi]]	# 这将作为解的输出，每一行是一个数据点，三个分量分别是t、x、y
	h = 0.1	# 在违法的边缘试探.jpg
	while sol[-1][0] < tmax:
		[x1,y1] = RK(h,dxdt,dydt,sol[-1][1],sol[-1][2])
		[x2,y2] = RK(h,dxdt,dydt,x1,y1)
		[x3,y3] = RK(2*h,dxdt,dydt,sol[-1][1],sol[-1][2])
		hcorrection = err/max(abs(x3-x2),abs(y3-y2))
		if hcorrection>=1:
			sol.append([sol[-1][0]+h,x1,y1])
		h = h*hcorrection**0.2
	return sol
# import math
# s = DSolve(xrate,yrate,0.8,0.8,err=1e-11,tmax=2)
# for a in s:
	# print(math.log(a[1])-a[1]+2/3*math.log(a[2])-4/3*a[2])

# 画图画图
import numpy
import matplotlib
matplotlib.use('PS')
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

for i in [0.8,1.0,1.2,1.4,1.6]:
	s= DSolve(xrate,yrate,i,i,err=1e-11,tmax=10)
	xlist = []; ylist = [];
	for a in s:
		xlist.append(a[1]);
		ylist.append(a[2]);
	cm = plt.cm.get_cmap('brg')
	tick = [[],[],[]]
	for tock2 in range(int(4*i+13)):
		tock = tock2/2
		for a in s:
			if a[0]>tock:
				for j in range(3):
					tick[j].append(a[j])
				break
	
	ax.plot(xlist,ylist,linewidth=.2)
	global b
	b = ax.scatter(tick[1],tick[2],c=tick[0],vmin=0,vmax=9,s=12,marker="*",cmap=cm)
cb = plt.colorbar(b)
cb.set_label('t')
plt.xlabel('x')
plt.ylabel('y')
fig.savefig("RK.eps",bbox_inches='tight')