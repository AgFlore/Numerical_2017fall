import math
import random
random.seed()
#	读取文件并存入数组
datafile=open('correlation-function.dat')
[samples,timeslices]=datafile.readline().split()
samples=int(samples)
timeslices=int(timeslices)
datraw=[[] for i in range(timeslices)]
for N in range(samples):
	for M in range(timeslices):
		[t,Re,Im]=datafile.readline().split()
		t=int(t)
		Re=float(Re)
		Im=float(Im)
		datraw[t].append(Re)
#	现在的datraw的第一个entry是0～63的timeslice，第二个entry是0～199的序号
#	作对称化
datmat=[datraw[i] for i in range(int(timeslices/2+1))]
for Nt in range(1,int(timeslices/2)):
	for i in range(samples):
		datmat[Nt][i]=(datraw[Nt][i]+datraw[timeslices-Nt][i])/2
#	现在的datmat的第一个entry是0～32的timeslice，第二个entry是0～199的序号

#	估计中心值
corcent = [sum(datmat[i])/samples for i in range(len(datmat))]
#	估计方差
corresi = [[(datmat[i][j]-corcent[i])**2 for j in range(samples)]for i in range(len(corcent))]
corsvar = [sum(corresi[i])/(samples-1) for i in range(len(corcent))]
zaoxinbi = [corsvar[i]**0.5/abs(corcent[i]) for i in range(len(corcent))]

#此处需要输出信噪比数值

#	用ratio法计算有效质量
meffrc=[math.log(corcent[i]/corcent[i+1]) for i in range(len(corcent)-1)]

#	Jackknife法估算meffrc的中心值
meffrjc=[]
for i in range(len(meffrc)):
	c=0
	for j in range(samples):
		datata=[datmat[i][:],datmat[i+1][:]]
		datata[0].pop(j)
		datata[1].pop(j)
		c+=math.log(sum(datata[0])/sum(datata[1]))
	meffrjc.append(c/samples)



#	Jackknife法估算meffrc的误差
meffrj=[]
for i in range(len(meffrc)):
	resid=0
	for j in range(samples):
		datata=[datmat[i][:],datmat[i+1][:]]
		datata[0].pop(j)
		datata[1].pop(j)
		meffff=math.log(sum(datata[0])/sum(datata[1]))
		resid+=(meffff-meffrjc[i])**2
	meffrj.append(((1-1/samples)*resid)**0.5)

#此处可以输出meffrc的误差meffrj


#	最小卡方拟合
chi2perdof = 1000000000
mint = maxt = mfitted = 0
for min in range(len(meffrc)-4):
	for max in range(min+3,len(meffrc)):
		cinterv = meffrc[min:max+1]
		vinterv = meffrj[min:max+1]
		temp1 = [cinterv[i]**2/vinterv[i]**2 for i in range(len(cinterv))]
		temp2 = [cinterv[i]/vinterv[i]**2 for i in range(len(cinterv))]
		mfit = sum(temp1)/sum(temp2)
		temp3 = [((cinterv[i]-mfit)/vinterv[i])**2 for i in range(len(cinterv))]
		chiperd = sum(temp3)/len(cinterv)
		if (chiperd < chi2perdof):
			chi2perdof = chiperd
			mfitted = mfit
			mint = min
			maxt = max

print ([mint,maxt,mfitted,chi2perdof])
print ([meffrc,meffrj])


#	用acosh法计算有效质量
meffcc=[math.acosh((corcent[i]+corcent[i+2])/(2*corcent[i+1])) for i in range(len(corcent)-2)]

#	Jackknife法估算meffcc的中心值
meffcjc=[]
for i in range(len(meffcc)):
	c=0
	for j in range(samples):
		datata=[datmat[i][:],datmat[i+1][:],datmat[i+2][:]]
		datata[0].pop(j)
		datata[1].pop(j)
		datata[2].pop(j)
		c += math.acosh((sum(datata[0])+sum(datata[2]))/(2*sum(datata[1])))
	meffcjc.append(c/samples)


#	Jackknife法估算meffcc的误差
meffcj=[]
for i in range(len(meffcc)):
	resid=0
	for j in range(samples):
		datata=[datmat[i][:],datmat[i+1][:],datmat[i+2][:]]
		datata[0].pop(j)
		datata[1].pop(j)
		datata[2].pop(j)
		meffff = math.acosh((sum(datata[0])+sum(datata[2]))/(2*sum(datata[1])))
		resid+=(meffff-meffcjc[i])**2
	meffcj.append(((1-1/samples)*resid)**0.5)

#此处可以输出meffcc的误差meffcj


#	最小卡方拟合
chi2perdof = 1000000000
mint = maxt = mfitted = 0
for min in range(len(meffcc)-4):
	for max in range(min+3,len(meffcc)):
		cinterv = meffcc[min:max+1]
		vinterv = meffcj[min:max+1]
		temp1 = [cinterv[i]**2/vinterv[i]**2 for i in range(len(cinterv))]
		temp2 = [cinterv[i]/vinterv[i]**2 for i in range(len(cinterv))]
		mfit = sum(temp1)/sum(temp2)
		temp3 = [((cinterv[i]-mfit)/vinterv[i])**2 for i in range(len(cinterv))]
		chiperd = sum(temp3)/len(cinterv)
		if (chiperd < chi2perdof):
			chi2perdof = chiperd
			mfitted = mfit
			mint = min
			maxt = max

print ([mint,maxt,mfitted,chi2perdof])
print ([meffcc,meffcj])


# Boot strap for covariance matrix
def CovMat (ca,cb):
	s=[(ca[i]-sum(ca)/len(ca))*(cb[i]-sum(cb)/len(cb)) for i in range(len(ca))]
	return (sum(s)/(len(ca)-1))

def BSSample (ta,tb,N=2000):
	bslist=[]
	for i in range(N):
		ca = []; cb =[]
		for whatever in range(samples):
			j = random.randrange(samples)
			ca.append(datmat[ta][j])
			cb.append(datmat[tb][j])
		bslist.append(CovMat(ca,cb)/(CovMat(ca,ca)*CovMat(cb,cb))**0.5)
	return bslist

s = BSSample(3,4)
c34 = sum(s)/2000
r34 = [(s[i]-c34)**2 for i in range (2000)]
v34 = (sum(r34)/1999)**0.5
s = BSSample(3,5)
c35 = sum(s)/2000
r35 = [(s[i]-c35)**2 for i in range (2000)]
v35 = (sum(r35)/1999)**0.5

print([c34,v34])
print([c35,v35])





