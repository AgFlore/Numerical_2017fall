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
print (zaoxinbi)