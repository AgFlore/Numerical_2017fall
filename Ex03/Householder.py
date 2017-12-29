import random
import timeit
random.seed()
def Householder (A):
	#	输入的方阵A，a[i][j]即是第i行第j列
	N = len(A)
	Q = [[0 for i in range(N)] for j in range(N)];
	for i in range(N):
		Q[i][i]=1
	R = A[:]
	for k in range(N):
		# 作为反射标的的单位矢量
		rotor = [R[i][k] for i in range(k,N)]
		# 算rotor的模长
		M = 0
		for r in rotor:
			M += r**2
		M = M**0.5
		rotor[0]+=M
		M = 0
		for r in rotor:
			M += r**2
		M = M**0.5
		for i in range(len(rotor)):
			rotor[i]=rotor[i]/M
		# 对Q和R分别作旋转。注意R中倒映的是列而Q中倒映的是行
		for i in range(k,N):
			roverlap = qoverlap = 0
			for j in range(k,N):
				roverlap += rotor[j-k]*R[j][i]
				qoverlap += rotor[j-k]*Q[i][j]
			for j in range(k,N):
				R[j][i] = 2*roverlap*rotor[j-k]-R[j][i]
				Q[i][j] = 2*qoverlap*rotor[j-k]-Q[i][j]
	return(Q+R)

#生成
A = [[2-2*random.random() for i in range(6)]for j in range(6)]
print(A)
print(Householder(A))
print(timeit.timeit('C=Householder([[2-2*random.random() for i in range(6)]for j in range(6)])',setup="from __main__ import Householder,A,random",number=10000))