import random
import timeit
random.seed()
def Givens (A):
	#	输入的方阵A，a[i][j]即是第i行第j列
	N = len(A)
	Q = [[0 for i in range(N)] for j in range(N)];
	for i in range(N):
		Q[i][i]=1
	R = A[:]
	#	旋转操作。嗯就是旋转操作。
	for j in range(N):
		for i in range(j+1,N):
			M = (R[j][j]**2+R[i][j]**2)**0.5
			[Q[j][j],Q[j][i],Q[i][j],Q[i][i]] = [(R[j][j]*Q[j][j]+R[i][j]*Q[i][j])/M,(R[j][j]*Q[j][i]-R[i][j]*Q[i][i])/M,(R[i][j]*Q[j][j]+R[j][j]*Q[i][j])/M,(-R[i][j]*Q[j][i]+R[j][j]*Q[i][i])/M]
			[R[j][j],R[j][i],R[i][j],R[i][i]] = [(R[j][j]*R[j][j]+R[i][j]*R[i][j])/M,(R[j][j]*R[j][i]+R[i][j]*R[i][i])/M,(-R[i][j]*R[j][j]+R[j][j]*R[i][j])/M,(-R[i][j]*R[j][i]+R[j][j]*R[i][i])/M]
	
	return(Q+R)

#生成
A = [[2-2*random.random() for i in range(6)]for j in range(6)]
print(A)
print(Givens(A))
print(timeit.timeit('C=Givens([[2-2*random.random() for i in range(6)]for j in range(6)])',setup="from __main__ import Givens,A,random",number=10000))
print(timeit.timeit('C=Givens([[2-2*random.random() for i in range(12)]for j in range(12)])',setup="from __main__ import Givens,A,random",number=10000))
print(timeit.timeit('C=Givens([[2-2*random.random() for i in range(18)]for j in range(18)])',setup="from __main__ import Givens,A,random",number=10000))