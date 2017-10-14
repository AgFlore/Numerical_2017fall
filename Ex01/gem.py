print("Applying Gauss Elimination Method to Hibert Matrix.")
N = int( input("Input the N you wish to stop at..."))
hilbert = [[1/(i+j-1) for i in range(1,N+1)] for j in range(1,N+1)]
f = open('GEMOutput.txt','w')
writ='Solving H_n x=b using GEM, where H_n is the Hilbert Matrix, and b=(1,1,...,1).'
f.write(writ)
triaflag=input('press y and enter if you want to have the reduced triangular matrix of H_n printed in the file. If you don\'t, enter right away...')
for m in range (2,N+1):
	a= [[1/(i+j-1) for i in range(1,m+1)]for j in range(1,m+1)]
	for i in range (0,m):
		a[i].append(1)#a便是m*(m+1)的增广矩阵
	#化上三角。i是化上三角所进行到的前线位置
	for i in range (0,m-1):
		#支点选择
		pivotinglist = [[abs(a[j][i]),j] for j in range (i,m)]
		pivot = max(pivotinglist,key=lambda item: item[0])[1]
		a[pivot],a[i]=a[i],a[pivot]
		#作初等变换将下三角项置零
		for j in range (i+1,m):
			coef=a[j][i]/a[i][i]
			for k in range (i,m+1):#k是列号，所以上限到m
				a[j][k]-=a[i][k]*coef
	#偷偷看一眼上三角成果如何
	if (triaflag=='y'):
		writ='\n\n\n\nFor n='+str(m)+', we have eliminated the matrix to a upper-triangular one. We shall print its transpostition here:'
		for i in range (0,m):
			writ+='\n'
			writ+=str(a[i])
		f.write(writ)
	#从上三角求出解
	x=[a[i][m]for i in range (0,m)]
	for i in range (1,m+1):
		for j in range (1,i):
			x[m-i]-=a[m-i][m-j]*x[m-j]
		x[m-i]=x[m-i]/a[m-i][m-i]
	writ='\n\nAnd we have solved x out. x='+str(x)+'is the solution to H_'+str(m)+'x=b.\n'
	f.write(writ)
print ('The results have been written to GEMOutput.txt. The program shall stop now.')
f.close()	
