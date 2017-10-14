print("Cholesky decomposition of Hilbert Matrix.")
N = int( input("Input the N you wish to stop at..."))
hilbert = [[1/(i+j-1) for i in range(1,N+1)] for j in range(1,N+1)]
f = open('choleskyOutput.txt','w')

chole = [[0 for i in range(1,N+1)]for j in range(1,N+1)]
chole[0][0]= hilbert[0][0]**0.5
for i in range(1,N):
	for j in range(0,i):
		hikhjk=0
		for k in range(0,j):
			hikhjk+=chole[i][k]*chole[j][k]
		chole[i][j]=(hilbert[i][j]-hikhjk)/chole[j][j]
		hik=hilbert[i][i]
	for k in range(0,i):
		hik-=(chole[i][k])**2
	chole[i][i]=hik**0.5
writ='The Cholesky decomposition of H_'+str(N)+' is: (printing in a left-triangular manner)'
for i in range (0,N):
	writ+=('\n')	
	writ+=(str(chole[i]))
f.write(writ)
#以下开始解方程
for m in range (2,N+1):
	x = [1 for i in range(1,m+1)]
	for i in range (0,m):
		for j in range (0,i):
			x[i]-=chole[i][j]*x[j]
		x[i]=x[i]/chole[i][i]
	writ='\n\n\nWe have solved H\'y=b first for n='+str(m)+'.\ny='+str(x)+'\nNow solving H_x=y...'
	f.write(writ)
	for i in range (1,m+1):
		for j in range (1,i):
			x[m-i]-=chole[m-j][m-i]*x[m-j]
		x[m-i]=x[m-i]/chole[m-i][m-i]
	f.write('\nx=')
	f.write (str(x))
	f.write ('\n')
	f.write ('is the solution for n=')
	f.write (str(m))
	f.write ('.')
print ('The results have been written to choleskyOutput.txt. The program shall stop now.')
f.close()




