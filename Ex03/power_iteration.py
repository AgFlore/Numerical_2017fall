def A (a):
	b = [2*a[i]-a[i-1]-a[(i+1)%(len(a))] for i in range(len(a))]
	return b

def P (a,b):
	c = [a[i]*b[i] for i in range(len(a))]
	return (sum(c))

a = [0 for i in range(10)]
a[0] = 1
for i in range(1000):
	a = A(a)
	M = P(a,a)
	a = [a[j]/M**0.5 for j in range(len(a))]

lmax=P(a,A(a))
print(lmax)
print(a)