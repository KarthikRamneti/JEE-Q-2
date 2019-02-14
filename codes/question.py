import numpy as np
import matplotlib.pyplot as plt
import subprocess
import shlex

def mid_pt(B,C):
	D = (B+C)/2
	return D
	
def norm_vec(AB):
	return np.matmul(omat,np.matmul(AB,dvec))
	
def line_intersect(AD,CF):
	n1 = norm_vec(AD)
	n2 = norm_vec(CF)
	N = np.vstack((n1,n2))
	p = np.zeros(2)
	p[0] = np.matmul(n1,AD[:,0])
	p[1] = np.matmul(n2,CF[:,0])
	return np.matmul(np.linalg.inv(N),p)
	

len = 100
lam_1 = np.linspace(0,1,len)
lam_2 = np.linspace(-6,6,len)
x_G = np.zeros((2,len))
x_C = np.zeros((2,len))
x_AB = np.zeros((2,len))
x_CA = np.zeros((2,len))
x_BC = np.zeros((2,len))

for i in range(len):
	A = np.array([2,5], dtype = float)
	B = np.array([4,-11], dtype = float)
	C_V = np.array([1,(-13)/7], dtype = float)
	C = np.array([lam_2[i],(-9*lam_2[i]-4)/7], dtype = float)
	D = mid_pt(B,C)
	E = mid_pt(A,C)
	F = mid_pt(A,B)
	D1 = mid_pt(B,C_V)
	E1 = mid_pt(A,C_V)


	AD = np.vstack((A,D)).T
	CF = np.vstack((C,F)).T
	AD1 = np.vstack((A,D1)).T
	CF1 = np.vstack((C_V,F)).T

	dvec = np.array([-1,1])
	omat = np.array([[0,1],[-1,0]])

	G = line_intersect(AD,CF)
	G1 = line_intersect(AD1,CF1)

	temp1 = A + lam_1[i]*(B-A)
	x_AB[:,i]= temp1
	temp2 = G
	x_G[:,i]= temp2
	temp3 = A + lam_1[i]*(C_V-A)
	x_CA[:,i]= temp3
	temp4 = C
	x_C[:,i]= temp4
	temp5 = B + lam_1[i]*(C_V-B)
	x_BC[:,i]= temp5

plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$AB$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$AB$')
plt.plot(x_G[0,:],x_G[1,:],label='$Locus$')

plt.plot(x_C[0,:],x_C[1,:],label='$C$')

plt.plot(A[0],A[1], 'o')
plt.text(A[0]*(1+0.1),A[1]*(1-0.1),'A(2,5)')
plt.plot(B[0],B[1], 'o')
plt.text(B[0]*(1-0.2),B[1]*(1),'B(4,-11)')
plt.plot(C_V[0],C_V[1], 'o')
plt.text(C_V[0]*(1-0.3),C_V[1]*(1-0.2),'C1')
plt.plot(G1[0],G1[1], 'o')
plt.text(G1[0]*(1+0.05),G1[1]*(1),'G1')


plt.xlabel('$x$')
plt.xlabel('$y$')
plt.legend(loc='best')
plt.grid()

plt.show()
