import math
import numpy as np
import scipy.sparse.linalg as lin
from operators import *
#import matplotlib
#import matplotlib.pyplot as plt

def eigenSys(n,h,a,numE):
	hamilPos,hamilNeg=isingHParitySplit(n,h,a)
	eigenPos=lin.eigsh(hamilPos,k=numE/2,which="SA",maxiter=10000)
	eigenNeg=lin.eigsh(hamilNeg,k=numE/2,which="SA",maxiter=10000)
	eigenNet=lin.eigsh(isingHInter(n,h,a),k=numE,which="SA",maxiter=10000)
	return (eigenPos[0],eigenNeg[0],eigenPos[0]-eigenNeg[0],[parityToFull(a,0) for a in np.transpose(eigenPos[1])],[parityToFull(a,1) for a in np.transpose(eigenNeg[1])],eigenNet[0],eigenNet[1])

def eigenSysCut(n,h,a,numE,cutoff):
	eCut=lambda x: eCutoff(x,cutoff)
	hamilPos,hamilNeg=map(eCut,JisingHParitySplit(n,h,a))
	eigenPos=lin.eigsh(hamilPos,k=numE/2,which="SA",maxiter=10000)
	eigenNeg=lin.eigsh(hamilNeg,k=numE/2,which="SA",maxiter=10000)
	eigenNet=lin.eigsh(isingHInter(n,h,a),k=numE,which="SA",maxiter=10000)
	return (eigenPos[0],eigenNeg[0],eigenPos[0]-eigenNeg[0],[parityToFull(a,0) for a in np.transpose(eigenPos[1])],[parityToFull(a,1) for a in np.transpose(eigenNeg[1])],eigenNet[0],eigenNet[1])


def bigEigenSys(n,h,a,numE):
	hamilPos,hamilNeg=isingHParitySplit(n,h,a)
	eigenPos=lin.eigsh(hamilPos,k=numE/2,which="LA",maxiter=10000)
	eigenNeg=lin.eigsh(hamilNeg,k=numE/2,which="LA",maxiter=10000)
	eigenNet=lin.eigsh(isingHInter(n,h,a),k=numE,which="LA",maxiter=10000,tol=1e-20)
	return (eigenPos[0],eigenNeg[0],eigenPos[0]-eigenNeg[0],[parityToFull(a,0) for a in np.transpose(eigenPos[1])],[parityToFull(a,1) for a in np.transpose(eigenNeg[1])],eigenNet[0],eigenNet[1])

def medEigenSys(n,h,a,numE):
	#in specifying sigma, we need to alter the mode
	hamilPos,hamilNeg=isingHParitySplit(n,h,a)
	eigenPos=lin.eigsh(hamilPos,k=numE/2,which="LM",maxiter=1000,sigma=0,mode='normal')
	eigenNeg=lin.eigsh(hamilNeg,k=numE/2,which="LM",maxiter=1000,sigma=0,mode='normal')
	eigenNet=lin.eigsh(isingHInter(n,h,a),k=numE,which="LM",maxiter=1000,sigma=0,mode='normal')
	return (eigenPos[0],eigenNeg[0],eigenPos[0]-eigenNeg[0],[parityToFull(a,0) for a in np.transpose(eigenPos[1])],[parityToFull(a,1) for a in np.transpose(eigenNeg[1])],eigenNet[0],eigenNet[1])

def fullPotts(n,f,phi,numE):
	return lin.eigsh(pottsH(n,f,phi),k=numE,which="SA",maxiter=10000)

def eigenPotts(n,f,phi,numE):
	data0=lin.eigsh(pottsHParitySplit(n,f,phi,0),k=numE,which="SA",maxiter=100000,return_eigenvectors=False)
	data1=lin.eigsh(pottsHParitySplit(n,f,phi,1),k=numE,which="SA",maxiter=100000,return_eigenvectors=False)
	data2=lin.eigsh(pottsHParitySplit(n,f,phi,2),k=numE,which="SA",maxiter=100000,return_eigenvectors=False)
	return (data0,data1,data2)

def triality(states):
	n=int(math.log(len(states[0]))/math.log(3))
	trial=trialityOp(n)
	return [np.conjugate(states[i]).dot(trial.dot(states[i])) for i in range(len(states))]
	

def sign(n):
	bits=[1 if digit=='1' else 0 for digit in bin(n)[2:]]
	net=sum((bits[i]*(1-bits[j]) for i in range(len(bits)) for j in range(i-1)))
	return (-1)**net

def spinToFermion(state):
	for i in range(len(state)):
		state[i]*=sign(i)
	return state

def parityToFull(state,parity): #Maps a spin-parity state onto a full state
	#The binary representation of the index of each basis element is the state, so the parity operation just expands the state vector to twice
	#its dimension and appends an extra binary digit on the end.
	ret=np.zeros(2*len(state))
	for i in range(len(state)):
		if parity==0: #Even
			end=sum([1 if digit=='1' else 0 for digit in bin(i)[2:]])%2
		else:
			end=(1+sum([1 if digit=='1' else 0 for digit in bin(i)[2:]]))%2
		ret[2*i+end]=state[i]
	return ret

def entanglementValues(state,numVals):
	m=np.zeros((math.sqrt(len(state)),math.sqrt(len(state))))
	for i in range(len(state)):
		m[int(i/math.sqrt(len(state)))][int(i%math.sqrt(len(state)))]=state[i]*sign(i)
	_,s,_=lin.svds(m,k=numVals,which="LA")
	return s

def entanglementSys(n,h,a,numS,numVals):
	x=np.zeros(numS*numVals)
	y=np.zeros(numS*numVals)
	sys=eigenSys(n,h,a,numS)
	for i in range(len(np.transpose(sys[4]))):
		a=0-2*np.log(entanglementValues(np.transpose(sys[4])[i],numVals))
		for j in range(len(a)):
			x[numVals*i+j]=i
			y[numVals*i+j]=a[j]
#	plt.scatter(x,y,s=1)
#	plt.xlabel('State Number')
#	plt.ylabel('Entanglement Eigenvalues')
#	plt.show()

def inverseExp(state):
	return np.conjugate(state).dot(inversionOp(int(math.log(len(state),2))).dot(state))

def splitByParity(states,energies):
	par=[(states[i],inverseExp(spinToFermion(states[i])),energies[i]) for i in range(len(energies))]
	even=[a[2] for a in par if a[1]>0]
	odd=[a[2] for a in par if a[1]<0]
	return (even,odd)
