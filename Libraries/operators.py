import math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as lin
from scipy.sparse import bsr_matrix
from scipy.sparse import dia_matrix

sx = np.array([[0, 1],[ 1, 0]])
sy = np.array([[0, -1j],[1j, 0]])
sz = np.array([[1, 0],[0, -1]])
tau=np.array([[0,0,1],[1,0,0],[0,1,0]])
tauD=np.ma.conjugate(np.transpose(tau))
sigma=np.array([[1,0,0],[0,-0.5+1j*math.sqrt(3)/2,0],[0,0,-0.5-1j*math.sqrt(3)/2]])
sigmaD=np.ma.conjugate(np.transpose(sigma))

def kronProd(data,n,d): #data should be a dict of (index,array). n gives the overall size of the space. d gives the size of the arrays.
	temp=[data.get(i,np.identity(d)) for i in range(n)]
	ret=bsr_matrix(np.identity(1))
	for i in range(n):
		ret=sp.kron(ret,temp[i])
	return ret

def quadTerm(a,b,ai,bi,n): #a,b are arrays, ai,bi are their indices in a list of identities, and n is the size of the space.
	hamil={}
	hamil.update({ai:a,bi:b})
	return kronProd(hamil,n,len(a))

def chi(n,i):
	ret=bsr_matrix(np.identity(1))
	for k in range(i-1):
		ret=sp.kron(ret,sx)
	ret=sp.kron(ret,sy)
	for k in range(n-i):
		if ret.shape[0]<math.pow(2,n):
			ret=sp.kron(ret,np.identity(2))
	return ret

def psi(n,i):
	ret=bsr_matrix(np.identity(1))
	for k in range(i-1):
		ret=sp.kron(ret,sx)
	ret=sp.kron(ret,sz)
	for k in range(n-i):
		if ret.shape[0]<math.pow(2,n):
			ret=sp.kron(ret,np.identity(2))
	return ret

def isingOp(n,i):
	if i%2==0:
		return chi(n,i/2)
	else:
		return psi(n,(i+1)/2)

def inversionOp(n):
	m=np.zeros((2**n,2**n))
	for i in xrange(2**n):
		m[i][sum(1<<(n-1-j) for j in range(n) if i>>j&1)]=1
	return m

def isingHInter(n,h,a): #use a=0 for the noninteracting case
	return -sum((quadTerm(sz,sz,i,i+1,n)+a*quadTerm(sx,sx,i,i+1,n) for i in range(n-1)))-h*sum((kronProd({i:sx},n,2) for i in range(n)))

def trialityOp(n):
	return kronProd(dict([(i,tauD) for i in range(n)]),n,3)

def isingHParitySplit(n,h,a):
	hamil=-sum((quadTerm(sx,sx,i,i+1,n-1)+a*quadTerm(sz,sz,i,i+1,n-1) for i in range(n-2)))-h*sum((kronProd({i:sz},n-1,2) for i in range(n-1)))
	hTerm=h*kronProd(dict((i,sz) for i in range(n-1)),n-1,2)
	jTerm=kronProd({n-2:sx},n-1,2)
	aTerm=a*kronProd(dict((i,sz) for i in range(n-2)),n-1,2)
	hamilPos=hamil-hTerm-jTerm-aTerm
	hamilNeg=hamil+hTerm-jTerm+aTerm
	return (hamilPos,hamilNeg)

def JisingHParitySplit(n,h,a):
	hamil=-sum((quadTerm(sz,sz,i,i+1,n-1)+a*quadTerm(sx,sx,i,i+1,n-1) for i in range(n-2)))-h*sum((kronProd({i:sx},n-1,2) for i in range(n-1)))
	hTerm=h*kronProd(dict((i,sx) for i in range(n-1)),n-1,2)
	jTerm=kronProd({n-2:sz},n-1,2)
	aTerm=a*kronProd(dict((i,sx) for i in range(n-2)),n-1,2)
	hamilPos=hamil-hTerm-jTerm-aTerm
	hamilNeg=hamil+hTerm-jTerm+aTerm
	return (hamilPos,hamilNeg)

def eCutoff(op,cutoff):
	minn=min(op.diagonal())
	diag=1*(op.diagonal()<(minn+cutoff))
	proj=dia_matrix(([diag],[0]),shape=(len(diag),len(diag)))
	proj=bsr_matrix(proj)
	return bsr_matrix(proj.dot(op.dot(proj)))


def pottsH(n,f,phi):
	phase=np.exp(-1j*phi)
	#t1=-f*sum((kronProd({i:tauD*phase+tau/phase},n,3) for i in range(n))) # This line allows there to be a phase on the f term as well
	t1=-f*sum((kronProd({i:tauD+tau},n,3) for i in range(n)))
	t2=-sum((quadTerm(sigmaD,sigma*phase,i,i+1,n) for i in range(n-1)))
	t3=-sum((quadTerm(sigma,sigmaD/phase,i,i+1,n) for i in range(n-1)))
	return t1+t2+t3

def pottsHParitySplit(n,f,phi,j):
	phase=np.exp(-1j*phi)
	#hamil=-f*sum((kronProd({i:sigmaD*phase+sigma/phase},n-1,3) for i in range(n-1))) This line allows there to be a phase on the f term as well
	hamil=-f*sum((kronProd({i:sigmaD+sigma},n-1,3) for i in range(n-1)))
	hamil=hamil-sum((quadTerm(tauD,tau*phase,i,i+1,n-1) for i in range(n-2)))
	hamil=hamil-sum((quadTerm(tau,tauD/phase,i,i+1,n-1) for i in range(n-2)))
	hamil=hamil-kronProd({n-2:tauD},n-1,3)*phase-kronProd({n-2:tau},n-1,3)/phase
	fact=1
	if j==1:
		fact=np.exp(1j*2*math.pi/3)
	elif j==2:
		fact=np.exp(1j*4*math.pi/3)
	#hamil=hamil-fact*f*(phase*kronProd(dict((i,sigma) for i in range(n-1)),n-1,3))
	#hamil=hamil-f*(kronProd(dict((i,sigmaD) for i in range(n-1)),n-1,3)/phase)/fact
	hamil=hamil-fact*f*(kronProd(dict((i,sigma) for i in range(n-1)),n-1,3))
	hamil=hamil-f*(kronProd(dict((i,sigmaD) for i in range(n-1)),n-1,3))/fact
	return hamil

def freeFermion(n,h): #This is done in the Majorana basis.
	hamil=np.zeros(2*n,2*n)
	for i in range(2*n-1):
		if i%2==0:
			hamil[i+1][i]=h
			hamil[i][i+1]=-h
		else:
			hamil[i+1][i]=1
			hamil[i][i+1]=-1
	hamil=1j*hamil
	return hamil

def freeFermionOperators(n,h):
	hamil=freeFermion(n,h)
	return np.linalg.eig(hamil)[1]

def majoranaFromFermion(n,h):
	vects=freeFermionOperators(n,h)
	return np.linalg.inv(np.transpose(vects))






























