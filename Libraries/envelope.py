import analysis
import numpy as np
from math import pi
from math import log
from math import exp
from functools import partial as part
import functools as func
import matplotlib.pyplot as plt
import matplotlib
import random
import time

numE=24

def spect(n,phi,f,norm):
	try:
		start=time.time()
		temp=analysis.eigenPotts(n,f,phi,numE)
		t0=np.sort(temp[0])[:numE]
		t1=np.sort(temp[1])[:numE]
		t2=np.sort(temp[2])[:numE]
		rms=np.sqrt(np.linalg.norm(t0-t1)**2+np.linalg.norm(t1-t2)**2+np.linalg.norm(t0-t2)**2)
		print n,log(abs(rms/norm)),time.time()-start
		return (n,phi,log(abs(rms/norm)))
	except:
		print 'err'
		return (0,0,0)

phis=np.arange(0,pi/6,2*pi/300)
norm=[spect(4,i,0.2,1) for i in phis]
data=[(4,i,0) for i in phis]
data=data+[spect(5,phis[i],0.2,exp(norm[i][2])) for i in range(len(phis))]
print 'done'
data=data+[spect(6,phis[i],0.2,exp(norm[i][2])) for i in range(len(phis))]
print 'done'
data=data+[spect(7,phis[i],0.2,exp(norm[i][2])) for i in range(len(phis))]
print 'done'
data=data+[spect(8,phis[i],0.2,exp(norm[i][2])) for i in range(len(phis))]
print 'done'
data=data+[spect(9,phis[i],0.2,exp(norm[i][2])) for i in range(len(phis))]
print 'done'
data=data+[spect(10,phis[i],0.2,exp(norm[i][2])) for i in range(len(phis))]
print 'done'
data=data+[spect(11,phis[i],0.2,exp(norm[i][2])) for i in range(len(phis))]
print 'done'
data=data+[spect(12,phis[i],0.2,exp(norm[i][2])) for i in range(len(phis))]
print 'done'
for i in range(len(phis)):
	y=[data[i+j*len(phis)][2] for j in range(8)]
	y=[y[j+1]-y[j] for j in range(len(y)-1)]
	print phis[i],y


data=zip(*data)
plt.scatter(data[0],data[2],c=data[1])
plt.colorbar()
plt.show()