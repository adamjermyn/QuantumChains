import analysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import random
import multiprocessing

f=0.01
numPhi=2000
phiMax=6.28
numE=3
n=7
cpus=8

def spect(phi):
	temp=analysis.eigenPotts(n,f,phiMax*phi/numPhi,numE)
	t0=np.sort(temp[0][0])[:numE]
	t1=np.sort(temp[1][0])[:numE]
	t2=np.sort(temp[2][0])[:numE]
	rms=np.sqrt(np.linalg.norm(t0-t1)**2+np.linalg.norm(t1-t2)**2+np.linalg.norm(t0-t2)**2)
	return [phi*phiMax/numPhi,rms]

pool=multiprocessing.Pool(processes=cpus)
data=pool.map(spect,range(numPhi))

data1=zip(*data)
data1[0]=list(data1[0])
data1[1]=list(data1[1])
plt.scatter(data1[0],data1[1])
plt.show()