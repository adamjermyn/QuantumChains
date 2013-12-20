import analysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import random
import multiprocessing

f=0.1
phi=0
numE=54
cpus=4

def spect(n):
	temp=analysis.eigenPotts(n,f,phi,min(3*numE,3**(n-2)-1))
	t0=np.sort(temp[0][0])[:numE]
	t1=np.sort(temp[1][0])[:numE]
	t2=np.sort(temp[2][0])[:numE]
	rms=np.sqrt(np.linalg.norm(t0-t1)**2+np.linalg.norm(t1-t2)**2+np.linalg.norm(t0-t2)**2)
	return [n,np.log(rms)]

pool=multiprocessing.Pool(processes=cpus)
data=pool.map(spect,[3+i for i in range(6)])

data1=zip(*data)
data1[0]=list(data1[0])
data1[1]=list(data1[1])
print data1
plt.scatter(data1[0],data1[1])
plt.show()