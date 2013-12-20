import analysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import random
import multiprocessing

f=0.5
numPhi=512
phiMax=6.28
numE=54
n=7
cpus=8

def spect(phi):
	temp=analysis.eigenPotts(n,f,phiMax*phi/numPhi,numE)
	return [(phi*phiMax/numPhi,j,k) for j,k in zip(temp[0],temp[2])]

pool=multiprocessing.Pool(processes=cpus)
data=pool.map(spect,range(numPhi))

data1=[]
for i in range(len(data)):
	data1=data1+data[i]
data1=zip(*data1)
data1[0]=list(data1[0])
data1[1]=list(data1[1])
data1[2]=list(data1[2])
data1[2]=np.round(np.angle(data1[2])*3/(2*3.141592),3)
plt.scatter(data1[0],data1[1],c=data1[2],cmap=matplotlib.cm.rainbow,lw = 0,s=12)
print data1[2]
plt.show()