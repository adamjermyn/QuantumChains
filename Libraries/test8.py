import analysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import random
import multiprocessing

phi=0
numF=10
maxF=0.2
numE=3
cpus=8

def spect6(f):
	temp=analysis.eigenPotts(6,maxF*f/numF,phi,3*numE)
	t0=np.sort(temp[0][0])[:numE]
	t1=np.sort(temp[1][0])[:numE]
	t2=np.sort(temp[2][0])[:numE]
	rms=np.sqrt(np.linalg.norm(t0-t1)**2+np.linalg.norm(t1-t2)**2+np.linalg.norm(t0-t2)**2)
	return [maxF*f/numF,rms]
def spect7(f):
	temp=analysis.eigenPotts(7,maxF*f/numF,phi,3*numE)
	t0=np.sort(temp[0][0])[:numE]
	t1=np.sort(temp[1][0])[:numE]
	t2=np.sort(temp[2][0])[:numE]
	rms=np.sqrt(np.linalg.norm(t0-t1)**2+np.linalg.norm(t1-t2)**2+np.linalg.norm(t0-t2)**2)
	return [maxF*f/numF,rms]
def spect8(f):
	temp=analysis.eigenPotts(8,maxF*f/numF,phi,3*numE)
	t0=np.sort(temp[0][0])[:numE]
	t1=np.sort(temp[1][0])[:numE]
	t2=np.sort(temp[2][0])[:numE]
	rms=np.sqrt(np.linalg.norm(t0-t1)**2+np.linalg.norm(t1-t2)**2+np.linalg.norm(t0-t2)**2)
	return [maxF*f/numF,rms]
def spect9(f):
	temp=analysis.eigenPotts(9,maxF*f/numF,phi,3*numE)
	t0=np.sort(temp[0][0])[:numE]
	t1=np.sort(temp[1][0])[:numE]
	t2=np.sort(temp[2][0])[:numE]
	rms=np.sqrt(np.linalg.norm(t0-t1)**2+np.linalg.norm(t1-t2)**2+np.linalg.norm(t0-t2)**2)
	return [maxF*f/numF,rms]
def spect10(f):
	temp=analysis.eigenPotts(10,maxF*f/numF,phi,3*numE)
	t0=np.sort(temp[0][0])[:numE]
	t1=np.sort(temp[1][0])[:numE]
	t2=np.sort(temp[2][0])[:numE]
	rms=np.sqrt(np.linalg.norm(t0-t1)**2+np.linalg.norm(t1-t2)**2+np.linalg.norm(t0-t2)**2)
	return [maxF*f/numF,rms]

pool=multiprocessing.Pool(processes=cpus)
data0=pool.map(spect6,range(numF))
data1=pool.map(spect7,range(numF))
data2=pool.map(spect8,range(numF))
data3=pool.map(spect9,range(numF))
data4=pool.map(spect10,range(numF))
data0=zip(*data0)
data1=zip(*data1)
data2=zip(*data2)
data3=zip(*data3)
data4=zip(*data4)
plt.subplot(511)
plt.scatter(data0[0],data0[1])
plt.subplot(512)
plt.scatter(data1[0],data1[1])
plt.subplot(513)
plt.scatter(data2[0],data2[1])
plt.subplot(514)
plt.scatter(data3[0],data3[1])
plt.subplot(515)
plt.scatter(data4[0],data4[1])
plt.show()