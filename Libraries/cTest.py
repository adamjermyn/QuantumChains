import analysis
import numpy as np
from math import pi
from functools import partial as part
import functools as func
import matplotlib.pyplot as plt
import matplotlib
import random

numE=3

def spect(n,phi,f):
	try:
		temp=analysis.eigenPotts(n,f,phi,numE)
		t0=np.sort(temp[0][0])[:numE]
		t1=np.sort(temp[1][0])[:numE]
		t2=np.sort(temp[2][0])[:numE]
		rms=np.sqrt(np.linalg.norm(t0-t1)**2+np.linalg.norm(t1-t2)**2+np.linalg.norm(t0-t2)**2)
		return [n,f,phi,rms]
	except:
		return [0,0,0,0]

phis=np.arange(0,2*pi,2*pi/600)
plt.subplot(221)
data=[spect(6,i,0.1) for i in phis]
data1=zip(*data)
plt.scatter(data1[2],data1[3])
plt.title('N=6')
plt.subplot(222)
data=[spect(7,i,0.1) for i in phis]
data1=zip(*data)
plt.scatter(data1[2],data1[3])
plt.title('N=7')
plt.subplot(223)
data=[spect(8,i,0.1) for i in phis]
data1=zip(*data)
plt.scatter(data1[2],data1[3])
plt.title('N=8')
plt.subplot(224)
data=[spect(9,i,0.1) for i in phis]
data1=zip(*data)
plt.scatter(data1[2],data1[3])
plt.title('N=9')
plt.show()