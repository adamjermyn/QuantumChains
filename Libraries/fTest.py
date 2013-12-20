import analysis
import numpy as np
from math import pi
from functools import partial as part
import functools as func
import matplotlib.pyplot as plt
import matplotlib
import random

numE=3

def spect0(n,phi,f):
	try:
		temp=analysis.eigenPotts(n,f,phi,numE)
		t0=np.sort(temp[0][0])[:numE]
		t1=np.sort(temp[1][0])[:numE]
		t2=np.sort(temp[2][0])[:numE]
		return [n,f,phi,t0[0]]
	except:
		return [0,0,0,0]
def spect1(n,phi,f):
	try:
		temp=analysis.eigenPotts(n,f,phi,numE)
		t0=np.sort(temp[0][0])[:numE]
		t1=np.sort(temp[1][0])[:numE]
		t2=np.sort(temp[2][0])[:numE]
		return [n,f,phi,t1[0]]
	except:
		return [0,0,0,0]
def spect2(n,phi,f):
	try:
		temp=analysis.eigenPotts(n,f,phi,numE)
		t0=np.sort(temp[0][0])[:numE]
		t1=np.sort(temp[1][0])[:numE]
		t2=np.sort(temp[2][0])[:numE]
		return [n,f,phi,t2[0]]
	except:
		return [0,0,0,0]


phis=np.arange(0,2*pi,0.01)
plt.subplot(221)
data=[spect0(7,i,0.1) for i in phis]
data1=zip(*data)
plt.scatter(data1[2],data1[3])
plt.subplot(222)
data=[spect1(7,i,0.1) for i in phis]
data1=zip(*data)
plt.scatter(data1[2],data1[3])
plt.subplot(223)
data=[spect2(7,i,0.1) for i in phis]
data1=zip(*data)
plt.scatter(data1[2],data1[3])
plt.show()