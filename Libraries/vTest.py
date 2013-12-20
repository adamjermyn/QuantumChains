import analysis
import numpy as np
from math import pi
from functools import partial as part
import functools as func
import matplotlib.pyplot as plt
import matplotlib
import random

numE=3

def spectDiff(n,phi,f):
	temp=analysis.eigenPotts(n,f,phi,numE)
	t0=np.sort(temp[0][0])[:numE/3]
	t1=np.sort(temp[1][0])[:numE/3]
	t2=np.sort(temp[2][0])[:numE/3]
	temp=np.concatenate((t0,t1,t2))
	temp=np.sort(temp)
	data=np.sort(analysis.fullPotts(n,f,phi,numE)[0])
	return [i-j for i,j in zip(data,temp)]

phis=np.arange(0,2*pi,0.1)
print [spectDiff(5,i,0.1) for i in phis]
