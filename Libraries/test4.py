import analysis
import numpy as np
import matplotlib.pyplot as plt
import random

numF=10
maxF=0.1
numE=54
n=8
data=[[(f*maxF/numF,j) for j in analysis.eigenPotts(n,maxF*f/numF,0,numE)[0]] for f in range(numF)]
data1=[]
for i in range(len(data)):
	data1=data1+data[i]
data1=zip(*data1)
data1[0]=list(data1[0])
data1[1]=list(data1[1])
plt.scatter(data1[0],data1[1],s=1)
plt.show()