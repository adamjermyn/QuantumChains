import analysis
import numpy as np
import matplotlib.pyplot as plt

def split(n):
	numI=100
	data=[[(0.001*i,np.log(abs(j))) for j in analysis.medEigenSys(n,0.1,0.001*i,10)[2]] for i in range(numI)]
	data1=[]
	for i in range(numI):
		data1=data1+data[i]
	data1=zip(*data1)
	return data1

data=[split(n+7) for n in range(4)]

plt.subplot(411)
plt.scatter(data[0][0],data[0][1],s=1)
plt.subplot(412)
plt.scatter(data[1][0],data[1][1],s=1)
plt.subplot(413)
plt.scatter(data[2][0],data[2][1],s=1)
plt.subplot(414)
plt.scatter(data[3][0],data[3][1],s=1)
plt.show()