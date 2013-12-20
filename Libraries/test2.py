import analysis
import numpy as np
import matplotlib.pyplot as plt

def split():
	numN=9
	data=[[(4+n,np.log(abs(j))) for j in analysis.medEigenSys(4+n,0.1,0.1,2)[2]] for n in range(numN)]
	data1=[]
	for i in range(numN):
		data1=data1+data[i]
	data1=zip(*data1)
	return data1


data=split()

plt.scatter(data[0],data[1])
plt.show()