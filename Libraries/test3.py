import analysis
import numpy as np
import matplotlib.pyplot as plt

def split(n):
	numI=10
	data=[]
	for i in range(numI):
		print i
		try:
			temp=[(0.01*i,np.log(abs(j))) for j in analysis.medEigenSys(n,0.1,0.01*i,2)[2]]
			data.append(temp)
		except:
			print "Convergence problem."
	data1=[]
	for i in range(len(data)):
		data1=data1+data[i]
	data1=zip(*data1)
	return data1

data=[]
for n in range(7):
	print n
	temp=split(n+4)
	if len(temp)>0:
		data.append(temp)
	else:
		data.append([(0),(0)])

print data

plt.subplot(711)
plt.title('N=5')
plt.scatter(data[0][0],data[0][1],s=1)
plt.xlim([0,0.1])
plt.ylim([-26,0])
plt.subplot(712)
plt.title('N=6')
plt.scatter(data[1][0],data[1][1],s=1)
plt.xlim([0,0.1])
plt.ylim([-26,0])
plt.subplot(713)
plt.title('N=7')
plt.scatter(data[2][0],data[2][1],s=1)
plt.xlim([0,0.1])
plt.ylim([-26,0])
plt.subplot(714)
plt.title('N=8')
plt.scatter(data[3][0],data[3][1],s=1)
plt.xlim([0,0.1])
plt.ylim([-26,0])
plt.subplot(715)
plt.title('N=9')
plt.scatter(data[4][0],data[4][1],s=1)
plt.xlim([0,0.1])
plt.ylim([-26,0])
plt.subplot(716)
plt.title('N=10')
plt.scatter(data[5][0],data[5][1],s=1)
plt.xlim([0,0.1])
plt.ylim([-26,0])
plt.subplot(717)
plt.title('N=11')
plt.scatter(data[6][0],data[6][1],s=1)
plt.xlim([0,0.1])
plt.ylim([-26,0])
plt.show()