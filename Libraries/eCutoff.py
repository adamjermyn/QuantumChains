import matplotlib.pyplot as plt
from analysis import *

def spectrum(nMin,nMax,h,cut):
	return np.array([(n,np.log(np.abs(eigenSysCut(n+nMin,h,0,4,2)[2][1]))) for n in range(nMax-nMin)])

ref=np.log(np.abs(eigenSys(10,0.4,0,4)[2][1]))
data1=np.log(np.abs([np.sort(eigenSysCut(10,0.4,0,4,cut)[2])[1] for cut in range(35)]))
print data1
data1=data1-ref
plt.scatter(range(35),data1)
plt.show()