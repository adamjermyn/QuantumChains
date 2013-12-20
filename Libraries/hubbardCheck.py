from analysis import *

def measure(data):
	return np.log(np.sqrt((data[0]-data[1])**2+(data[0]-data[2])**2+(data[1]-data[2])**2))

data=[measure(eigenPotts(3+n,0.005,0,3)) for n in range(6)]
data=np.array([i-data[0] for i in data])
print data