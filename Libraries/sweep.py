import analysis
import operators
import numpy as np
import cPickle as pickle

def numE(n):
	return min(2**n-2,50)


data={}
f=open('data','r+b')

try:
	data=pickle.load(f)
except:
	pass

#do stuff

nn=[i+4 for i in range(10)]
hh=np.linspace(0,1,num=21)
aa=np.linspace(0,1,num=21)

for n in nn:
	for h in hh:
		for a in aa:
			data[(n,h,a)]=analysis.eigenSys(n,h,a,numE(n))
			print (n,h,a)

f=open('data','wb')
pickle.dump(data,f)
