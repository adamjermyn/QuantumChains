from operators import *
np.set_printoptions(threshold=np.nan,linewidth=200)
d=np.real(pottsH(3,0.1,0.4).toarray())
print d+np.identity(len(d))*4*np.cos(0.4)