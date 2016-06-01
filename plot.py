import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
import scipy.sparse as spr

def plotting():
	t = np.loadtxt('./files/t.dat')
	p = np.loadtxt('./files/p.dat')
	fu = np.loadtxt('./files/fu.dat')
	ia = np.loadtxt('./files/ia.dat') -1
	ja = np.loadtxt('./files/ja.dat') -1
	arr = np.loadtxt('./files/arr.dat')
	M = spr.csr_matrix((arr,(ia,ja)))
	M = M.toarray()
	
	
	x = np.linalg.solve(M,fu)
	
	tr = Triangulation(p[:,0],p[:,1],triangles = t-1)
	ax = Axes3D(plt.gcf())
	ax.plot_trisurf(tr,x)
	return 0
