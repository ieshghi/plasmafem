import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation

def plotting(m):
	t = np.loadtxt('t.dat')
	p = np.loadtxt('p.dat')
	fu = np.loadtxt('fu.dat')
	
	tr = Triangulation(p[:,0],p[:,1],triangles = t[:m,:]-1)
	ax = Axes3D(plt.gcf())
	ax.plot_trisurf(tr,x)
	return 0
