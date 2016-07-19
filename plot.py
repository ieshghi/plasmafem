import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
import scipy.sparse as spr

def plotspec(ret=False):
        a = np.loadtxt('files/deriv.dat')
        b = a[:,2] + 1j*a[:,3]
        f = a[:,0]
        d = a[:,1]
        x = a[:,4]
        plt.plot(x,d)
        plt.plot(x,np.real(b))
        if ret==True :
            return b,f,d
        else:
            return 0


def conver(slope=-2):
        conv = np.loadtxt('files/conv.dat')
        c = conv.T[1]
        conv = conv.T[0]
        plt.loglog(c,conv)
        plt.loglog(c,c**(slope))
        return 0


def plotting(ex = 0):
	t = np.loadtxt('./files/t.dat')
	p = np.loadtxt('./files/p.dat')
	forsol = np.loadtxt('./files/x.dat')
	exact = np.loadtxt('./files/exact.dat')
	tr = Triangulation(p[:,0],p[:,1],triangles = t-1)
	
	if ex == 0:
		ax = Axes3D(plt.gcf())
		ax.plot_trisurf(tr,forsol)
	else:
		ax = Axes3D(plt.gcf())
		ax.plot_trisurf(tr,exact)
                
	return 0

def plotsolve():
	t = np.loadtxt('./files/t.dat')
	p = np.loadtxt('./files/p.dat')
	fu = np.loadtxt('./files/fu.dat')
	ia = np.loadtxt('./files/ia.dat') -1
	ja = np.loadtxt('./files/ja.dat') -1
	arr = np.loadtxt('./files/arr.dat')
	M = spr.csr_matrix((arr,(ia,ja)))
	M = M.toarray()
	pysol = np.linalg.solve(M,fu)
	tr = Triangulation(p[:,0],p[:,1],triangles = t-1)
	ax = Axes3D(plt.gcf())
	ax.plot_trisurf(tr,pysol)
	return 0
	
def comparepy(rat=0,dif=0):
	t = np.loadtxt('./files/t.dat')
	p = np.loadtxt('./files/p.dat')
	fu = np.loadtxt('./files/fu.dat')
	ia = np.loadtxt('./files/ia.dat') -1
	ja = np.loadtxt('./files/ja.dat') -1
	arr = np.loadtxt('./files/arr.dat')
	forsol = np.loadtxt('./files/x.dat')
	M = spr.csr_matrix((arr,(ia,ja)))
	M = M.toarray()
	pysol = np.linalg.solve(M,fu)
	tr = Triangulation(p[:,0],p[:,1],triangles = t-1)
	diff = abs(forsol-pysol)
	ratio = abs(forsol/pysol)

	if rat == 1:
			ax1 = Axes3D(plt.gcf())
			ax1.plot_trisurf(tr,ratio)
	if dif == 1:
			ax2 = Axes3D(plt.gcf())
			ax2.plot_trisurf(tr,diff)

def errorcheck(rat=0,dif=0):
	t = np.loadtxt('./files/t.dat')
	p = np.loadtxt('./files/p.dat')
	forsol = np.loadtxt('./files/x.dat')
	exact = np.loadtxt('./files/exact.dat')
	tr = Triangulation(p[:,0],p[:,1],triangles = t-1)
	diff = forsol-exact
	ratio = forsol/exact
	ratio[ratio==0] = np.nan
	if rat==1:
		ax = Axes3D(plt.gcf())
		ax.plot_trisurf(tr,ratio)
		
	if dif==1:
		ax = Axes3D(plt.gcf())
		ax.plot_trisurf(tr,forsol-exact)
	
	return diff,ratio
