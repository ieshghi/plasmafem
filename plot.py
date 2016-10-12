import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
import scipy.sparse as spr

def plottest(choice):
	t = np.loadtxt('./infiles/t.txt')
	p = np.loadtxt('./infiles/p.txt')
        sol = np.loadtxt('./files/sol.dat')
        solx = np.loadtxt('./files/solx.dat')
        soly = np.loadtxt('./files/soly.dat')
        diffx = np.loadtxt('./files/diffx.dat')
        diffy = np.loadtxt('./files/diffy.dat')
        ratx = np.loadtxt('./files/ratx.dat')
        raty = np.loadtxt('./files/raty.dat')
        ex = np.loadtxt('./files/ex.dat')
        ey = np.loadtxt('./files/ey.dat')
        exa = np.loadtxt('./files/exact.dat')
	boundx = np.loadtxt('./files/boundx.dat')
	boundy = np.loadtxt('./files/boundy.dat')
        tr = Triangulation(p[:,0],p[:,1],triangles = t-1)       
        if choice == 's':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,sol) 
        elif choice == 'exact':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,exa)
        elif choice == 'sx':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,solx)
        elif choice == 'sy':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,soly)
        elif choice == 'dx':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,diffx)
        elif choice == 'dy':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,diffy)
        elif choice == 'rx':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,ratx)
        elif choice == 'ry':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,raty)
        elif choice == 'ex':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,ex)
        elif choice == 'ey':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,ey)
        elif choice == 'boundx':
            N = boundx.size
            x = np.linspace(0,1,N)
            plt.plot(x,boundx)
            plt.show()
        elif choice == 'boundy':
            N = boundy.size
            x = np.linspace(0,1,N)
            plt.plot(x,boundy)
            plt.show()



def conver(slope=-2):
        conv = np.loadtxt('files/conv.dat')
        c = conv.T[1]
        conv = conv.T[0]
        plt.loglog(c,conv)
        plt.loglog(c,c**(slope))
        return 0


def plotting(ex = 0):
	t = np.loadtxt('./infiles/t.dat')
	p = np.loadtxt('./infiles/p.dat')
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

