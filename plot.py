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
            ax.plot_trisurf(tr,sol,cmap=plt.cm.CMRmap) 
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
        elif choice == 'exact':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,exa)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
        elif choice == 'sx':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,solx)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
        elif choice == 'sy':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,soly)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
        elif choice == 'dx':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,diffx)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
        elif choice == 'dy':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,diffy)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
        elif choice == 'rx':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,ratx)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
        elif choice == 'ry':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,raty)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
        elif choice == 'ex':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,ex)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
        elif choice == 'ey':
            ax = Axes3D(plt.gcf())
            ax.plot_trisurf(tr,ey)
            ax.set_xlabel('X axis')
            ax.set_ylabel('Y axis')
            ax.set_zlabel('Z axis')
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
        fig = plt.figure()
        conv = np.loadtxt('files/convsol.dat')
        c = conv.T[1]
        conv = conv.T[0]
        ax = fig.add_subplot(111)
#        ax.loglog(c,conv,label='Result')
#       ax.loglog(c,c**(slope),label=='Expected')
        
        ax.plot(np.log(c),np.log(conv),label='Result')
        ax.plot(np.log(c),np.log(c**(slope)),label='Expected')
        ax.set_title('Linear elements $\psi$ convergence')
            
        ax.set_xlabel('log(N)')
        ax.set_ylabel('log(Error)')
        ax.legend()
        ax.plot()  
        

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

