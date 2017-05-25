import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
import scipy.sparse as spr

def plottest(choice):
  t=np.loadtxt('./infiles/t.txt')
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
  exbx = np.loadtxt('./files/exbx.dat')
  exby = np.loadtxt('./files/exby.dat')
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
      x = boundx[:,1]
      y = boundx[:,0]
      plt.scatter(x,y,color='red')
      plt.show()
  elif choice == 'boundy':
      x = boundy[:,1]
      y = boundy[:,0]
      plt.scatter(x,y,color='red')
      plt.show()
  elif choice == 'exbx':
      x = exbx[:,1]
      y = exbx[:,0]
      plt.scatter(x,y)
      plt.show()
  elif choice == 'exby':
      x = exby[:,1]
      y = exby[:,0]
      plt.scatter(x,y)
      plt.show()



def conver(slope=-2,der=0):
  fig = plt.figure()
  if der==0:
      conv = np.loadtxt('files/convsol.dat')
  elif der=='x':
      conv = np.loadtxt('files/convx.dat')
  elif der=='y':
      conv = np.loadtxt('files/convy.dat')
  elif der=='xx':
      conv = np.loadtxt('files/convxx.dat')
  elif der=='xy':
      conv = np.loadtxt('files/convxy.dat')
  elif der=='yy':
      conv = np.loadtxt('files/convyy.dat')
  c = conv[:,1]
  conv = conv[:,0]
  ax = fig.add_subplot(111)
#  ax.loglog(c,conv,label='Result')
#  ax.loglog(c,c**(slope),label=='Expected')
  
  ax.plot(np.log10(c),np.log10(conv),label='Result')
  ax.plot(np.log10(c),slope*np.log10(c),label='Expected')
  ax.set_title('Linear elements $u$ convergence')
      
  ax.set_xlabel('log10(N)')
  ax.set_ylabel('log10(Error)')
  ax.legend()
  ax.plot()  
  

  return 0
def conv(slope=-2):
  fig = plt.figure()
  conv = np.loadtxt('files/conv.dat')
  print(conv)
  c = np.array([10.,20.,40.,80.])
  ax = fig.add_subplot(111)
#  ax.loglog(c,conv,label='Result')
#       ax.loglog(c,c**(slope),label=='Expected')
  
  ax.loglog(c,conv)
  ax.loglog(c,c**(slope))
  ax.set_title('Linear elements $\psi$ convergence')
      
  ax.set_xlabel('log10(N)')
  ax.set_ylabel('log10(Error)')
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

