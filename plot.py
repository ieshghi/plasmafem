import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
import scipy.sparse as spr

def boundmesh():
    a = np.loadtxt('files/upx.txt')
    errx = a[:,0]
    erry = a[:,1]
    theta = a[:,2]*np.pi
    r = a[:,3]
    p = np.loadtxt('infiles/p.txt')
    t = np.loadtxt('infiles/t.txt')-1
    b = np.loadtxt('infiles/b.txt')-1

    xloc = r*np.cos(theta)+1
    yloc = r*np.sin(theta)
    
    bigloc = [errx>5e-3]
    plt.scatter(xloc[bigloc],yloc[bigloc],c='r',marker='x')

    q = Triangulation(p[:,0],p[:,1],triangles=t)
    plt.triplot(q)

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

def trackconv():
  fig = plt.figure()
  conv = np.loadtxt('track.dat')
  fem = abs(conv[:,0])
  error = abs(conv[:,1])
  edge = conv[:,2]
  ax = fig.add_subplot(111)
  
  ax.plot((-1)*np.log10(edge),np.log10(error),label='Minimum track')
  ax.plot((-1)*np.log10(edge),np.log10(fem),label='Finite elements')
  ax.plot((-1)*np.log10(edge),(-2)*np.log10(edge[0])-(-2)*np.log10(edge)+np.log10(fem[0]),label='Second order')

  ax.set_title('Newton convergence')
      
  ax.set_xlabel('-log10(h)')
  ax.set_ylabel('log10(Error)')
  ax.legend()
  ax.plot()

def conver(slope=-2,der=0,offset=0):
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
  elif der == 'bound_y':
      conv = np.loadtxt('files/convuy.dat')
  elif der == 'bound_x':
      conv = np.loadtxt('files/convux.dat')
  c = conv[:,1]
  conv = conv[:,0]
  ax = fig.add_subplot(111)
#  ax.loglog(c,conv,label='Result')
#  ax.loglog(c,c**(slope),label=='Expected')
  
  ax.plot((-1)*np.log10(c),np.log10(conv),'rx',label='Result')
  ax.plot((-1)*np.log10(c),-slope*np.log10(c),label='Expected')

  ax.set_title('Linear elements convergence')
      
  ax.set_xlabel('-log10(h)')
  ax.set_ylabel('log10(Error)')
  ax.legend()
  ax.plot()  
  

  return 0


