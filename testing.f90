program testing
  use mesh
  use curvestuff
  implicit none
  real *8:: ep,del,kap,infi,findif,d1,d2,d3,stdev
  real *8, dimension(:),allocatable::solx,soly,sol,ubx,uby,ex,ey,exa,ratx,raty,diffx,diffy,exbx,exby
  real *8, dimension(:,:),allocatable::p
  integer,dimension(:),allocatable::b
  integer,dimension(:,:),allocatable::t
  integer::i,n,nb

  ep = 0.32d0
  del = 0.33d0
  kap = 1.7d0
  infi = 1.0d0*1e-14
  findif = 1.0d0*1e-3
  call switchpars(ep,del,kap,d1,d2,d3)

  call derpois(d1,d2,d3,infi,findif,solx,soly,sol,p,t,b,ubx,uby)

  n = size(sol)
  nb = size(b)

  allocate(ex(n),ey(n),exa(n),diffx(n),ratx(n),diffy(n),raty(n),exbx(nb),exby(nb)) !store exact solutions, for comparison

  do i=1,n
    exa(i) = exact(p(i,1),p(i,2),d1,d2,d3)
    ex(i) = exactx(p(i,1),p(i,2),d1,d2,d3)
    ey(i) = exacty(p(i,1),p(i,2),d1,d2,d3)
    diffx(i) = ex(i)-solx(i)
    diffy(i) = ey(i)-soly(i)
    ratx(i) = solx(i)/ex(i)
    raty(i) = soly(i)/ey(i)
  enddo

  do i=1,nb
    exbx(i) = exactx(p(b(i),1),p(b(i),2),d1,d2,d3)
    exby(i) = exacty(p(b(i),1),p(b(i),2),d1,d2,d3)
  enddo

  stdev = 0.0d0
  open(1,file='files/convsol.dat',position='append')
    do i=1,n
      stdev = stdev + ((exa(i)-sol(i))**2)/n
    enddo
    stdev = sqrt(stdev)
    write(1,*)  stdev,int(sqrt(float(n)))
  close(1) 

  stdev = 0.0d0
  open(1,file='files/convx.dat',position='append')
    do i=1,n
      stdev = stdev + ((ex(i)-solx(i))**2)/n
    enddo
    write(1,*)  stdev,int(sqrt(float(n)))
  close(1) 

  stdev = 0.0d0
  open(1,file='files/convy.dat',position='append')
    do i=1,n
      stdev = stdev + ((ey(i)-soly(i))**2)/n
    enddo
    write(1,*)  stdev,int(sqrt(float(n)))
  close(1) 

  open(1,file='files/sol.dat')
    do i=1,size(sol)
      write(1,*) sol(i)
    enddo
  close(1) 
  
  open(1,file='files/solx.dat')
    do i=1,size(solx)
      write(1,*) solx(i)
    enddo
  close(1) 
  
  open(1,file='files/soly.dat')
    do i=1,size(soly)
      write(1,*) soly(i)
    enddo
  close(1) 

  open(1,file='files/ex.dat')
    do i=1,size(ex)
      write(1,*) ex(i)
    enddo
  close(1) 

  open(1,file='files/ey.dat')
    do i=1,size(ey)
      write(1,*) ey(i)
    enddo
  close(1) 

  open(1,file='files/exact.dat')
    do i=1,size(exa)
      write(1,*) exa(i)
    enddo
  close(1) 

  open(1,file='files/diffx.dat')
    do i=1,size(diffx)
      write(1,*) diffx(i)
    enddo
  close(1) 
  
  open(1,file='files/diffy.dat')
    do i=1,size(diffy)
      write(1,*) diffy(i)
    enddo
  close(1) 

  open(1,file='files/ratx.dat')
    do i=1,size(ratx)
      write(1,*) ratx(i)
    enddo
  close(1) 

  open(1,file='files/raty.dat')
    do i=1,size(raty)
      write(1,*) raty(i)
    enddo
  close(1) 

  open(1,file='files/boundx.dat')
    do i=1,size(ubx)
      write(1,*) ubx(i)
    enddo
  close(1) 

  open(1,file='files/boundy.dat')
    do i=1,size(uby)
      write(1,*) uby(i)
    enddo
  close(1) 

  open(1,file='files/exbx.dat')
    do i=1,size(exbx)
      write(1,*) exby(i)
    enddo
  close(1) 
  
  open(1,file='files/exby.dat')
    do i=1,size(exbx)
      write(1,*) exby(i)
    enddo
  close(1) 



endprogram testing
