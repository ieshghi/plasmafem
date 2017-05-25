program testing
  use mesh
  use curvestuff
  implicit none
  real *8:: ep,del,kap,infi,findif,d1,d2,d3,stdev
  real *8, dimension(:),allocatable::solx,soly,sol,ubx,uby,ex,ey,exa,exx,eyy,exy,ratx,&
          raty,diffx,diffy,exbx,exby,areas,solxx,solxy,solyy
  real *8, dimension(:,:),allocatable::p
  integer,dimension(:),allocatable::b
  integer,dimension(:,:),allocatable::t
  integer::i,n,nb,nt

  ep = 0.32d0
  del = 0.33d0
  kap = 1.7d0
  infi = 1.0d0*1e-14
  findif = 1.0d0*1e-3
  call switchpars(ep,del,kap,d1,d2,d3)
  call dderpois(d1,d2,d3,infi,findif,solx,soly,solxx,solxy,solyy,sol,p,t,areas)
  n = size(sol)
  nb = size(b)
  nt = size(t(:,1))

  allocate(ex(n),ey(n),exa(n),exx(n),eyy(n),exy(n),diffx(n),ratx(n),diffy(n),raty(n),exbx(nb),exby(nb)) !store exact solutions, for comparison

  do i=1,n
    exa(i) = exact(p(i,1),p(i,2),d1,d2,d3)
    ex(i) = exactx(p(i,1),p(i,2),d1,d2,d3)
    ey(i) = exacty(p(i,1),p(i,2),d1,d2,d3)
    exx(i) = exactxx(p(i,1),p(i,2),d1,d2,d3)
    exy(i) = exactxy(p(i,1),p(i,2),d1,d2,d3)
    eyy(i) = exactyy(p(i,1),p(i,2),d1,d2,d3)
    diffx(i) = ex(i)-solx(i)
    diffy(i) = ey(i)-soly(i)
    ratx(i) = solx(i)/ex(i)
    raty(i) = soly(i)/ey(i)
  enddo

  stdev = 0.0d0
  open(1,file='files/convxx.dat',position='append')
    do i=1,n
      stdev = stdev + (areas(i)*(exx(i)-solxx(i))**2)/n
    enddo
    stdev = sqrt(stdev)
    write(1,*)  stdev,int(sqrt(float(nt)))
  close(1) 

  stdev = 0.0d0
  open(1,file='files/convyy.dat',position='append')
    do i=1,n
      stdev = stdev + (areas(i)*(exy(i)-solxy(i))**2)/n
    enddo
    stdev = sqrt(stdev)
    write(1,*)  stdev,int(sqrt(float(nt)))
  close(1) 

  stdev = 0.0d0
  open(1,file='files/convxy.dat',position='append')
    do i=1,n
      stdev = stdev + (areas(i)*(eyy(i)-solyy(i))**2)/n
    enddo
    stdev = sqrt(stdev)
    write(1,*)  stdev,int(sqrt(float(nt)))
  close(1) 

  stdev = 0.0d0
  open(1,file='files/convsol.dat',position='append')
    do i=1,n
      stdev = stdev + (areas(i)*(exa(i)-sol(i))**2)/n
    enddo
    stdev = sqrt(stdev)
    write(1,*)  stdev,int(sqrt(float(nt)))
  close(1) 

  stdev = 0.0d0
  open(1,file='files/convx.dat',position='append')
    do i=1,n
      stdev = stdev + (areas(i)*(ex(i)-solx(i))**2)/n
    enddo
    stdev = sqrt(stdev)
    write(1,*)  stdev,int(sqrt(float(nt)))
  close(1) 

  stdev = 0.0d0
  open(1,file='files/convy.dat',position='append')
    do i=1,n
      stdev = stdev + (areas(i)*(ey(i)-soly(i))**2)/n
    enddo
    stdev = sqrt(stdev)
    write(1,*)  stdev,int(sqrt(float(nt)))
  close(1) 


endprogram testing

