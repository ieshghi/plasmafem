program testing
  use mesh
  use curvestuff
  implicit none
  real *8:: infi,findif,c,d1,d2,d3,d4,gam,stdev,edge,x,y
  real *8, dimension(:),allocatable::solx,soly,sol,ex,ey,exa,exx,eyy,exy,&
          areas,solxx,solxy,solyy,psi,psix,psiy,psixx,psixy,psiyy
  real *8, dimension(:,:),allocatable::p,p2
  integer,dimension(:),allocatable::b,b2
  integer,dimension(:,:),allocatable::t,t2
  integer::i,n,nb,nt
    
  findif = 1e-3
  infi = 1e-14

  call distmesh(p2,t2,b2,d1,d2,d3,d4,c,gam)
  call dderpois(infi,findif,solx,soly,solxx,solxy,solyy,sol,p,t,areas)
  n = size(sol)
  nb = size(b)
  nt = size(t(:,1))

  allocate(ex(n),ey(n),exa(n),exx(n),eyy(n),exy(n),psi(n),psix(n),psixx(n),psiy(n),psiyy(n),psixy(n)) !store exact solutions, for comparison

  do i=1,n
    x = p(i,1)
    y = p(i,2)
    
    psi(i) = sol(i)*sqrt(x)
    psix(i) = sqrt(x)*solx(i)+sol(i)/(2*sqrt(x))
    psiy(i) = sqrt(x)*soly(i)
    psixx(i) = solx(i)/sqrt(x) + sqrt(x)*solxx(i) - sol(i)/(4.0d0*x**(3.0d0/2.0d0))
    psixy(i) = sqrt(x)*solxy(i) + soly(i)/(2*sqrt(x))
    psiyy(i) = sqrt(x)*solyy(i)

    exa(i) = exact(x,y,c,d1,d2,d3,d4)
    ex(i) = exactx(x,y,c,d1,d2,d3,d4)
    ey(i) = exacty(x,y,c,d1,d2,d3,d4)
    exx(i) = exactxx(x,y,c,d1,d2,d3,d4)
    exy(i) = exactxy(x,y,c,d1,d2,d3,d4)
    eyy(i) = exactyy(x,y,c,d1,d2,d3,d4)
  enddo
  
  open(unit=1,file='infiles/h.txt')
    read(1,*) edge
  close(1)

  stdev = 0.0d0
!  open(1,file='paper/triplot/psi.txt')
!    do i=1,n
!      write(1,*) psi(i)
!    enddo
!  close(1) 
  
  open(1,file='files/convxx.dat',position='append')
    stdev = rell2(psixx,exx,areas)
    write(1,*)  stdev,edge!int(sqrt(float(nt)))
  close(1) 

  open(1,file='files/convyy.dat',position='append')
    stdev = rell2(psiyy,eyy,areas)
    write(1,*)  stdev,edge!int(sqrt(float(nt)))
  close(1) 

  open(1,file='files/convxy.dat',position='append')
    stdev = rell2(psixy,exy,areas)
    write(1,*)  stdev,edge!int(sqrt(float(nt)))
  close(1) 

  open(1,file='files/convsol.dat',position='append')
    stdev = rell2(psi,exa,areas)
    write(1,*)  stdev,edge!int(sqrt(float(nt)))
  close(1) 

  open(1,file='files/convx.dat',position='append')
    stdev = rell2(psix,ex,areas)
    write(1,*)  stdev,edge!int(sqrt(float(nt)))
  close(1) 

  open(1,file='files/convy.dat',position='append')
    stdev = rell2(psiy,ey,areas)
    write(1,*)  stdev,edge!int(sqrt(float(nt)))
  close(1) 

  open(1,file='files/convex.dat',position='append')
    stdev = 0
    do i = 1,n
        stdev = stdev + exx(i)**2
    enddo
    stdev = sqrt(stdev)
    write(1,*)  stdev,edge!int(sqrt(float(nt)))
  close(1) 

  open(1,file='files/nodenom.dat',position='append')
    stdev = 0
    do i = 1,n
        stdev = stdev + (psixx(i)-exx(i))**2
    enddo
    stdev = sqrt(stdev)
    write(1,*)  stdev,edge!int(sqrt(float(nt)))
  close(1) 

endprogram testing

