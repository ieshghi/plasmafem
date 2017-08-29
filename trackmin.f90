program trackmin
  use mesh
  use curvestuff
  implicit none
  real *8:: eps,del,kap,infi,findif,d1,d2,d3,c,halinfi
  real *8,dimension(2)::extr
  real *8, dimension(:),allocatable::solx,soly,sol,areas,solxx,solxy,solyy
  real *8, dimension(:,:),allocatable::p,ps
  integer,dimension(:),allocatable::b,bs
  integer,dimension(:,:),allocatable::t,ts
  integer::i
  
  call distmesh(ps,ts,bs,eps,del,kap,c) !get the parameters of our tokamak
  infi = 1.0d0*1e-14
  halinfi = 1.0d0*1e-1 
  !ask Antoine
  findif = 1.0d0*1e-3
  call dderpois(infi,findif,solx,soly,solxx,solxy,solyy,sol,p,t,areas) ! calculate first and second derivatives
  call switchpars(eps,del,kap,d1,d2,d3)
  extr = newton2d(halinfi,solx,soly,solxx,solyy,solxy,p,t,areas) !2d rootfinding method, returns location of the minimum
  ! (minx,miny)
  write(*,*) 'Result : x = ',extr(1),', y = ',extr(2)
  write(*,*) 'Expected : x = ',sqrt(2*d2/(c/2 + 4*d3)),', y = ',0
  open(1,file='track.dat',position='append')
    write(1,*) extr(1),extr(2),c
  close(1)
contains

function newton2d(infi,fx,fy,fxx,fyy,fxy,p,t,areas)
  implicit none
  real *8, dimension(:)::fx,fy,fxx,fyy,fxy,areas
  real *8, dimension(:,:)::p
  integer,dimension(:,:)::t
  real *8:: infi,extx,exty,error,xguess,yguess,xnew,ynew
  real *8,dimension(2)::newton2d
  real *8, dimension(2,2):: jac
  real *8, dimension(:,:),allocatable::centr
  integer:: n,nt,i
  
  nt = size(t(:,1))
  allocate(centr(nt,2))
  do i=1,nt
    centr(i,1) = (p(t(i,1),1) + p(t(i,2),1) + p(t(i,3),1))/3.0d0
    centr(i,2) = (p(t(i,1),2) + p(t(i,2),2) + p(t(i,3),2))/3.0d0
  enddo

  xguess = 0.5d0
  yguess = 0.5d0
  error = 1000
  do while (error>infi)
    !evaluate jacobian @ (xguess,yguess) (elements of jacobian are fxx,fxy,fyx,fyy)
    jac(1,1) = interp2d(p,t,fxx,xguess,yguess,centr,areas)
    jac(2,2) = interp2d(p,t,fyy,xguess,yguess,centr,areas)
    jac(1,2) = interp2d(p,t,fxy,xguess,yguess,centr,areas)
    jac(2,1) = interp2d(p,t,fxy,xguess,yguess,centr,areas)
   
    call inv2(jac) 
    !invert jacobian (write explicit 2x2 invertor)

    xnew = xguess - jac(1,1)*interp2d(p,t,fx,xguess, yguess,centr,areas) - jac(1,2)*interp2d(p,t,fy,xguess,yguess,centr,areas)
    ynew = xguess - jac(2,1)*interp2d(p,t,fx,xguess, yguess,centr,areas) - jac(2,2)*interp2d(p,t,fy,xguess,yguess,centr,areas)
    !calculate new xguess_n+1 = xguess_n - inv(jac)*[fx(xguess,yguess) , fy (xguess,yguess)]^T

    error = sqrt((xnew-xguess)**2+(ynew-yguess)**2)
    write(*,*) 'Newton error = ',error
    xguess = xnew
    yguess = ynew
  enddo  
  newton2d = (/xnew,ynew/)
endfunction newton2d

function interp2d(p,t,f,x,y,centr,areas)
  implicit none
  real *8, dimension(:)::f,areas
  real *8, dimension(:,:)::p
  integer,dimension(:,:)::t
  real *8:: x,y,interp2d,d,dcurr,a1,a2,a3
  real *8, dimension(:,:)::centr
  integer :: ind,i
  
  d = 100000
  ind = 0
  do i=1,size(centr(:,1))
    dcurr = sqrt((x-centr(i,1))**2 + (y-centr(i,2))**2)
    if(dcurr<d)then
      ind = i
      d = dcurr
    endif
  enddo

  a1 = tricross(p(t(ind,2),1)-p(t(ind,3),1),p(t(ind,2),2)-p(t(ind,3),2),x-p(t(ind,2),1),y-p(t(ind,2),2))
  a2 = tricross(p(t(ind,1),1)-p(t(ind,3),1),p(t(ind,1),2)-p(t(ind,3),2),x-p(t(ind,1),1),y-p(t(ind,1),2))
  a3 = tricross(p(t(ind,2),1)-p(t(ind,1),1),p(t(ind,2),2)-p(t(ind,1),2),x-p(t(ind,2),1),y-p(t(ind,2),2))

  interp2d = 1.0d0/areas(ind)*(f(t(ind,1))*a1 + f(t(ind,2))*a2 + f(t(ind,3))*a3)
endfunction interp2d

function tricross(ux,uy,vx,vy)
  implicit none
  real *8::ux,uy,vx,vy,tricross

  tricross = 0.5d0*abs(ux*vy-vx*uy)
endfunction tricross

subroutine inv2(mat)
  implicit none
  real *8,dimension(2,2)::mat,newmat
  real *8::det

  det = mat(1,1)*mat(2,2)-mat(2,1)*mat(1,2)

  newmat(1,1) = (-1.0d0)/det*mat(1,1)
  newmat(2,2) = (-1.0d0)/det*mat(2,2)
  newmat(1,2) = 1.0d0/det*mat(2,1)
  newmat(2,1) = 1.0d0/det*mat(1,2)

  mat = newmat
endsubroutine inv2


endprogram trackmin
