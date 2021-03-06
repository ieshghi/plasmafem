program trackmin
  use mesh
  use curvestuff
  implicit none
  integer,parameter :: seed = 86456
  real *8::eps,del,kap,infi,findif,d1,d2,d3,d4,c,halinfi,femerror,edge,val,known1,known2,gam
  real *8::extr(2),naive(2),h,x,y
  real *8, dimension(:),allocatable::solx,soly,sol,areas,solxx,solxy,solyy,errs,psi,psix,psiy,psixx,psixy,psiyy
  real *8, dimension(:,:),allocatable::p,ps
  integer,dimension(:),allocatable::b,bs
  integer,dimension(:,:),allocatable::t,ts
  integer::i,minpos,n
  
  call distmesh(ps,ts,bs,d1,d2,d3,d4,c,gam) !get the parameters of our tokamak
  infi = 1.0d0*1e-14
  halinfi = 1.0d0*1e-8
  findif = 1.0d0*1e-3
  call dderpois(infi,findif,solx,soly,solxx,solxy,solyy,sol,p,t,areas) ! calculate first and second derivatives
  n = size(sol)
  allocate(psi(n),psix(n),psiy(n),psixx(n),psixy(n),psiyy(n))
  do i=1,n
    x = p(i,1)
    y = p(i,2)
    
    psi(i) = sol(i)*sqrt(x)
    psix(i) = sqrt(x)*solx(i)+sol(i)/(2*sqrt(x))
    psiy(i) = sqrt(x)*soly(i)
    psixx(i) = solx(i)/sqrt(x) + sqrt(x)*solxx(i) - sol(i)/(4.0d0*x**(3.0d0/2.0d0))
    psixy(i) = sqrt(x)*solxy(i) + soly(i)/(2*sqrt(x))
    psiyy(i) = sqrt(x)*solyy(i)

  enddo
 
  open(2,file='./infiles/h.txt')
    read(2,*) h
  close(2)


!  extr = newton2d(halinfi,solx,soly,solxx,solyy,solxy,p,t) !2d rootfinding method, returns location of the minimum
  extr = newton2d(halinfi,psix,psiy,psixx,psiyy,psixy,p,t)
  minpos = minloc(sqrt(psix**2+psiy**2),dim=1)
  naive = (/p(minpos,1),p(minpos,2)/)
  open(1,file='track.dat',position='append')
    write(1,*) extr(1),extr(2),naive(1),naive(2),c,h
  close(1)
contains

function newton2d(infi,fx,fy,fxx,fyy,fxy,p,t)
    implicit none
    real *8, dimension(:)::fx,fxx,fy,fxy,fyy
    real *8::p(:,:)
    integer::t(:,:)
    real *8:: infi,error,xguess,xnew,yguess,ynew
    real *8::newton2d(2),jac(2,2)
    integer:: n,nt,i,minpos

    minpos = minloc(abs(fx),dim=1)
    xguess = p(minpos,1)
    yguess = 0.0d0
    error = 100

    do while(error>infi)
        jac(1,1) = interp(p,t,fxx,xguess,yguess,infi)
        jac(2,2) = interp(p,t,fyy,xguess,yguess,infi)
        jac(2,1) = interp(p,t,fxy,xguess,yguess,infi)
        jac(1,2) = jac(2,1)
       
        xnew = xguess - jac(1,1)*interp(p,t,fx,xguess,yguess,infi) - jac(1,2)*interp(p,t,fy,xguess,yguess,infi)
        ynew = yguess - jac(2,1)*interp(p,t,fx,xguess,yguess,infi) - jac(2,2)*interp(p,t,fy,xguess,yguess,infi)
        
        error = sqrt((xnew-xguess)**2 + (ynew-yguess)**2)

        xguess = xnew
        yguess = ynew
    enddo

    newton2d = (/xguess,yguess/)

endfunction newton2d

function inv2(mat)
    implicit none
    real *8::mat(2,2),inv2(2,2),det

    det = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

    inv2(1,1) = 1./det*mat(2,2)
    inv2(2,2) = 1./det*mat(1,1)
    inv2(2,1) = -1./det*mat(2,1)
    inv2(1,2) = -1./det*mat(1,2)
endfunction inv2
    
    
function newton(infi,fx,fxx,p,t)
  implicit none
  real *8, dimension(:)::fx,fxx
  real *8, dimension(:,:)::p
  integer,dimension(:,:)::t
  real *8:: infi,error,xguess,xnew
  real *8::newton
  integer:: n,nt,i,minpos

  minpos = minloc(abs(fx),dim=1)
  xguess = p(minpos,1)
  error = 100
  do while(error>infi)
    xnew = xguess - interp(p,t,fx,xguess,0.0d0,infi)/interp(p,t,fxx,xguess,0.0d0,infi)

    error = abs(xnew-xguess)
    write(*,*) xguess
    xguess = xnew
  enddo  

  newton = xguess
endfunction newton

function interp(p,t,f,x,y,errtol)
  implicit none
  real *8, dimension(:)::f
  real *8, dimension(:,:)::p
  integer,dimension(:,:)::t
  real *8:: x,y,interp,a,b,c,x1,x2,x3,y1,y2,y3,z1,z2,z3,errtol
  integer :: ind,i
  
  ind = 0

  do i=1,size(t(:,1))
    if (checktri(p,t,i,x,y,errtol))then
        ind = i
    endif
  enddo

  x1 = p(t(ind,1),1)
  x2 = p(t(ind,2),1)
  x3 = p(t(ind,3),1)
  y1 = p(t(ind,1),2)
  y2 = p(t(ind,2),2)
  y3 = p(t(ind,3),2)
  z1 = f(t(ind,1))
  z2 = f(t(ind,2))
  z3 = f(t(ind,3))

  a = (y2-y1)*(z3-z1)-(y3-y1)*(z2-z1)
  b = (z2-z1)*(x3-x1)-(z3-z1)*(x2-x1)
  c = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)

  interp = z1-(a*(x-x1)+b*(y-y1))/c
endfunction interp

function checktri(p,t,i,x,y,errtol)
  implicit none
  integer::i
  logical::checktri
  real *8, dimension(:,:)::p
  integer, dimension(:,:)::t
  real *8:: x,y,ax,ay,bx,by,cx,cy,ar1,ar2,ar3,totar,errtol
  
  ax = p(t(i,1),1)
  ay = p(t(i,1),2)
  bx = p(t(i,2),1)
  by = p(t(i,2),2)
  cx = p(t(i,3),1)
  cy = p(t(i,3),2)
  
  ar1 = triarea(x,y,bx,by,cx,cy)
  ar2 = triarea(x,y,ax,ay,cx,cy)
  ar3 = triarea(x,y,ax,ay,bx,by)
  totar = triarea(ax,ay,bx,by,cx,cy)
  
  checktri = ((ar1+ar2+ar3)<=totar+errtol)
endfunction checktri

function triarea(ax,ay,bx,by,cx,cy)
  implicit none
  real *8:: ax,ay,bx,by,cx,cy,triarea,det

  det = abs(((ax-cx)*(by-cy))-((bx-cx)*(ay-cy)))
  triarea = det/2.0d0

endfunction triarea

endprogram trackmin
