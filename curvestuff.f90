module curvestuff
contains

function upper(theta,tarc) !this function needs to be checked! not sure if it's doing its work
  implicit none
  real *8,dimension(:)::tarc
  real *8::theta
  integer::n,upper,i
  n = size(tarc)
  upper = 0
  do i = 1,n
    if(tarc(i)>theta .and. upper==0) then
  upper=i
    endif
  enddo
  if (upper==0) then
    upper = 1
  endif
end function upper
function lower(theta,tarc) !this function needs to be checked! not sure if it's doing its work
  implicit none
  real *8,dimension(:)::tarc
  real *8::theta
  integer::n,lower,i
  n = size(tarc)
  lower = 0
  do i = 1,n
    if(tarc(i)<theta) then
  lower=i
    endif
  enddo
  if (lower==0) then
    lower = 1
  endif
end function lower

subroutine derpois(d1,d2,d3,infi,findif,solx,soly,sol,p,t,b,ubx,uby) !solves poisson equation with first derivatives to second order error.
  !also important (less so) once done debugging, outputting ubx,uby, and b is unnecessary
  use mesh
  implicit none
  real *8,dimension(:,:),allocatable::gn,p,tran
  real *8,dimension(:),allocatable::tarc,uh,xin,yin,dx,dy,ddx,ddy,rarc,upx,upy,uhn,un,upn,ux,uy,ubx,uby,sol,solx,soly
  !for debugging
  real *8,dimension(:),allocatable::fux,fuy,fun,fut
  !\for debugging
  real *8::d1,d2,d3,pi,ds,eps,del,kap,l,infi,findif
  real *8,dimension(2)::that,nhat,der,dder
  real *8,dimension(2,2)::flipmat
  real *8,dimension(7)::args
  real *8::det,temp
  integer,dimension(:),allocatable::b
  integer,dimension(:,:),allocatable::t
  integer::i,j,k,n,m,bsize
  complex *16,dimension(:),allocatable::cxarr

  call distmesh(p,t,b,eps,del,kap) !we import the arrays describing the finite element decomposition of the tokamak
  call switchpars(eps,del,kap,d1,d2,d3)
  args = (/d1,d2,d3,0.7d0,infi,1.0d0,0.0d0/)!arguments for findr
  call fftgen(5000,args,tran)!generate the fft file
  bsize = size(b)
  n = 6*size(b) !we want the size of our edge decomposition to be comparable to that of the fem, but maybe more accurate
  pi = 4.0d0*atan(1.0d0)

  allocate(xin(n),yin(n),dx(n),dy(n),ddx(n),ddy(n),gn(n,n))
  !debugging
  allocate(fux(n),fuy(n),fun(n),fut(n))
  !\debugging
  allocate(rarc(n),uh(n),cxarr(n),uhn(n),un(n),upn(n),ux(n),uy(n),ubx(bsize),uby(bsize))
  call arcparam(0.0d0,2.0d0*pi,tarc,ds,n,l,d1,d2,d3,infi,findif,tran) !generate a parametrisation of the boundary. tarc is the array
! of angles which give equidistant points along the boundary

  do i = 1,n
    der = dxdy(tarc(i),d1,d2,d3,findif,infi,tran) !array of first derivatives at those points
    dder = ddxddy(tarc(i),d1,d2,d3,findif,infi,tran) !array of second derivatives
    rarc(i) = findr_fft(tarc(i),tran) !array of radii away from the center of the tokamak (1,0)
    xin(i) = 1.0d0 + rarc(i)*cos(tarc(i)) !x coordinates
    yin(i) = rarc(i)*sin(tarc(i))! y coordinates
    dx(i) = der(1) !put first derivatives in a size 2 array
    dy(i) = der(2)
    ddx(i) = dder(1) !same for second derivatives
    ddy(i) = dder(2)
  enddo

  call getgnmat(gn,xin,yin,dx,dy,ddx,ddy,n) !as the name says, solves for g_n
  call gradyoupee(upx,upy,d1,d2,d3,ds,n,m,sol,infi,findif,tran) !we have the gradient of u^p. 
  call solveyouh(gn,xin,yin,dx,dy,upx,upy,uh,n,ds) ! solves for u^h

  do i =1,n
    cxarr(i) = cmplx(uh(i),0.0d0,kind=16)
  enddo
  
  open(1,file='st.txt') !for debugging

  call specder(0.0d0,2.0d0*pi,n,cxarr,uhn) !spectral derivative of u^h gives us u^h_t, which is equal to u^h_n

  do i = 1,n
    nhat = (/(-1.0d0)*dy(i)/sqrt(dx(i)**2+dy(i)**2),dx(i)/sqrt(dx(i)**2+dy(i)**2)/)
    that = (/dx(i)/sqrt(dx(i)**2+dy(i)**2),dy(i)/sqrt(dx(i)**2+dy(i)**2)/)
    upn(i) = upx(i)*nhat(1)+upy(i)*nhat(2)!dotting the gradient with nhat gives normal derivative
    un(i) = uhn(i) + upn(i)
    det = nhat(1)*that(2)-nhat(2)*that(1) !same for tangential derivative
    ux(i) = 1.0d0/det*(un(i)*that(2)-0*nhat(2)) !the zero comes from the fact that we know u_t to be 0
    uy(i) = 1.0d0/det*(0*nhat(1)-un(i)*that(1))
    
    !debugging
    xin(i) = 1.0d0 + rarc(i)*cos(tarc(i)) !x coordinates
    yin(i) = rarc(i)*sin(tarc(i))! y coordinates
    fux(i) = exactx(xin(i),yin(i),d1,d2,d3)
    fuy(i) = exacty(xin(i),yin(i),d1,d2,d3)
    fun(i) = fux(i)*nhat(1)+fuy(i)*nhat(2)
    fut(i) = fux(i)*that(1)+fuy(i)*that(2)
    fut(i) = (fun(i)-upn(i))/uhn(i)

    write(1,*) upn(i),uhn(i)    
    !\debugging
    
  enddo
!  write(1,*) sum(fut)/(max(1,size(fut))), sqrt(float(size(p(:,1))))
  
  do i = 1,bsize !we linearly interpolate (along theta) the values of ux and uy on the boundary to the vertices of the relevant triangles
    temp = atan2(p(b(i),2),p(b(i),1)-1.0d0) !find the angle at which point i is along the boundary
    if(temp<0) then
  temp = 2.0d0*pi+temp !if angle is negative, express in between 0 and 2pi instead
    endif
    j = upper(temp,tarc) !in our array of angles, find the index which is right above our point
    k = lower(temp,tarc) !find the one right below

    ubx(i) = interp1d(temp,tarc(k),tarc(j),ux(k),ux(j)) !interpolate x derivative boundary
    uby(i) = interp1d(temp,tarc(k),tarc(j),uy(k),uy(j)) !same for y

  enddo
  
  write(*,*) ('taking derivatives...')
  call firstder(d1,d2,d3,solx,p,t,b,ubx,0)
  call firstder(d1,d2,d3,soly,p,t,b,uby,1)
  
  close(1)

end subroutine derpois


subroutine solveyouh(gn,xin,yin,dx,dy,upx,upy,uh,n,ds) !solves linear system for u^h
  implicit none
  integer::n,i,j,info
  integer,dimension(n)::pvt
  real *8::ds,norm
  real *8,dimension(n)::xin,yin,dx,dy,rnx,rny,upx,upy,uh,ones,rhs
  real *8,dimension(n,n)::lhs,gn
  real *8,dimension(2)::ut
  complex *16,dimension(n)::pot,potn,sigma,mu
  complex *16,dimension(n,2)::grad  

  do i = 1,n !in this loop, we build all the arrays necessary for l2dacquadwrapl to run
    ones(i) = 1.0d0 
    norm = sqrt(dx(i)**2+dy(i)**2)
    rnx(i) = (1.0d0)*dy(i)/norm !x-component of normal derivative vector
    rny(i) = (-1.0d0)*dx(i)/norm !y-component
    ut = (/dx(i)/norm,dy(i)/norm/) !tangent unit vector
    mu(i) = cmplx(0.0D0,0.0D0,kind=16) !the mu element in this integration is zero
    sigma(i) = cmplx(upx(i)*ut(1)+upy(i)*ut(2),kind=16)
    do j=1,n !this nested loop builds the left hand side matrix, which should be 1/2*eye(n) + ds*g + ds^2*ones(n,n)
  if(i==j) then
    lhs(i,j) = 0.5d0 + ds*gn(i,j) + ds**2
  else
    lhs(i,j) = ds*gn(i,j) + ds**2
  endif
    enddo
  enddo

  call l2dacquadwrapl(xin,yin,ones,rnx,rny,ds,n,1,sigma,0,mu,4,1,4,-1,pot,potn,grad) !call integration routine

  do i = 1,n
    rhs(i) = (-1.0d0)*real(pot(i))
    pvt(i) = 0
  enddo

  call dgesv(n,1,lhs,n,pvt,rhs,n,info) !solve linear system
    
  do i = 1,n
    uh(i) = rhs(i)
  enddo

endsubroutine solveyouh


subroutine getgnmat(gn,xin,yin,dx,dy,ddx,ddy,n) !xin,yin should come from the arcparam subroutine, n is their length
    implicit none
    integer::n,i,j
    real *8, dimension(n,2)::xp
    real *8, dimension(n)::xin,yin,dx,dy,ddx,ddy
    real *8, dimension(n,n)::gn
    real *8, dimension(2)::gxgy,x
    real *8::pi
    pi = 4.0d0*atan(1.0d0)

    do i=1,n
    xp(i,1) = xin(i)
    xp(i,2) = yin(i)
    do j=1,n
      gn(i,j)=0.0d0
    enddo
    enddo
    do i=1,n
    do j=1,n
    x(1) = xin(i)
    x(2) = yin(i)
    gn(i,j) = ((-1.0d0)*dy(j)*gx(x,(/xp(j,1),xp(j,2)/)) &
  + dx(j)*gy(x,(/xp(j,1),xp(j,2)/)))/sqrt(dx(j)**2 + dy(j)**2)
  enddo
     gn(i,i) =  (ddx(i)*dy(i) - ddy(i)*dx(i))/(4.0d0*pi*(dx(i)**2+dy(i)**2)**(3.0d0/2.0d0))
  enddo
end subroutine getgnmat


function gx(x,xp)
    implicit none
    real *8,dimension(2)::x,xp
    real *8::gx,pi
    pi = 4.0d0*atan(1.0d0)
    
    gx =(1.0d0/(2.0d0*pi))*(xp(1)-x(1))*((x(1)-xp(1))**2+(x(2)-xp(2))**2)**(-1)
end function gx

function gy(x,xp)
    implicit none
    real *8,dimension(2)::x,xp
    real *8::gy,pi
    pi = 4*atan(1.0d0)
    
    gy=(1.0d0/(2.0d0*pi))*(xp(2)-x(2))*((x(1)-xp(1))**2+(x(2)-xp(2))**2)**(-1)
end function gy


subroutine gradyoupee(upx,upy,d1,d2,d3,ds,nb,m,x,infi,findif,tran) !computes u^p on the boundary of the tokamak using qbx-fmm integration methods.
    use mesh
    implicit none
    integer::n,m,i,nb
    real *8, dimension(:,:)::tran
    real *8, dimension(:,:), allocatable::srcloc,targloc,targnorm
    real *8, dimension(:), allocatable::srcval,psol,x,y,tarc,r,upx,upy
    complex *16, dimension(:), allocatable::pot
    real *8:: d1,d2,d3,pi,ds,l,infi,findif
    real *8,dimension(7)::args
    real *8,dimension(2)::der

    pi = 4.0d0*atan(1.0d0)
    call poissolve(d1,d2,d3,srcloc,x,srcval)
    n = size(srcval)
    m = nb
    call arcparam(0.0d0,2.0d0*pi,tarc,ds,m,l,d1,d2,d3,infi,findif,tran)
    allocate(targloc(2,m),targnorm(2,m),pot(m),r(m),upx(m),upy(m))

    args = (/d1,d2,d3,0.7d0,infi,1.0d0,0.0d0/)
    do i = 1,m
    r(i) = findr_fft(tarc(i),tran)
    targloc(1,i) = 1.0d0 + r(i)*cos(tarc(i))
    targloc(2,i) = r(i)*sin(tarc(i))
    der = dxdy(tarc(i),d1,d2,d3,findif,infi,tran)
    targnorm(1,i) = der(2)/sqrt(der(1)**2+der(2)**2)
    targnorm(2,i) = (-1.0d0)*der(1)/sqrt(der(1)**2+der(2)**2)
    enddo  

    call l2dacquadwrap(srcloc,srcval,targloc,targnorm,n,m,6,-1,pot)

  do i = 1,m
    upx(i) = (-1.0d0)*real(pot(i)) !this might have to be revised, depending on the convention for the normal direction (in/out)
    upy(i) = (1.0d0)*imag(pot(i))
  enddo

end subroutine gradyoupee

subroutine specder(xmin,xmax,n,input,deriv)  !takes spectral derivatives of order n of a function evaluated at n points.
    use, intrinsic :: iso_c_binding
    implicit none
    include '/usr/include/fftw3.f03'
    type(c_ptr) :: plan
    integer :: n,i
    integer *8::plan_forward,plan_backward
    complex(c_double_complex), dimension(n) :: input,output,input2,output2
    real *8,dimension(n)::x,k,der,deriv
    real *8,parameter::pi =3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862
    real *8::xmin,xmax,dx

    call dfftw_plan_dft_1d_(plan_forward,n, input,output,fftw_forward,fftw_estimate)
    call dfftw_execute_(plan_forward)
    call dfftw_destroy_plan_(plan_forward)

    do i = 1,n/2
    k(i) = i-1.0d0
    end do

    k(n/2+1) = 0.0d0

    do i = n/2+2,n
    k(i) = (-1.0d0)*n+i-1.0d0
    end do

    do i=1,n
    input2(i) =2.0d0*pi/(xmax-xmin)*k(i)*cmplx(0.0D0,1.0D0,kind=16)*output(i)/n
    end do

    call dfftw_plan_dft_1d_(plan_backward,n, input2, output2,fftw_backward,fftw_estimate)
    call dfftw_execute_(plan_backward, input2, output2)
    call dfftw_destroy_plan_(plan_backward)

    do i = 1,n
	  deriv(i) = real(output2(i))
    end do
end subroutine specder

subroutine arcparam(a,b,tarc,darc,n,l,d1,d2,d3,infi,findif,tran) !provides n evenly spaced points along a curve parametrised by r,theta between theta = a and theta = b
  implicit none
  integer::n,i,j
  real *8, dimension(:,:)::tran
  real *8 a,b,darc,l,tinit,tfguess,tfupdate,currerr,ds,d1,d2,d3,infi,findif
  real *8,dimension(2)::der
  real *8,dimension(:),allocatable::t,w,tarc
  call lgmap(t,w,a,b,1)

  l = 0.0d0
  do i = 1,size(t)
    der = dxdy(t(i),d1,d2,d3,findif,infi,tran)
    l = l + sqrt(der(1)**2+der(2)**2)*w(i)
  end do
  write(*,*) l
  darc = l/n
  allocate(tarc(n))
  do i = 1,n
    tarc(i) = 0.0d0
  end do
  tarc(1) = a

  do i = 2,n
    tinit = tarc(i-1)
    tfguess = tinit+darc
    tfupdate = tfguess
    currerr=1.0d0
    do while(abs(currerr) > findif) !infi
  deallocate(t,w)
  tfguess = tfupdate
  call lgmap(t,w,tinit,tfguess,0)
  ds = 0.0d0
  do j=1,size(t)
    der = dxdy(t(j),d1,d2,d3,findif,infi,tran)
    ds = ds + sqrt(der(1)**2+der(2)**2)*w(j)
  end do
  currerr = ds-darc
  der = dxdy(tfguess,d1,d2,d3,findif,infi,tran)
  tfupdate = tfguess - currerr/sqrt(der(1)**2+der(2)**2)
    end do
    tarc(i) = tfguess
  end do
end subroutine arcparam


subroutine lgmap(x,w,a,b,mode) !does 16-point or 1000-point gauss-legendre quadrature, mapped from [-1,1] to [a,b]
  implicit none
  integer::nr,i,j,ios,mode
    integer, parameter :: maxrecs = 1000000
  real *8,dimension(:),allocatable::x,w,xt,wt
  real *8::a,b
  character(len=4) :: filename
    character(len=1) :: junk
  if (mode == 0) then
    filename = 'sixt'
  else if (mode == 1) then
    filename = 'thou'
  else if (mode ==2) then
    filename = 'teth'
  else
    write(*,*) 'mode must be 0 or 1'
  end if


  nr = 0
  open(unit=1,file='legendre/'//filename//'_w.txt')
    do j=1,maxrecs
    read(1,*,iostat=ios) junk
   if (ios /= 0) exit
    if (j == maxrecs) then
    stop
    endif
    nr = nr + 1
    enddo
    rewind(1)
    allocate(x(nr),w(nr),xt(nr),wt(nr))
    do j=1,nr
    read(1,*) w(j)
    end do
    close(1)
  open(unit = 1,file = 'legendre/'//filename//'_x.txt')
    do j=1,nr
  read(1,*) x(j)
    end do
    close(1)
  
  do i = 1,nr
    xt(i) = (a*(1.0d0-x(i))+b*(1.0d0+x(i)))/2.0d0
    wt(i) = ((b-a)/2.0d0)*w(i)
  end do
  do i = 1,nr
    x(i) = xt(i)
    w(i) = wt(i)
  end do
end subroutine lgmap

function ddxddy(theta,d1,d2,d3,infi,rerror,tran)
  implicit none
  real *8::infi,rerror
  real *8,dimension(:,:)::tran
  real *8 ::theta,d1,d2,d3,tp,tm,dx,dy
  real *8, dimension(2)::ddxddy,vec
  real *8,dimension(7)::args
  args = (/d1,d2,d3,0.7d0,rerror,1.0d0,0.0d0/)

  if(infi<1e-8)then
    infi=1e-8
  endif

  !tp = theta + infi
  !tm = theta - infi
  
  dx = 0.0d0
  dx = (-1.0d0/12.0d0)*findr_fft(theta+2.0d0*infi,tran)*cos(theta+2.0d0*infi) 
  dx = dx + (-1.0d0/12.0d0)*findr_fft(theta-2.0d0*infi,tran)*cos(theta-2.0d0*infi)
  dx = dx + (4.0d0/3.0d0)*findr_fft(theta+infi,tran)*cos(theta+infi)
  dx = dx + (4.0d0/3.0d0)*findr_fft(theta-infi,tran)*cos(theta-infi)
  dx = dx + (-5.0d0/2.0d0)*findr_fft(theta,tran)*cos(theta)
  dx = dx/(infi*infi)

  dy = 0.0d0
  dy = (-1.0d0/12.0d0)*findr_fft(theta+2.0d0*infi,tran)*sin(theta+2.0d0*infi) 
  dy = dy + (-1.0d0/12.0d0)*findr_fft(theta-2.0d0*infi,tran)*sin(theta-2.0d0*infi)
  dy = dy + (4.0d0/3.0d0)*findr_fft(theta+infi,tran)*sin(theta+infi)
  dy = dy + (4.0d0/3.0d0)*findr_fft(theta-infi,tran)*sin(theta-infi)
  dx = dx + (-5.0d0/2.0d0)*findr_fft(theta,tran)*sin(theta)
  dy = dy/(infi*infi)
  
!  ddxddy(1) = (findr(tp,args)*cos(tp)-2*findr_fft(theta)*cos(theta)+findr(tm)*cos(tm))/(infi**2)
 !   ddxddy(2) = (findr(tp,args)*sin(tp)+findr(tm)*sin(tm)-2*findr_fft(theta)*sin(theta))/(infi**2)

  ddxddy(1)=dx
  ddxddy(2)=dy

end function ddxddy

function dxdy(theta,d1,d2,d3,infi,rerror,tran) !takes dx/dtheta and dy/dtheta @ theta on a tokamak defined by eps,del,kap
  implicit none
  real *8::infi,rerror
  real *8, dimension(:,:)::tran
  real *8 ::theta,d1,d2,d3,dx,dy
  real *8, dimension(2)::dxdy
  real *8, dimension(7)::args
  args = (/d1,d2,d3,0.7d0,rerror,1.0d0,0.0d0/)
  
  dx = 0.0d0
  dx = (-1.0d0/280.0d0)*findr_fft(theta+4.0d0*infi,tran)!*dcos(theta+4.0d0*infi) 
  dx = dx + (1.0d0/280.0d0)*findr_fft(theta-4.0d0*infi,tran)!*dcos(theta-4.0d0*infi)
  dx = dx + (4.0d0/105.0d0)*findr_fft(theta+3.0d0*infi,tran)!*dcos(theta+3.0d0*infi)
  dx = dx + (-4.0d0/105.0d0)*findr_fft(theta-3.0d0*infi,tran)!*dcos(theta-3.0d0*infi)
  dx = (-1.0d0/5.0d0)*findr_fft(theta+2.0d0*infi,tran)!*dcos(theta+2.0d0*infi) 
  dx = dx + (1.0d0/5.0d0)*findr_fft(theta-2.0d0*infi,tran)!*dcos(theta-2.0d0*infi)
  dx = dx + (4.0d0/5.0d0)*findr_fft(theta+infi,tran)!*dcos(theta+infi)
  dx = dx + (-4.0d0/5.0d0)*findr_fft(theta-infi,tran)!*dcos(theta-infi)
  dx = dx/infi

  dx = dx*cos(theta)-findr_fft(theta,tran)*sin(theta)
  
  dy = 0.0d0
  dy = (-1.0d0/280.0d0)*findr_fft(theta+4.0d0*infi,tran)!*dsin(theta+4.0d0*infi) 
  dy = dy + (1.0d0/280.0d0)*findr_fft(theta-4.0d0*infi,tran)!*dsin(theta-4.0d0*infi)
  dy = dy + (4.0d0/105.0d0)*findr_fft(theta+3.0d0*infi,tran)!*dsin(theta+3.0d0*infi)
  dy = dy + (-4.0d0/105.0d0)*findr_fft(theta-3.0d0*infi,tran)!*dsin(theta-3.0d0*infi)
  dy = (-1.0d0/5.0d0)*findr_fft(theta+2.0d0*infi,tran)!*dsin(theta+2.0d0*infi) 
  dy = dy + (1.0d0/5.0d0)*findr_fft(theta-2.0d0*infi,tran)!*dsin(theta-2.0d0*infi)
  dy = dy + (4.0d0/5.0d0)*findr_fft(theta+infi,tran)!*dsin(theta+infi)
  dy = dy + (-4.0d0/5.0d0)*findr_fft(theta-infi,tran)!*dsin(theta-infi)
  dy = dy/infi
  
  dy = dy*sin(theta)+findr_fft(theta,tran)*cos(theta)

!  dxdy(1) = (findr_fft(theta+infi)*cos(theta+infi)-findr_fft(theta-infi)*cos(theta-infi))/(2*infi)
!  dxdy(2) = (findr_fft(theta+infi)*sin(theta+infi)-findr_fft(theta-infi)*sin(theta-infi))/(2*infi)
  
  dxdy(1) = dx
  dxdy(2) = dy

end function dxdy

function findr_fft(theta,tran)
  implicit none
  real *8::theta,a
  real *8, dimension(:,:)::tran
  real *8::findr_fft,ret
  integer::m,i

  m = size(tran(:,1))  
  ret = 0.0d0
  do i=1,m
    a = tran(i,1)*theta
    ret = ret + (tran(i,2))*cos(a)
  enddo

  findr_fft = ret
endfunction findr_fft

subroutine fftgen(n,args,tran)
  use,intrinsic::iso_c_binding
  use mesh
  implicit none
  include '/usr/include/fftw3.f03'
  type (c_ptr) :: plan
  integer *4::i,j,m,mid,n
  integer *8::plan_forward,plan2
  integer *8,dimension(:),allocatable::kf,savedk
  real *8, dimension(:), allocatable:: tarc
  real *8, dimension(:),allocatable::savedr
  real *8, dimension(:,:),allocatable::tran
  complex *16,dimension(n)::rhat,rarc,rhat2,rarc2
  real *8::pi,r0
  real *8,dimension(7)::args
  pi = 4.0d0*atan(1.0d0)
  tarc = linspace(0.0d0,2.0d0*pi,int(n,kind=4))
 

  !call findr n times, store array of values r(theta)
  do i=1,n
    rarc(i) = findr(tarc(i),args)*cmplx(1.0D0,0.0D0,kind=16)
  enddo
  !fourier transform r(theta) = sum(rhat*e^(ik*theta))   
  plan = fftw_plan_dft_1d(n,rarc,rhat,fftw_forward,fftw_estimate)
  call fftw_execute_dft(plan,rarc,rhat)
  call fftw_destroy_plan(plan)
  !change k array, put the second half before the first one
 
  mid = n/2+1
  allocate(kf(n))

  do i=1,n
    if(i<=mid+1)then
      kf (i) = i*1.0d0-1.0d0
    else
      kf (i) = i*1.0d0-1.0d0-n*1.0d0 
    endif 
  enddo

  !keep larger modes
  m = 0
  do i=1,n
    if(abs(rhat(i))>5e-12)then
      m = m+1
    endif
  enddo
  allocate(savedr(m),savedk(m),tran(m,2))

  j=1
  do i=1,n
    if(abs(rhat(i))>5e-12)then
      savedr(j) = real(rhat(i))
      savedk(j) = kf(i)
      j=j+1
    endif
  enddo
  !store rhat in a file

  open(1,file='files/fft.txt')
  write(1,*) m
  write(*,*) m
  do i=1,m
    tran(i,1) = savedk(i)
    tran(i,2) = savedr(i)/n
    write(1,*) savedk(i),savedr(i)/n
  enddo
  close(1)
endsubroutine fftgen

function findr(theta,args)!d1,d2,d3,theta,guess,errtol,centx,centy) !newton's method on a tokamak defined by eps,del,kap
  implicit none
  integer::i
  real *8,dimension(7)::args
  real *8::d1,d2,d3,theta,findr,guess,currerr,errtol,centx,centy,rp,rm,dfr,one,pi
  real *8,dimension(2)::uvec,x,center,xp,xm
  d1 = args(1)
  d2 = args(2)
  d3 = args(3)
  guess = args(4)
  errtol = args(5)
  centx = args(6)
  centy = args(7)
  one = 1.0d0
  center(1) = centx
  center(2) = centy
  findr = guess
  
  uvec(1) = dcos(theta)
  uvec(2) = dsin(theta)
  
  do i = 1,2
    x(i) = findr*uvec(i) + center(i)
  end do 
  
  currerr = tokam(x(1),x(2),one,d1,d2,d3)
  
  do while (abs(currerr) > errtol)
    rp = findr+errtol
    rm = findr-errtol
    do i = 1,2
  xp(i) = rp*uvec(i)
  xm(i) = rm*uvec(i)
    end do
  
    dfr=(tokam(centx+xp(1),centy+xp(2),one,d1,d2,d3)-tokam(centx+xm(1),centy+xm(2),one,d1,d2,d3))/(2.0d0*errtol)
    findr = findr - currerr/dfr    
    do i=1,2
  x(i) = findr*uvec(i) + center(i)
    end do
    currerr = tokam(x(1),x(2),one,d1,d2,d3)
  end do
end function findr

function tokam(x,y,c,d1,d2,d3) !tokamak function
  implicit none
  real *8 :: x,y,c,d1,d2,d3,tokam

  tokam = c/8.0d0*x**4 + d2*x**2 + d3*(x**4-4.0d0*x**2*y**2)+d1
end function tokam

subroutine switchpars(ep,de,kap,d1,d2,d3) !switches from eps,del,kap to d1,d2,d3
  implicit none
  integer::info
  real *8::ep,de,kap,d1,d2,d3,const,one
  real *8,dimension(3,3)::a,at
  real *8,dimension(3)::b,pvt
  one = 1.0d0
  const = (-1.0d0)*(1.0d0/8.0d0)

  a=reshape((/one,(1.0d0+ep)**2,(1.0d0+ep)**4,one,(1.0d0-ep)**2,(1.0d0-ep)**4,one&
  ,(1.0d0-de*ep)**2,(1.0d0-de*ep)**4-4.0d0*(1.0d0-de*ep)**2*kap**2.0d0*ep**2/),shape(a))
  at = transpose(a)
  b = (/const*(1.0d0+ep)**4,const*(1.0d0-ep)**4,const*(1.0d0-de*ep)**4/)

  call dgesv(3,1,at,3,pvt,b,3,info)
  d1 = b(1)
  d2 = b(2)
  d3 = b(3)

end subroutine switchpars

function interp1d(x,x1,x2,y1,y2)
  implicit none
  real *8::x,x1,x2,y1,y2,interp1d

  interp1d = y1*(1.0d0-(x-x1)/(x2-x1))+y2*((x-x1)/(x2-x1))
end function interp1d
end module curvestuff
