module curvestuff
contains
        
subroutine gssolve_wrapper(c,p,t,b,srcloc,answer,srcval,areas,bound,order,sol) !As the name suggests, this is a wrapper for the
!        Grad-Shafranov solver routine, to be used in the context of Gradyoupee (that function needs the outputs srcloc, srcval
!        This solves an equation of the form \Delta u(x,y) = F(x,y,u), either through iteration (order=1) or through a direct solve
!        assuming that F(x,y,u) = g(x,y)*u + h(x,y)
  use mesh
  use functions
  implicit none
  real *8,parameter::small = 1e-14 !When we use fixed point iteration (only for the initial solution), we want it to be to this
!  accuracy
  real *8::error,c
  integer:: order,nt,i,n
  real *8, dimension(:),allocatable::guess,answer
  real *8, dimension(:,:)::p
  integer, dimension(:,:)::t
  integer, dimension(:)::b
  real *8, dimension(:,:), allocatable::srcloc
  real *8, dimension(:)::bound
  real *8, dimension(:),optional::sol
  real *8, dimension(:), allocatable::srcval,areas,x
  
  n = size(p(:,1))
  allocate(guess(n),answer(n))
  error = 1000 !we set the error to be big so our while loop will start properly
  if(order==1)then !order=1 is the case for the initial solve, where we have to use iteration
    do i=1,n
      guess(i) = 1 !we set our initial guess for u to be coarse.
    enddo
    do while(error>small)
      if(allocated(srcloc).eqv..true.) then
        deallocate(srcloc,srcval,areas,x) !The gssolve routine takes these as unallocated arrays, so we need to deallocate them
        !every time
      endif
    
      call gssolve(c,p,t,b,srcloc,x,srcval,areas,bound,order,guess) !Solves given the guess, and the boundary "bound".
    
      error=0 !Initialise the error
      do i=1,n
        error = error + (x(i)-guess(i))**2 !compute the L-2 norm of the difference between the functions
        guess(i) = x(i)
      enddo
      error = sqrt(error)
      write(*,*) "Error = ",error
    enddo
  elseif(order==2)then !Order=2 is when we are solving the first derivatives. In this case we can do a direct solve, as we already
!          have the solution from before (sol). We input this in the guess array, but it's not a guess. It's the actual solution
!          from a previous solve.
    do i=1,n
      guess(i) = sol(i)
    enddo
    call gssolve(c,p,t,b,srcloc,x,srcval,areas,bound,order,guess) !Call the solver
  endif

  answer = x !Spit out the answer
endsubroutine gssolve_wrapper

subroutine dderpois(infi,findif,solx,soly,solxx,solxy,solyy,sol,p,t,b,areas,ux_old,uy_old,ubx,uby,tarc,rarc,upx_old,upy_old) !Master subroutine. Computes G-S solution, and its
!       derivatives. !REMOVE UN,UBN AFTER DEBUGGING
  use mesh
  use functions
  implicit none
  real *8,dimension(:,:),allocatable::gn,p
  complex *16, dimension(:,:),allocatable::tran
  real *8,dimension(:),allocatable::tarc,uh,xin,yin,dx,dy,ddx,ddy,rarc,upx,upy,uhn,uxt,ubxx,ubxy
  real *8,dimension(:),allocatable::un,upn,ux,uy,ubx,uby,sol,solx,soly,areas,solxx,solyy,solxy,x,ubn
  real *8, dimension(:),allocatable::ux_old,uy_old,upx_old,upy_old
  real *8::d1,d2,d3,d4,c,pi,ds,eps,del,kap,l,infi,findif,gam
  real *8,dimension(2)::that,nhat,der,dder
  real *8,dimension(2,2)::flipmat
  real *8,dimension(10)::args
  real *8::det,temp
  integer,dimension(:),allocatable::b
  integer,dimension(:,:),allocatable::t
  integer::i,j,k,n,m,bsize,ssize
  complex *16,dimension(:),allocatable::cxarr

  write(*,*) "Starting solve"

  call distmesh(p,t,b,d1,d2,d3,d4,c,gam) !we import the arrays describing the finite element decomposition of the tokamak
  args = (/d1,d2,d3,d4,0.3d0,1.0d0*1e-14,1.0d0,0.0d0,c,gam/)!arguments for findr
  
  write(*,*) "Generating boundary FFT"

  call fftgen(10000,args,tran)!generate the fft array
  bsize = size(b)
  n = 3*size(b) !we want 3x oversampling in our boundary
  pi = 4.0d0*atan(1.0d0) !FORTRAN needs us to tell is what pi is, even though it knows the arctan() function...
  allocate(gn(n,n))
  allocate(xin(n),yin(n),dx(n),dy(n),ddx(n),ddy(n))
  allocate(rarc(n),uh(n),cxarr(n),uhn(n),upn(n),ubx(bsize),uby(bsize),uxt(n),ubxx(bsize),ubxy(bsize)) !REMOVE UBN
  write(*,*) "Generating boundary parametrisation"
  call arcparam(0.0d0,2.0d0*pi,tarc,ds,n,l,d1,d2,d3,d4,c,gam,infi,findif,tran) !generate a parametrisation of the boundary. tarc is the array
! of angles which give equidistant points along the boundary
  do i = 1,n
    der = dxdy(tarc(i),d1,d2,d3,d4,c,gam,findif,infi,tran) !array of first derivatives at those points
    dder = ddxddy(tarc(i),d1,d2,d3,d4,c,gam,findif,infi,tran) !array of second derivatives
    rarc(i) = findr_fft(tarc(i),tran) !array of radii away from the center of the tokamak (1,0)
    xin(i) = 1.0d0 + rarc(i)*cos(tarc(i)) !x coordinates
    yin(i) = rarc(i)*sin(tarc(i))! y coordinates
    dx(i) = der(1) !put first derivatives in a size 2 array
    dy(i) = der(2)
    ddx(i) = dder(1) !same for second derivatives
    ddy(i) = dder(2)
  enddo

  call getgnmat(gn,xin,yin,dx,dy,ddx,ddy,n) !as the name says, solves for the G matrix (described in the paper) 

  !CHANGE un,ubn BACK AFTER DEBUGGING

  call derpois(d1,d2,d3,d4,c,gam,infi,findif,solx,soly,ux,uy,sol,p,t,b,rarc,tarc,&
      xin,yin,dx,dy,ddx,ddy,gn,tran,ubx,uby,ds,l,un,ubn,upx_old,upy_old) !Master
!  first-derivative subroutine. Spits out the initial solution and the first derivatives

  allocate(ux_old(size(ux)),uy_old(size(uy)))
  ux_old = ux
  uy_old = uy

  ssize = size(sol)
  allocate(solyy(ssize))
  !Solve for solxx,solxy
  do i =1,n
    cxarr(i) = cmplx(ux(i),0.0D00,kind=16) !We need to store arrays inside complex arrays for the spectral derivative code to take
!    them in properly
  enddo
  call specder(0.00d0,l,n,cxarr,uxt) !Take spectral derivative of the x derivative boundary to get the tangential component of
!  second derivative boundary

  write(*,*) "Solving with first derivative boundary..." 

  call gradyoupee(upx,upy,d1,d2,d3,d4,c,gam,p,t,b,tarc,n,x,infi,findif,tran,areas,ubx,2,sol) !Get gradient of u^p
  call modsolveyouh(gn,xin,yin,dx,dy,upx,upy,uh,n,ds,uxt)!Modified version of the code which solves the integral equation for U^h,
!  which takes in the tangential derivative on the boundary, since that modifies the answer. 

  do i =1,n
    cxarr(i) = cmplx(uh(i),0.0D00,kind=16)
  enddo

  call specder(0.0d0,l,n,cxarr,uhn) !spectral derivative of U^h gives us U^h_t, which is equal to u^h_n

  do i = 1,n
    nhat = (/(-1.0d0)*dy(i)/sqrt(dx(i)**2+dy(i)**2),dx(i)/sqrt(dx(i)**2+dy(i)**2)/)
    that = (/dx(i)/sqrt(dx(i)**2+dy(i)**2),dy(i)/sqrt(dx(i)**2+dy(i)**2)/)
    upn(i) = upx(i)*nhat(1)+upy(i)*nhat(2)!dotting the gradient with nhat gives normal derivative
    un(i) = uhn(i) + upn(i)
    det = nhat(1)*that(2)-nhat(2)*that(1) !convert from tangential coords to x,y
    ux(i) = 1.0d0/det*(un(i)*that(2)-uxt(i)*nhat(2)) 
    uy(i) = 1.0d0/det*(uxt(i)*nhat(1)-un(i)*that(1))
  enddo
  
  do i = 1,bsize !we linearly interpolate (along theta) the values of ux and uy on the boundary to the vertices of the relevant triangles
    temp = atan2(p(b(i),2),p(b(i),1)-1.0d0) !find the angle at which point i is along the boundary
    if(temp<0) then
        temp = 2.0d0*pi+temp !if angle is negative, express in between 0 and 2pi instead
    endif
    j = upper(temp,tarc) !in our array of angles, find the index which is right above our point
    k = lower(temp,tarc) !find the one right below
    ubxx(i) = interp1d(temp,tarc(k),tarc(j),ux(k),ux(j)) !interpolate x derivative boundary
    ubxy(i) = interp1d(temp,tarc(k),tarc(j),uy(k),uy(j)) !same for y
  enddo
 
  write(*,*) ('Taking second derivatives...')
  call weirdder(c,solxx,p,t,b,ubxx,2,sol,solx)!Direct-solve for the derivatives
  call weirdder(c,solxy,p,t,b,ubxy,3,sol,soly)!

  !solve solyy
  do i=1,ssize
    solyy(i) = gsrhs(sol(i),p(i,1),p(i,2),c) - solxx(i) !we already know the answer for the yy derivative
  enddo
endsubroutine dderpois

subroutine derpois(d1,d2,d3,d4,c,gam,infi,findif,solx,soly,ux,uy,sol,p,t,b,rarc,tarc,xin,yin,dx,dy,ddx,& !REMOVE UN,UBN AFTER DEBUGGING
                ddy,gn,tran,ubx,uby,ds,l,un,ubn,upx,upy) !solves poisson equation with first derivatives to second order error.
  use mesh
  use functions
  implicit none
  real *8,dimension(:,:)::gn,p
  complex *16,dimension(:,:)::tran
  real *8,dimension(:)::tarc,rarc,dx,dy,ddx,ddy
  real *8,dimension(:),allocatable::uh,xin,yin,upx,upy,uhn,un,upn,ux,uy,ubx,uby,sol,solx,soly,areas,bound,ubn
  real *8::d1,d2,d3,d4,c,gam,pi,ds,eps,del,kap,l,infi,findif
  real *8,dimension(2)::that,nhat,der,dder
  real *8,dimension(2,2)::flipmat
  real *8,dimension(10)::args
  real *8::det,temp
  integer,dimension(:)::b
  integer,dimension(:,:)::t
  integer::i,j,k,n,m,bsize
  complex *16,dimension(:),allocatable::cxarr

  bsize = size(b)
  n = 3*size(b) !we want the size of our edge decomposition to be comparable to that of the fem, but maybe more accurate
  pi = 4.0d0*atan(1.0d0)

  allocate(uh(n),cxarr(n),uhn(n),un(n),upn(n),ux(n),uy(n),bound(bsize),ubn(bsize))
  
  do i=1,bsize
    bound(i) = 0.0d0
  enddo

  write(*,*) "Initial solve. Iterating...."

  call gradyoupee(upx,upy,d1,d2,d3,d4,c,gam,p,t,b,tarc,n,sol,infi,findif,tran,areas,bound,1) !we have the gradient of u^p. 
  call solveyouh(gn,xin,yin,dx,dy,upx,upy,uh,n,ds) ! solves for u^h  
  do i =1,n
    cxarr(i) = cmplx(uh(i),0.0D00,kind=16)
  enddo
  
  call specder(0.0d0,l,n,cxarr,uhn) !spectral derivative of U^h gives us U^h_t, which is equal to little_u^h_n

  do i = 1,n
    nhat = (/(-1.0d0)*dy(i)/sqrt(dx(i)**2+dy(i)**2),dx(i)/sqrt(dx(i)**2+dy(i)**2)/)
    that = (/dx(i)/sqrt(dx(i)**2+dy(i)**2),dy(i)/sqrt(dx(i)**2+dy(i)**2)/)
    upn(i) = upx(i)*nhat(1)+upy(i)*nhat(2)!dotting the gradient with nhat gives normal derivative
    un(i) = uhn(i) + upn(i)
    det = nhat(1)*that(2)-nhat(2)*that(1) !same for tangential derivative
    ux(i) = 1.0d0/det*(un(i)*that(2)-0*nhat(2)) !the zero comes from the fact that we know u_t to be 0
    uy(i) = 1.0d0/det*(0*nhat(1)-un(i)*that(1))
  enddo
  
  do i = 1,bsize !we linearly interpolate (along theta) the values of ux and uy on the boundary to the vertices of the relevant triangles
    temp = atan2(p(b(i),2),p(b(i),1)-1.0d0) !find the angle at which point i is along the boundary
    if(temp<0) then
        temp = 2.0d0*pi+temp !if angle is negative, express in between 0 and 2pi instead
    endif
    j = upper(temp,tarc) !in our array of angles, find the index which is right above our point
    k = lower(temp,tarc) !find the one right below
    
    ubn(i) = interp1d(temp,tarc(k),tarc(j),un(k),un(j)) !REMOVE AFTER DEBUGGING
    ubx(i) = interp1d(temp,tarc(k),tarc(j),ux(k),ux(j)) !interpolate x derivative boundary
    uby(i) = interp1d(temp,tarc(k),tarc(j),uy(k),uy(j)) !same for y
  enddo

  write(*,*) ('Taking first derivatives...')
  
  call weirdder(c,solx,p,t,b,ubx,0,sol)
  call weirdder(c,soly,p,t,b,uby,1,sol)
end subroutine derpois

subroutine modsolveyouh(gn,xin,yin,dx,dy,upx,upy,uh,n,ds,uxt) !solves linear system for u^h
  implicit none
  integer::n,i,j,info
  integer,dimension(n)::pvt
  real *8::ds,norm
  real *8,dimension(n)::xin,yin,dx,dy,rnx,rny,upx,upy,uh,ones,rhs,uxt
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
    sigma(i) = cmplx(upx(i)*ut(1)+upy(i)*ut(2)-uxt(i),0.0D0,kind=16)
    do j=1,n !this nested loop builds the left hand side matrix, which should be 1/2*eye(n) + ds*g + ds^2*ones(n,n)
      if(i==j) then
        lhs(i,j) = 0.5d0 + ds*gn(i,j) + ds**2
      else
        lhs(i,j) = ds*gn(i,j) + ds**2
      endif
    enddo
  enddo

  call l2dacquadwrapl(xin,yin,ones,rnx,rny,ds,n,1,sigma,0,mu,2,1,1,-1,pot,potn,grad) !call integration routine

  do i = 1,n
    rhs(i) = (-1.0d0)*real(pot(i))
    pvt(i) = 0
  enddo

  call dgesv(n,1,lhs,n,pvt,rhs,n,info) !solve linear system
    
  do i = 1,n
    uh(i) = rhs(i)
  enddo
endsubroutine modsolveyouh

subroutine solveyouh(gn,xin,yin,dx,dy,upx,upy,uh,n,ds) !solves linear system for u^h
  implicit none
  integer::n,i,j,info
  integer,dimension(n)::pvt
  real *8::ds,norm
  real *8,dimension(n)::xin,yin,dx,dy,rnx,rny,upx,upy,uh,ones,rhs
  real *8,dimension(n,n)::lhs,gn
  real *8,dimension(2)::that
  !real *8,dimension(2,n)::that
  complex *16,dimension(n)::pot,potn,sigma,mu
  complex *16,dimension(n,2)::grad  

  !ones(:) = 1.0d0
  !norm(:) = sqrt(dx**2+dy**2)
  !rnx(:) = dy(:)/norm(:)
  !rny(:) = -dx(:)/norm(:)
  !that(1,:) = dx(:)/norm(:)
  !that(2,:) = dy(:)/norm(:)
  !mu(:) = cmplx(0,kind=16)
  !sigma(:) = cmplx(upx(:)*that(1,:)+upy(:)*that(2,:),kind=16)
  !lhs(:,:) = ds*gn(:,:)+ds**2
  !do i = 1,n
  !do j = 1,n
  !      lhs(i,j) = lhs(i,j) +0.5d0
  !enddo
  !enddo

  do i = 1,n !in this loop, we build all the arrays necessary for l2dacquadwrapl to run
    ones(i) = 1.0d0 
    norm = sqrt(dx(i)**2+dy(i)**2)
    rnx(i) = (1.0d0)*dy(i)/norm !x-component of normal derivative vector
    rny(i) = (-1.0d0)*dx(i)/norm !y-component
    that = (/dx(i)/norm,dy(i)/norm/) !tangent unit vector
    mu(i) = cmplx(0.0D0,0.0D0,kind=16) !the mu element in this integration is zero
    sigma(i) = cmplx(upx(i)*that(1)+upy(i)*that(2),kind=16)
    do j=1,n !this nested loop builds the left hand side matrix, which should be 1/2*eye(n) + ds*g + ds^2*ones(n,n)
      if(i==j) then
        lhs(i,j) = 0.5d0 + ds*gn(i,j) + ds**2
      else
        lhs(i,j) = ds*gn(i,j) + ds**2
      endif
    enddo
  enddo

  call l2dacquadwrapl(xin,yin,ones,rnx,rny,ds,n,1,sigma,0,mu,2,1,1,-1,pot,potn,grad) !call integration routine

  do i = 1,n
    rhs(i) = (-1.0d0)*real(pot(i))
    pvt(i) = 0
  enddo

  call dgesv(n,1,lhs,n,pvt,rhs,n,info) !solve linear system
    
  do i = 1,n
    uh(i) = rhs(i)
  enddo
endsubroutine solveyouh

subroutine getgnmat(gn,xin,yin,dx,dy,ddx,ddy,n)! This is a purely geometric definition for the G matrix. It is all described in the
!        paper 
    implicit none
    integer::n,i,j,k
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
      gn(i,i) = (ddx(i)*dy(i) - ddy(i)*dx(i))/(4.0d0*pi*(dx(i)**2+dy(i)**2)**(3.0d0/2.0d0))
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

subroutine gradyoupee(upx,upy,d1,d2,d3,d4,c,gam,p,t,b,tarc,m,x,infi,findif,tran,areas,bound,order,sol) !computes u^p on the boundary of the tokamak using qbx-fmm integration methods.
    use mesh
    use functions
    implicit none
    integer::n,m,i,nb,order,d
    complex *16, dimension(:,:)::tran
    real *8, dimension(:,:), allocatable::srcloc,targloc,targnorm
    real *8, dimension(:)::tarc,bound
    real *8, dimension(:),optional::sol
    real *8, dimension(:), allocatable::srcval,psol,x,y,r,upx,upy,areas,upx_test,upy_test,upx_corr,upy_corr
    complex *16, dimension(:), allocatable::pot,pot1,pot2
    complex *16, dimension(:,:),allocatable::grad1,grad2
    real *8:: d1,d2,d3,d4,c,gam,pi,l,infi,findif,mistake_x,mistake_y,edge
    real *8,dimension(2)::der
    real *8, dimension(:,:)::p
    integer, dimension(:,:)::t
    integer, dimension(:)::b

    pi = 4.0d0*atan(1.0d0)
    
    if(order==1)then !order = 1 is the initial solve
      call gssolve_wrapper(c,p,t,b,srcloc,x,srcval,areas,bound,order) !Solve the Grad Shafranov equation 
    elseif(order==2)then ! order = 2 is the solution for first derivatives
      call gssolve_wrapper(c,p,t,b,srcloc,x,srcval,areas,bound,order,sol)
    endif
    n = size(srcval)
    allocate(targloc(2,m),targnorm(2,m),pot(m),r(m),upx(m),upy(m))
    !DEBUG
    allocate(pot1(m),pot2(m),grad1(2,m),grad2(2,m))
    !DEBUG
    do i = 1,m
      r(i) = findr_fft(tarc(i),tran)
      targloc(1,i) = 1.0d0 + r(i)*cos(tarc(i))!targloc is the array of points where we want our integrated grad u^p
      targloc(2,i) = r(i)*sin(tarc(i))
      der = dxdy(tarc(i),d1,d2,d3,d4,c,gam,findif,infi,tran)
      targnorm(1,i) = der(2)/sqrt(der(1)**2+der(2)**2) !the array of normal vectors at those points
      targnorm(2,i) = (-1.0d0)*der(1)/sqrt(der(1)**2+der(2)**2)
    enddo  

   call l2dacquadwrap(srcloc,srcval,targloc,targnorm,n,m,2,-1,pot) !Call the fmm-qbx integrator
!  call l2dacquadwrap(srcloc,srcval,targloc,targnorm,n,m,2,-1,pot1,pot2,grad1,grad2) !Call the fmm-qbx integrator
  do i = 1,m
    upx(i) = (-1.0d0)*real(pot(i)) 
    upy(i) = (1.0d0)*imag(pot(i))
  enddo
!  upx = real(pot1)
!  upy = (-1)*imag(pot1)
  !DEBUGGING
 ! pot = 0
 ! if (order==1)then
 !   allocate(upx_test(m),upy_test(m),upx_corr(m),upy_corr(m))
 !   d = 30
 !   do i = 1,n
 !       srcval(i) = areas(i)*gruptest_rhs(srcloc(1,i),srcloc(2,i),d)
 !   enddo
 !   call l2dacquadwrap(srcloc,srcval,targloc,targnorm,n,m,2,-1,pot)
 !   upx_test = (-1.0d0)*real(pot)
 !   upy_test = imag(pot)
 !   do i = 1,m
 !       upx_corr(i) = gruptest_upx(targloc(1,i),targloc(2,i),d)
 !       upy_corr(i) = gruptest_upy(targloc(1,i),targloc(2,i),d)
 !   enddo
 !   
 !   mistake_x = rell2_notri(upx_test,upx_corr)
 !   mistake_y = rell2_notri(upy_test,upy_corr)
 !   open(1,file='files/upx_test.txt')
 !   do i = 1,m
 !       write(1,*) upx_test(i),upx_corr(i),upy_test(i),upy_corr(i)
 !   enddo
 !   close(1)
 !   open(1,file='infiles/h.txt')
 !       read(1,*) edge
 !   close(1)
 !   open(2,file='files/mistakes.txt',position='append')
 !       write(2,*) edge,mistake_x, mistake_y
 !   close(2)
 ! endif
   !DEBUGGING  
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

subroutine arcparam(a,b,tarc,darc,n,l,d1,d2,d3,d4,c,gam,infi,findif,tran) !provides n evenly spaced points along a curve parametrised by r,theta between theta = a and theta = b
  implicit none
  integer::n,i,j
  complex *16, dimension(:,:)::tran
  real *8 a,b,c,darc,l,tinit,tfguess,tfupdate,currerr,ds,d1,d2,d3,d4,gam,infi,findif
  real *8,dimension(2)::der
  real *8,dimension(:),allocatable::t,w,tarc
  call lgmap(t,w,a,b,1)
  l = 0.0d0
  do i = 1,size(t)
    der = dxdy(t(i),d1,d2,d3,d4,c,gam,findif,infi,tran)
    l = l + sqrt(der(1)**2+der(2)**2)*w(i)
  end do
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
    do while(abs(currerr) > infi) !infi
      deallocate(t,w)
      tfguess = tfupdate
      call lgmap(t,w,tinit,tfguess,0)
      ds = 0.0d0
      do j=1,size(t)
        der = dxdy(t(j),d1,d2,d3,d4,c,gam,findif,infi,tran)
        ds = ds + sqrt(der(1)**2+der(2)**2)*w(j)
      end do
      currerr = ds-darc
      der = dxdy(tfguess,d1,d2,d3,d4,c,gam,findif,infi,tran)
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

function ddxddy(theta,d1,d2,d3,d4,c,gam,infi,rerror,tran)
  implicit none
  real *8::infi,rerror
  complex *16,dimension(:,:)::tran
  real *8 ::theta,d1,d2,d3,d4,c,gam,tp,tm,dx,dy,ddx
  real *8, dimension(2)::ddxddy,vec
  real *8,dimension(10)::args
  args = (/d1,d2,d3,d4,0.3d0,rerror,1.0d0,0.0d0,c,gam/)

  if(infi<1e-8)then
    infi=1e-8
  endif
  
  dx = 0.0d0
  dx = (-1.0d0/280.0d0)*findr_fft(theta+4.0d0*infi,tran)!*cos(theta+4.0d0*infi) 
  dx = dx + (1.0d0/280.0d0)*findr_fft(theta-4.0d0*infi,tran)!*cos(theta-4.0d0*infi)
  dx = dx + (4.0d0/105.0d0)*findr_fft(theta+3.0d0*infi,tran)!*cos(theta+3.0d0*infi)
  dx = dx + (-4.0d0/105.0d0)*findr_fft(theta-3.0d0*infi,tran)!*cos(theta-3.0d0*infi)
  dx = dx +(-1.0d0/5.0d0)*findr_fft(theta+2.0d0*infi,tran)!*cos(theta+2.0d0*infi) 
  dx = dx + (1.0d0/5.0d0)*findr_fft(theta-2.0d0*infi,tran)!*cos(theta-2.0d0*infi)
  dx = dx + (4.0d0/5.0d0)*findr_fft(theta+infi,tran)!*cos(theta+infi)
  dx = dx + (-4.0d0/5.0d0)*findr_fft(theta-infi,tran)!*cos(theta-infi)
  dx = dx/infi

  ddx = 0.0d0
  ddx = (-1.0d0/560.0d0)*findr_fft(theta+4.0d0*infi,tran)!*cos(theta+4.0d0*infi) 
  ddx = ddx + (-1.0d0/560.0d0)*findr_fft(theta-4.0d0*infi,tran)!*cos(theta-4.0d0*infi)
  ddx = ddx + (8.0d0/315.0d0)*findr_fft(theta+3.0d0*infi,tran)!*cos(theta+3.0d0*infi)
  ddx = ddx + (8.0d0/315.0d0)*findr_fft(theta-3.0d0*infi,tran)!*cos(theta-3.0d0*infi)
  ddx = ddx + (-1.0d0/5.0d0)*findr_fft(theta+2.0d0*infi,tran)!*cos(theta+2.0d0*infi) 
  ddx = ddx + (-1.0d0/5.0d0)*findr_fft(theta-2.0d0*infi,tran)!*cos(theta-2.0d0*infi)
  ddx = ddx + (8.0d0/5.0d0)*findr_fft(theta+infi,tran)!*cos(theta+infi)
  ddx = ddx + (8.0d0/5.0d0)*findr_fft(theta-infi,tran)!*cos(theta-infi)
  ddx = ddx + (-205.0d0/72.0d0)*findr_fft(theta,tran)!*cos(theta-infi)
  ddx = ddx/(infi*infi)

  ddxddy(1) = ddx*cos(theta)-2*dx*sin(theta)-findr_fft(theta,tran)*cos(theta)
  
  ddxddy(2) = ddx*sin(theta)+2*dx*cos(theta)-findr_fft(theta,tran)*sin(theta)
  
end function ddxddy

function dxdy(theta,d1,d2,d3,d4,c,gam,infi,rerror,tran) !takes dx/dtheta and dy/dtheta @ theta on a tokamak defined by eps,del,kap
  implicit none
  real *8::infi,rerror
  complex *16, dimension(:,:)::tran
  real *8 ::theta,d1,d2,d3,d4,c,gam,dx,dy
  real *8, dimension(2)::dxdy
  real *8, dimension(10)::args
  args = (/d1,d2,d3,d4,0.3d0,rerror,1.0d0,0.0d0,c,gam/)
  
  dx = 0.0d0
  dx = (-1.0d0/280.0d0)*findr_fft(theta+4.0d0*infi,tran)!*cos(theta+4.0d0*infi) 
  dx = dx + (1.0d0/280.0d0)*findr_fft(theta-4.0d0*infi,tran)!*cos(theta-4.0d0*infi)
  dx = dx + (4.0d0/105.0d0)*findr_fft(theta+3.0d0*infi,tran)!*cos(theta+3.0d0*infi)
  dx = dx + (-4.0d0/105.0d0)*findr_fft(theta-3.0d0*infi,tran)!*cos(theta-3.0d0*infi)
  dx = dx +(-1.0d0/5.0d0)*findr_fft(theta+2.0d0*infi,tran)!*cos(theta+2.0d0*infi) 
  dx = dx + (1.0d0/5.0d0)*findr_fft(theta-2.0d0*infi,tran)!*cos(theta-2.0d0*infi)
  dx = dx + (4.0d0/5.0d0)*findr_fft(theta+infi,tran)!*cos(theta+infi)
  dx = dx + (-4.0d0/5.0d0)*findr_fft(theta-infi,tran)!*cos(theta-infi)
  dx = dx/infi

  dxdy(1) = dx*cos(theta)-findr_fft(theta,tran)*sin(theta)
  
  dxdy(2) = dx*sin(theta)+findr_fft(theta,tran)*cos(theta)

end function dxdy

function findr_fft(theta,tran)
  implicit none
  real *8::theta,a
  complex *16, dimension(:,:)::tran
  complex *16::ret
  real *8::findr_fft
  integer::m,i

  m = size(real(tran(:,1)))
  ret = cmplx(0.0D0,0.0D0,kind=16)
  do i=1,m
    a = real(tran(i,1))*theta
    ret = ret + (tran(i,2))*cmplx(cos(a),sin(a),kind=16)
  enddo

  findr_fft = abs(ret)
endfunction findr_fft

subroutine fftgen(n,args,tran)
  use,intrinsic::iso_c_binding
  use mesh
  use functions
  implicit none
  include '/usr/include/fftw3.f03'
  type (c_ptr) :: plan
  integer *4::i,j,m,mid,n
  integer *8::plan_forward,plan2
  integer *8,dimension(:),allocatable::kf,savedk
  real *8, dimension(:), allocatable:: tarc
  complex *16, dimension(:),allocatable::savedr
  complex *16, dimension(:,:),allocatable::tran
  complex *16,dimension(n)::rhat,rarc,rhat2,rarc2
  real *8::pi,r0
  real *8,dimension(10)::args
  pi = 4.0d0*atan(1.0d0)
  tarc = linspace2(0.0d0,2.0d0*pi,int(n,kind=4))

  !call findr n times, store array of values r(theta)
  do i=1,n
    rarc(i) = abs(findr(tarc(i),args))*cmplx(1.0D0,0.0D0,kind=16)
  enddo
  !fourier transform r(theta) = sum(rhat*e^(ik*theta))   
  plan = fftw_plan_dft_1d(n,rarc,rhat,fftw_forward,fftw_estimate)
  call fftw_execute_dft(plan,rarc,rhat)
  call fftw_destroy_plan(plan)
  !change k array, put the second half before the first one
  mid = n/2+1
  allocate(kf(n))

  do i=1,n
!    rhat(i) = rhat(i)
    if(i<=mid+1)then
      kf (i) = i*1.0d0-1.0d0
    else
      kf (i) = i*1.0d0-1.0d0-n*1.0d0 
    endif 
  enddo

  !keep larger modes
  m = 0
  j=0
  do i=1,n
    if(abs(real(rhat(i)))>5e-12)then  
      m = m+1
    endif
  enddo
  allocate(savedr(m),savedk(m),tran(m,2))
  j=1
  do i=1,n
    if(abs(real(rhat(i)))>5e-12)then
      savedr(j) = rhat(i)
      savedk(j) = kf(i)
      j=j+1
    endif
  enddo
  do i=1,m
    tran(i,1) = cmplx(savedk(i),0.0D0,kind=16)
    tran(i,2) = savedr(i)/n
  enddo
endsubroutine fftgen

function findr(theta,args)!d1,d2,d3,theta,guess,errtol,centx,centy) !newton's method on a tokamak defined by eps,del,kap
  implicit none
  integer::i
  real *8,dimension(10)::args
  real *8::d1,d2,d3,d4,c,gam,theta,findr,guess,currerr,errtol,centx,centy,rp,rm,dfr,one,pi
  real *8,dimension(2)::uvec,x,center,xp,xm
  d1 = args(1)
  d2 = args(2)
  d3 = args(3)
  d4 = args(4)
  guess = args(5)
  errtol = args(6)
  centx = args(7)
  centy = args(8)
  c = args(9)
  gam = args(10)
  center(1) = centx
  center(2) = centy
  findr = guess
  
  uvec(1) = cos(theta)
  uvec(2) = sin(theta)
  
  do i = 1,2
    x(i) = findr*uvec(i) + center(i)
  end do 
  
  currerr = tokam(x(1),x(2),d1,d2,d3,d4,c,gam)
  
  do while (abs(currerr) > errtol)
    rp = findr+errtol
    rm = findr-errtol
    do i = 1,2
  xp(i) = rp*uvec(i)
  xm(i) = rm*uvec(i)
    end do
  
    dfr=(tokam(centx+xp(1),centy+xp(2),d1,d2,d3,d4,c,gam)-tokam(centx+xm(1),centy+xm(2),d1,d2,d3,d4,c,gam))/(2.0d0*errtol)
    findr = findr - currerr/dfr    
    do i=1,2
  x(i) = findr*uvec(i) + center(i)
    end do
    currerr = tokam(x(1),x(2),d1,d2,d3,d4,c,gam)
  end do
end function findr

function tokam(x,y,d1,d2,d3,d4,c,gam) !tokamak function
  implicit none
  real *8 :: x,y,d1,d2,d3,d4,c,gam,tokam

!  tokam = c/8.0d0*x**4 + d2*x**2 + d3*(x**4-4.0d0*x**2*y**2)+d1
  tokam = x**4/8+c*(0.5d0*x**2*log(x)-x**4/8)+d1+d2*x**2+d3*(x**4-4*x**2*y**2)+d4*y
end function tokam

subroutine switchpars(ep,de,kap,d1,d2,d3,d4,c,gam) !switches from eps,del,kap to d1,d2,d3
  implicit none
  integer::info
  real *8::ep,de,kap,d1,d2,d3,d4,c,gam,const,one
  real *8,dimension(3,3)::a,at
  real *8,dimension(3)::b,pvt
  one = 1.0d0
  const = (-1.0d0)*(c/8.0d0)

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
function upper(theta,tarc) !finds first index in ordered array tarc which is larger than theta
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
function lower(theta,tarc) !finds last index in ordered array tarc which is smaller than theta
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
end module curvestuff
