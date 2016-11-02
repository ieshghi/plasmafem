module curvestuff
contains

function upper(theta,tarc) !THIS FUNCTION NEEDS TO BE CHECKED! NOT SURE IF IT'S DOING ITS WORK
	implicit none
	real(kind=8),dimension(:)::tarc
	real(kind=8)::theta
	integer::N,upper,i
	N = size(tarc)
	upper = 0
	do i = 1,n
		if(tarc(i)>theta .and. upper==0) then
			upper=i
		endif
	enddo
	if (upper==0) then
		upper = n
	endif
end function upper
function lower(theta,tarc) !THIS FUNCTION NEEDS TO BE CHECKED! NOT SURE IF IT'S DOING ITS WORK
	implicit none
	real(kind=8),dimension(:)::tarc
	real(kind=8)::theta
	integer::N,lower,i
	N = size(tarc)
	lower = 0
	do i = 1,n
		if(tarc(i)<theta) then
			lower=i
		endif
	enddo
	if (lower==0) then
		lower = n
	endif
end function lower

subroutine derpois(eps,del,kap,infi,solx,soly,sol,p,t,b,ubx,uby) !Solves poisson equation with first derivatives to second order error.
	!IMPORTANT: for the moment, it is necessary to import the correct values of eps,del,kap in derpois and in distmesh.
	!Otherwise we are working with two different tokamaks. This can be fixed later.
	!ALSO IMPORTANT (less so) once done debugging, outputting ubx,uby, and b is unnecessary
	use mesh
	implicit none
        real(kind=8),dimension(:,:),allocatable::Gn,p
        real(kind=8),dimension(:),allocatable::tarc,uh,xin,yin,dx,dy,ddx,ddy,rarc,upx,upy,uhn,un,upn,ux,uy,ubx,uby,sol,solx,soly
        real(kind=8)::pi,ds,eps,del,kap,L,infi
        real(kind=8),dimension(2)::that,nhat,der,dder
	real(kind=8),dimension(2,2)::flipmat
        real(kind=8),dimension(7)::args
	real(kind=8)::det,temp
	integer,dimension(:),allocatable::b
	integer,dimension(:,:),allocatable::t
        integer::i,j,k,N,m
	complex(kind=8),dimension(:),allocatable::cxarr
	
	
	call distmesh(p,t,b,eps,del,kap) !we import the arrays describing the finite element decomposition of the tokamak
	N = 2*size(b) !we want the size of our edge decomposition to be comparable to that of the FEM, but maybe more accurate
        pi = 4*atan(1.0)


        allocate(xin(N),yin(N),dx(N),dy(N),ddx(N),ddy(N),Gn(N,N))
	allocate(rarc(N),uh(N),cxarr(N),uhn(N),un(N),upn(N),ux(N),uy(N),ubx(N/2),uby(N/2))
        args = (/eps,del,kap,real(0.7,kind=8),real(infi,kind=8),real(1.0,kind=8),real(0.0,kind=8)/)
	
        call arcparam(real(0.0,kind=8),2*pi,tarc,ds,N,L,eps,del,kap) !generate a parametrisation of the boundary. Tarc is the array
! of angles which give equidistant points along the boundary

        do i = 1,N
                der = dxdy(tarc(i),eps,del,kap) !array of first derivatives at those points
                dder = ddxddy(tarc(i),eps,del,kap) !array of second derivatives
                rarc(i) = findr_fast(tarc(i),args) !array of radii away from the center of the tokamak (1,0)
                xin(i) = 1.0 + rarc(i)*cos(tarc(i)) !x coordinates
                yin(i) = rarc(i)*sin(tarc(i))! y coordinates
                dx(i) = der(1) !put first derivatives in a size 2 array
                dy(i) = der(2)
                ddx(i) = dder(1) !same for second derivatives
                ddy(i) = dder(2)
        enddo

        call getgnmat(Gn,xin,yin,dx,dy,ddx,ddy,N) !as the name says, solves for G_n
        call gradyoupee(upx,upy,eps,del,kap,ds,N,m,sol) !we have the gradient of U^p. 
	call solveyouh(Gn,xin,yin,dx,dy,upx,upy,uh,N,ds) ! Solves for U^h

	do i =1,N
		cxarr(i) = complex(uh(i),0.0)
	enddo

	call specder(0.0,real(2*pi,kind=4),N,cxarr,uhn) !spectral derivative of U^h gives us U^h_t, which is equal to u^h_n
	do i = 1,N
		nhat = (/(-1)*dy(i)/sqrt(dx(i)**2+dy(i)**2),dx(i)/sqrt(dx(i)**2+dy(i)**2)/)
		that = (/dx(i)/sqrt(dx(i)**2+dy(i)**2),dy(i)/sqrt(dx(i)**2+dy(i)**2)/)
		upn(i) = upx(i)*nhat(1)+upy(i)*nhat(2)!dotting the gradient with nhat gives normal derivative
		un(i) = uhn(i) + upn(i)
		det = nhat(1)*that(2)-nhat(2)*that(1) !same for tangential derivative
		ux(i) = 1.0/det*(un(i)*that(2)-0*nhat(2)) !the zero comes from the fact that we know u_t to be 0
		uy(i) = 1.0/det*(0*nhat(1)-un(i)*that(1))
	enddo

	do i = 1,N/2 !we linearly interpolate (along theta) the values of ux and uy on the boundary to the vertices of the relevant triangles
		temp = atan2(p(b(i),2),p(b(i),1)-1.0d0) !find the angle at which point i is along the boundary
		if(temp<0) then
			temp = 2*pi+temp !if angle is negative, express in between 0 and 2pi instead
		endif
		j = upper(temp,tarc) !In our array of angles, find the index which is right above our point
		k = lower(temp,tarc) !Find the one right below
		ubx(i) = interp1d(temp,tarc(k),tarc(j),ux(k),ux(j)) !interpolate x derivative boundary
		uby(i) = interp1d(temp,tarc(k),tarc(j),uy(k),uy(j)) !same for y
	enddo

	
	write(*,*) ('Taking derivatives...')
	call firstder(solx,p,t,b,ubx)
	call firstder(soly,p,t,b,uby)
	


end subroutine derpois


subroutine solveyouh(Gn,xin,yin,dx,dy,upx,upy,uh,N,ds) !solves linear system for U^h
	implicit none
	integer::N,i,j,info
	integer,dimension(N)::pvt
	real(kind=8)::ds
	real(kind=8),dimension(N)::xin,yin,dx,dy,rnx,rny,upx,upy,uh,ones,zeros,rhs
	real(kind=8),dimension(N,N)::lhs,gn
	real(kind=8),dimension(2)::ut
	complex(kind=8),dimension(N)::pot,potn,sigma,mu
	complex(kind=8),dimension(N,2)::grad	


	do i = 1,N
		ones(i) = 1.0
		zeros(i) = 0.0
		rnx(i) = dy(i)/sqrt(dx(i)**2+dy(i)**2)
		rny(i) = (-1)*dx(i)/sqrt(dx(i)**2+dy(i)**2)
		ut = (/dx(i),dy(i)/)/sqrt(dx(i)**2+dy(i)**2)
		mu(i) = cmplx(0.0,0.0)
		sigma(i) = cmplx(upx(i)*ut(1)+upy(i)*ut(2))
		do j=1,N
			if(i==j) then
				lhs(i,j) = 0.5 + ds*gn(i,j) + ds**2
			else
				lhs(i,j) = ds*gn(i,j) + ds**2
			endif
		enddo
	enddo

	call l2dacquadwrapl(xin,yin,ones,rnx,rny,ds,N,1,sigma,0,mu,4,1,4,-1,pot,potn,grad)


	do i = 1,N
		rhs(i) = (-1)*real(pot(i))
	enddo

	call dgesv(N,1,lhs,N,pvt,rhs,N,info)
       
	

	do i = 1,N
		uh(i) = rhs(i)
	enddo

end subroutine solveyouh


subroutine getgnmat(Gn,xin,yin,dx,dy,ddx,ddy,N) !xin,yin should come from the arcparam subroutine, N is their length
        implicit none
        integer::N,i,j
        real(kind=8), dimension(N,2)::xp
        real(kind=8), dimension(N)::xin,yin,dx,dy,ddx,ddy
        real(kind=8), dimension(N,N)::Gn
        real(kind=8), dimension(2)::gxgy,x
        real(kind=8)::pi
        pi = 4*atan(1.0)

        do i=1,N
                xp(i,1) = xin(i)
                xp(i,2) = yin(i)
        enddo
        do i=1,N
                do j=1,N
                        x(1) = xin(i)
                        x(2) = yin(i)
                        Gn(i,j) = ((-1)*dy(j)*Gx(x,(/xp(j,1),xp(j,2)/)) + dx(j)*Gy(x,(/xp(j,1),xp(j,2)/)))/sqrt(dx(j)**2 + dy(j)**2)
                enddo
               	Gn(i,i) =  (ddx(i)*dy(i) - ddy(i)*dx(i))/(4*pi*(dx(i)**2+dy(i)**2)**(3.0/2))
        enddo
end subroutine getgnmat


function gx(x,xp)
        implicit none
        real(kind=8),dimension(2)::x,xp
        real(kind=8)::gx,pi
        pi = 4*atan(1.0)
        
        gx =(1/(2*pi))*(xp(1)-x(1))*((x(1)-xp(1))**2+(x(2)-xp(2))**2)**(-1)
end function gx

function gy(x,xp)
        implicit none
        real(kind=8),dimension(2)::x,xp
        real(kind=8)::gy,pi
        pi = 4*atan(1.0)
        
        gy=(1/(2*pi))*(xp(2)-x(2))*((x(1)-xp(1))**2+(x(2)-xp(2))**2)**(-1)
end function gy


subroutine gradyoupee(upx,upy,eps,del,kap,ds,nb,m,x) !Computes u^p on the boundary of the tokamak using QBX-FMM integration methods.
        use mesh
        implicit none
        integer::n,m,i,nb
        real(kind=8), dimension(:,:), allocatable::srcloc,targloc,targnorm
        real(kind=8), dimension(:), allocatable::srcval,psol,x,y,tarc,r,upx,upy
        complex(kind=8), dimension(:), allocatable::pot
        real(kind=8):: del,eps,kap,pi,ds,l
        real(kind=8),dimension(7)::args
        real(kind=8),dimension(2)::der

        pi = 4*atan(1.0)
	call poissolve(srcloc,x,srcval)
	n = size(srcval)
        m = 4*nb
        call arcparam(real(0.0,kind=8),2*pi,tarc,ds,m,l,eps,del,kap)
        allocate(targloc(2,m),targnorm(2,m),pot(m),r(m),upx(m),upy(m))

        args = (/eps,del,kap,real(0.7,kind=8),real(1e-10,kind =8),real(1.0,kind=8),real(0.0,kind=8)/)
        do i = 1,m
                r(i) = findr_fast(tarc(i),args)
                targloc(1,i) = 1 + r(i)*cos(tarc(i))
                targloc(2,i) = r(i)*sin(tarc(i))
                der = dxdy(tarc(i),eps,del,kap)
                targnorm(1,i) = der(2)/sqrt(der(1)**2+der(2)**2)
                targnorm(2,i) = (-1)*der(1)/sqrt(der(1)**2+der(2)**2)
        enddo   
        call l2dacquadwrap(srcloc,srcval,targloc,targnorm,n,m,4,-1,pot)

	do i = 1,m
		upx(i) = real(pot(i)) !This might have to be revised, depending on the convention for the normal direction (in/out)
		upy(i) = (-1)*imag(pot(i))
	enddo
end subroutine gradyoupee



SUBROUTINE specder(xmin,xmax,N,input,deriv)  !Takes spectral derivatives of order N of a function evaluated at N points.
        USE, INTRINSIC :: iso_c_binding
        IMPLICIT NONE
        INCLUDE '/usr/include/fftw3.f03'
        type(C_PTR) :: plan
        integer :: N,i
        integer ::M
        complex(C_DOUBLE_COMPLEX), dimension(N) :: input,output,input2,output2
        real(kind=8),dimension(N)::x,k,der,deriv
        real,parameter::pi =3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862
        real::xmin,xmax,dx

	M = N

        dx = (xmax-xmin)/(N)

        do i = 1,N
                x(i) = xmin + dx*(i-1)
        end do

        plan = fftw_plan_dft_1d(N, input,output, FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan, input, output)
        call fftw_destroy_plan(plan)

        do i = 1,N/2
                k(i) = i-1
        end do

        k(N/2+1) = 0.0

        do i = N/2+2,N
                k(i) = (-1)*N+i-1
        end do

        do i=1,N
                input2(i) =2*pi/(xmax-xmin)*k(i)*cmplx(0.0,1.0)*output(i)/N
        end do

        plan = fftw_plan_dft_1d(N, input2, output2, FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan, input2, output2)
        call fftw_destroy_plan(plan)

        do i = 1,N
		deriv(i) = real(output2(i))
        end do
END SUBROUTINE specder


subroutine arcparam(a,b,tarc,darc,N,L,eps,del,kap) !Provides N evenly spaced points along a curve parametrised by r,theta between theta = a and theta = b
	implicit none
	integer::N,i,j
	real(kind=8) a,b,darc,L,tinit,tfguess,tfupdate,currerr,ds,eps,del,kap
	real(kind=8),dimension(2)::der
	real(kind=8),dimension(:),allocatable::t,w,tarc
	call lgmap(t,w,a,b,1)
	L = 0
	do i = 0,size(t)
		der = dxdy(t(i),eps,del,kap)
		L = L + sqrt(der(1)**2+der(2)**2)*w(i)
	end do
	darc = L/N
	allocate(tarc(N))
	do i = 1,N
		tarc(i) = 0
	end do
	tarc(1) = a

	do i = 2,N
		tinit = tarc(i-1)
		tfguess = tinit+darc
		tfupdate = tfguess
		currerr=1
		do while(abs(currerr) > 1e-6)
			deallocate(t,w)
			tfguess = tfupdate
			call lgmap(t,w,tinit,tfguess,0)
			ds = 0
			do j=1,size(t)
				der = dxdy(t(j),eps,del,kap)
				ds = ds + sqrt(der(1)**2+der(2)**2)*w(j)
			end do
			currerr = ds-darc
			der = dxdy(tfguess,eps,del,kap)
			tfupdate = tfguess - currerr/sqrt(der(1)**2+der(2)**2)
		end do
		tarc(i) = tfguess
	end do

end subroutine arcparam


subroutine lgmap(x,w,a,b,mode) !does 16-point or 1000-point gauss-legendre quadrature, mapped from [-1,1] to [a,b]
	implicit none
	integer::NR,i,J,ios,mode
        INTEGER, PARAMETER :: maxrecs = 1000000
	real(kind=8),dimension(:),allocatable::x,w,xt,wt
	real(kind=8)::a,b
	CHARACTER(LEN=4) :: filename
        CHARACTER(LEN=1) :: junk
	if (mode == 0) then
		filename = 'sixt'
	else if (mode == 1) then
		filename = 'thou'
	else
		write(*,*) 'mode must be 0 or 1'
	end if


	NR = 0
	OPEN(UNIT=1,FILE='legendre/'//filename//'_w.txt')
                DO J=1,maxrecs
                        READ(1,*,IOSTAT=ios) junk
 			IF (ios /= 0) EXIT
                        IF (J == maxrecs) THEN
                                STOP
                        ENDIF
                        NR = NR + 1
                ENDDO
                REWIND(1)
		ALLOCATE(x(NR),w(NR),xt(NR),wt(NR))
                DO J=1,NR
                        READ(1,*) w(j)
                END DO
                CLOSE(1)
	open(unit = 1,file = 'legendre/'//filename//'_x.txt')
		do J=1,NR
			read(1,*) x(J)
		end do
		close(1)
	
	do i = 1,NR
		xt(i) = (a*(1-x(i))+b*(1+x(i)))/2
		wt(i) = ((b-a)/2)*w(i)
	end do
	do i = 1,NR
		x(i) = xt(i)
		w(i) = wt(i)
	end do
end subroutine lgmap

function ddxddy(theta,eps,del,kap)
	implicit none
	real(kind=8),parameter::infi = 1e-5
        real(kind=8) ::theta,eps,del,kap,tp,tm
        real(kind=8), dimension(2)::ddxddy,vec
        real(kind=8),dimension(7)::args
        args = (/eps,del,kap,real(0.7,kind=8),real(infi,kind=8),real(1.0,kind=8),real(0.0,kind=8)/)

	tp = theta + infi
	tm = theta - infi

	ddxddy(1) = (findr_fast(tp,args)*cos(tp)-2*findr_fast(theta,args)*cos(theta)+findr_fast(tm,args)*cos(tm))/(infi**2)
        ddxddy(2) = (findr_fast(tp,args)*sin(tp)+findr_fast(tm,args)*sin(tm)-2*findr_fast(theta,args)*sin(theta))/(infi**2)

end function ddxddy

function dxdy(theta,eps,del,kap) !takes dx/dtheta and dy/dtheta @ theta on a tokamak defined by eps,del,kap
	implicit none
	real(kind=8),parameter::infi = 1e-8
	real(kind=8) ::theta,eps,del,kap
	real(kind=8), dimension(2)::dxdy
	real(kind=8),dimension(7)::args
	

	args = (/eps,del,kap,real(0.7,kind=8),real(infi,kind=8),real(1.0,kind=8),real(0.0,kind=8)/)
	
	dxdy(1) = (findr_fast(theta+infi,args)*cos(theta+infi)-findr_fast(theta-infi,args)*cos(theta-infi))/(2*infi)
	dxdy(2) = (findr_fast(theta+infi,args)*sin(theta+infi)-findr_fast(theta-infi,args)*sin(theta-infi))/(2*infi)

end function dxdy

function findr_fast(theta,args) !easier to call that findr itself, takes an array instead of a bunch of individual args
	implicit none
	real(kind=8),dimension(7)::args
	real(kind=8)::eps,del,kap,guess,theta,errtol,centx,centy
	real(kind=8)::findr_fast
	eps = args(1)
	del = args(2)
	kap = args(3)
	guess = args(4)
	errtol = args(5)
	centx = args(6)
	centy = args(7)

	findr_fast = findr(eps,del,kap,theta,guess,errtol,centx,centy)

end function findr_fast

function findr(eps,del,kap,theta,guess,errtol,centx,centy) !newton's method on a tokamak defined by eps,del,kap
	implicit none
	integer::i
	real(kind=8)::eps,del,kap,theta,findr,guess,currerr,errtol,centx,centy,rp,rm,dfr,one,pi
	real(kind=8),dimension(2)::uvec,x,center,xp,xm
	one = 1.0
	center(1) = centx
	center(2) = centy
	findr = guess
	
	uvec(1) = cos(theta)
	uvec(2) = sin(theta)
	
	do i = 1,2
		x(i) = findr*uvec(i) + center(i)
	end do 
	
	currerr = tokam(x(1),x(2),one,eps,del,kap)
	
	do while (abs(currerr) > errtol)
		rp = findr+errtol
		rm = findr-errtol
		do i = 1,2
			xp(i) = rp*uvec(i)
			xm(i) = rm*uvec(i)
		end do
			
		dfr=(tokam(centx+xp(1),centy+xp(2),one,eps,del,kap)-tokam(centx+xm(1),centy+xm(2),one,eps,del,kap))/(2*errtol)
		findr = findr - currerr/dfr    
		do i=1,2
			x(i) = findr*uvec(i) + center(i)
		end do
		currerr = tokam(x(1),x(2),one,eps,del,kap)
	end do
end function findr

function tokam(x,y,c,eps,del,kap) !tokamak function
	implicit none
	real(kind=8) :: x,y,c,eps,del,kap,d1,d2,d3,tokam

	call switchpars(eps,del,kap,d1,d2,d3)
	tokam = c/8*x**4 + d2*x**2 + d3*(x**4-4*x**2*y**2)+d1
end function tokam

function draw(x,y) !just a weird shape, for testing purposes
	implicit none
	real(kind=8) :: x,y,draw
	draw = x**2 + y**2 - sin(2*y)
endfunction draw

subroutine switchpars(ep,de,kap,d1,d2,d3) !switches from eps,del,kap to d1,d2,d3
	implicit none
	integer::info
	real(kind=8)::ep,de,kap,d1,d2,d3,const,one
	real(kind=8),dimension(3,3)::a,at
	real(kind=8),dimension(3)::b,pvt
	one = 1.0
	const = (-1)*(1.0/8)

	a=reshape((/one,(1+ep)**2,(1+ep)**4,one,(1-ep)**2,(1-ep)**4,one,(1-de*ep)**2,(1-de*ep)**4-4*(1-de*ep)**2*kap**2*ep**2/),shape(a))
	at = transpose(a)
	b = (/const*(1+ep)**4,const*(1-ep)**4,const*(1-de*ep)**4/)

	call dgesv(3,1,at,3,pvt,b,3,info)
	d1 = b(1)
	d2 = b(2)
	d3 = b(3)

end subroutine switchpars

function interp1d(x,x1,x2,y1,y2)
	implicit none
	real(kind=8)::x,x1,x2,y1,y2,interp1d

	interp1d = y1*(1.0d0-(x-x1)/(x2-x1))+y2*((x-x1)/(x2-x1))
end function interp1d
end module curvestuff
