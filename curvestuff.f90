module curvestuff
contains

subroutine arcparam(a,b,tarc,darc,N,L,eps,del,kap)
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
		do while(abs(currerr) > 1e-8)
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


subroutine lgmap(x,w,a,b,mode)
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


function dxdy(theta,eps,del,kap)
	implicit none
	real(kind=8),parameter::infi = 1e-10
	real(kind=8) ::theta,eps,del,kap
	real(kind=8), dimension(2)::dxdy
	real(kind=8),dimension(7)::args
	
	args = (/eps,del,kap,real(0.7,kind=8),real(infi,kind=8),real(1.0,kind=8),real(0.0,kind=8)/)

	dxdy(1) = (findr_fast(theta+infi,args)*cos(theta+infi)-findr_fast(theta-infi,args)*cos(theta-infi))/(2*infi)
	dxdy(2) = (findr_fast(theta+infi,args)*sin(theta+infi)-findr_fast(theta-infi,args)*sin(theta-infi))/(2*infi)
		

end function dxdy

function findr_fast(theta,args)
	implicit none
	real(kind=8),dimension(7)::args
	real(kind=8)::eps,del,kap,guess,theta,errtol,centx,centy
	real::findr_fast

	eps = args(1)
	del = args(2)
	kap = args(3)
	guess = args(4)
	errtol = args(5)
	centx = args(6)
	centy = args(7)

	findr_fast = findr(eps,del,kap,theta,guess,errtol,centx,centy)


end function findr_fast

function findr(eps,del,kap,theta,guess,errtol,centx,centy)
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

function tokam(x,y,c,eps,del,kap)
	implicit none
	real(kind=8) :: x,y,c,eps,del,kap,d1,d2,d3,tokam

	call switchpars(eps,del,kap,d1,d2,d3)
	tokam = c/8*x**4 + d2*x**2 + d3*(x**4-4*x**2*y**2)+d1
end function tokam

subroutine switchpars(ep,de,kap,d1,d2,d3)
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

end module curvestuff
