program test
	implicit none
	real(kind=8),dimension(10,2)::srcloc
	real(kind=8),dimension(10)::srcval
	real(kind=8),dimension(5,2)::targloc
	real(kind=8),dimension(5,2)::targnorm
	complex(kind=16),dimension(5)::pot
	integer:: i,n,m,ntj

	n=10
	m=5
	ntj=3

	do i=1,n
		srcloc(i,1)=i
		srcloc(i,2)=2*i
		srcval(i)=10
	enddo

	do i=1,m
		targloc(i,1)=2.0/i
		targloc(i,2)=3*i
		targnorm(i,1)=4*i
		targnorm(i,2)=5
	enddo

	call l2dacquadwrap(srcloc,srcval,targloc,targnorm,n,m,ntj,-1,pot)

	write(*,*) pot
	


end program test
