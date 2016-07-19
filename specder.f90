PROGRAM specder
	USE mesh
	USE, INTRINSIC :: iso_c_binding
	IMPLICIT NONE
	INCLUDE '/usr/include/fftw3.f03'
	type(C_PTR) :: plan
	integer::i
	integer, parameter:: N=1000
	integer(kind=8), parameter::M=1000
	complex(C_DOUBLE_COMPLEX), dimension(N) :: input,output,input2,output2
	real,dimension(N)::x,k,der
	real,parameter::pi =3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862
	real::xmin,xmax,dx
	xmax = 6*pi
	xmin = 0.0
	dx = (xmax-xmin)/(N)
	do i = 1,N
		x(i) = xmin + dx*(i-1)
	end do
	do i = 1,N
		input(i) = exp(sin(3.5*pi*x(i)))
		der(i) = 3.5*pi*cos(3.5*pi*x(i))*exp(sin(3.5*pi*x(i)))
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
	
	open(1,file='./files/deriv.dat') !write stuff to file before memory release, so that it can be plotted in python
        	do i = 1,N
                	write(1,*) real(input(i)),der(i),real(output2(i)),imag(output2(i)),x(i)
                end do
        close(1)

	
END PROGRAM specder
