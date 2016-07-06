MODULE mesh
	implicit none

	CONTAINS
	FUNCTION exact(x,y) !known solution to Poisson equation (for debugging purposes)
		implicit none
		real,parameter::pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986
		real(kind=8),intent(in):: x,y
		real(kind=8) exact

	!	exact = x*(1-x)*y*(1-y)
		
		exact = (-1)*0.5*(1.0/pi)*(1.0/pi)*(sin(pi*x)*sin(pi*y))
		END FUNCTION exact	
	FUNCTION foo(x,y) !right side of the equation
		implicit none
		real,parameter::pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986
		real(kind=8),intent(in):: x,y
		real(kind=8):: foo
		
	!	foo = 2*(x*x-x+y*y-y)
		
		foo = sin(pi*x)*sin(pi*y)
		END FUNCTION foo
			
	FUNCTION  linspace(a,b,n) !equivalent of python linspace
        implicit none
        real, intent(in):: a,b !start and endpoint
        integer(kind=8), intent(in):: n !number of elements
        integer:: i ! loop variable
        real:: dx, linspace(n)
        dx = (b-a)/(n-1) !spacing between x's
        do i = 1,n
                linspace(i)=a+(i-1)*dx ! fill the output array
        end do
		END FUNCTION linspace
	
	FUNCTION threedet(m) !Find the determinant of a 3x3 matrix
		implicit none
		real, dimension(3,3):: m
		real:: threedet

		threedet =  m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + m(1,3)*m(2,1)*m(3,2)
		threedet = threedet - m(1,3)*m(2,2)*m(3,1) - m(1,1)*m(2,3)*m(3,2) - m(1,2)*m(2,1)*m(3,3)
		END FUNCTION threedet	
	
	SUBROUTINE buildmesh (nx,ny,p,t,b) !makes triangle mesh for unit square
	!nx: number of gridpoints in x-direction
	!ny: number of gridpoints in y-direction
	!p: (2,N = vertex number) takes empty allocatable array as input, outputs positions of vertices
	!t: (3,NT = triangle number) takes empty allocatable array as input, outputs triangles
	!b: (1,NE = edge number) array of points on the edge
		implicit none
		integer,dimension(:,:),allocatable::t
		real,dimension(:,:),allocatable::p
		integer,dimension(:),allocatable::b
		integer(kind=8)::nx,ny,i,j !i,j are loop variables
		real,dimension(nx)::x !arrays of x and y positions
		real,dimension(ny)::y
		real::w,z
		!w=1
		!z=2
		x = linspace(0.0,1.0,nx)
		y = linspace(0.0,1.0,ny)

		allocate(p(nx*ny,2),t((nx-1)*(ny-1)*2,3),b((nx+ny)*2-4)) !vertex number = nx*ny, triangle number = (nx-1)*(ny-1)*2	
		
		do i = 1,nx*ny !we first fill in the point coordinates using the values in the x,y arrays we built earlier
			p(i,1)=x(modulo(i-1,nx)+1)
			p(i,2)=y(((i-1)/nx)+1)
		end do
		
		do i = 1,(nx-1)*(ny-1) !making the t array is slightly more bizarre. The first loop fills the first half of the array with the triangles which have a side facing downwards, then the second loop does the other half.
			t(i,1) = i + (i-1)/(nx-1)  
			t(i,2) = i + (i-1)/(nx-1) +1 
			t(i,3) = i + (i-1)/(nx-1) +nx  
			
			j = i + (nx-1)*(ny-1)
			t(j,1) = i + nx + (i-1)/(nx-1)
			t(j,2) = i + nx + (i-1)/(nx-1)+1
			t(j,3) = i + (i-1)/(nx-1)+1
		end do 

		j=1
		do i = 1,(nx*ny) !Now we fill in the b vector. The various if statements select the points along the boundary and assign them to the vector
			if (i<=nx) then
					b(j)=i
					j=j+1
			else if (modulo(i,nx)==1) then
					b(j)=i
					j=j+1
			else if (modulo(i,nx)==0) then
					b(j)=i
					j=j+1
			else if (nx*ny - i <= nx) then
					b(j)=i
					j=j+1
			end if
		end do
				
		END SUBROUTINE buildmesh
	
	SUBROUTINE threeinv(M) !Inverts 3x3 matrix. M is the matrix we want to invert
			real,dimension(3,3)::m
			real::det,a,b,c,d,e,f,g,h,i
			det = threedet(m)

			a = M(2,2)*M(3,3)-M(2,3)*M(3,2)
			b = (-1)*(M(2,1)*M(3,3)-M(2,3)*M(3,1))
			c = M(2,1)*M(3,2)-M(2,2)*M(3,1)
			d = (-1)*(M(1,2)*M(3,3)-M(1,3)*M(3,2))
			e = M(1,1)*M(3,3)-M(1,3)*M(3,1)
			f = (-1)*(M(1,1)*M(3,2)-M(1,2)*M(3,1))
			g = M(1,2)*M(2,3)-M(1,3)*M(2,2)
			h = (-1)*(M(1,1)*M(2,3)-M(1,3)*M(2,1))
			i = M(1,1)*M(2,2)-M(1,2)*M(2,1)

			M(1,1) = a
		    M(1,2) = d
			M(1,3) = g
			M(2,1) = b
			M(2,2) = e
			M(2,3) = h
			M(3,1) = c
			M(3,2) = f
			M(3,3) = i
			M = M*(1/det)
		END SUBROUTINE threeinv
	

	SUBROUTINE distmesh(p,t,b)
                implicit none
                INTEGER, PARAMETER :: maxrecs = 100000
                INTEGER :: J, NR, ios
                CHARACTER(LEN=100) :: inputfile
                CHARACTER(LEN=1) :: junk
                integer,dimension(:,:),allocatable::t
                real(kind=8),dimension(:,:),allocatable::p
                real,dimension(3)::temp
                integer,dimension(:),allocatable::b
                NR = 0
        OPEN(UNIT=1,FILE='infiles/p.txt')
                DO J=1,maxrecs
                        READ(1,*,IOSTAT=ios) junk
                        IF (ios /= 0) EXIT
                        IF (J == maxrecs) THEN
                                STOP
                        ENDIF
                        NR = NR + 1
                ENDDO
                REWIND(1)
                ALLOCATE(p(NR,2))
                DO J=1,NR
                        READ(1,*) p(j,1),p(j,2)
                END DO
                CLOSE(1)
                NR = 0
        OPEN(UNIT=1,FILE='infiles/t.txt')
                DO J=1,maxrecs
                        READ(1,*,IOSTAT=ios) junk
                        IF (ios /= 0) EXIT
                        IF (J == maxrecs) THEN
                                STOP
                        ENDIF
                        NR = NR + 1
                ENDDO
                REWIND(1)
                ALLOCATE(t(NR,3))
                DO J=1,NR
                        READ(1,*) temp(1),temp(2),temp(3)
			t(j,1) = int(temp(1))
                        t(j,2) = int(temp(2))
                        t(j,3) = int(temp(3))
                END DO
                CLOSE(1)
                NR = 0
        OPEN(UNIT=1,FILE='infiles/b.txt')
                DO J=1,maxrecs
                        READ(1,*,IOSTAT=ios) junk
                        IF (ios /= 0) EXIT
                        IF (J == maxrecs) THEN
                                STOP
                        ENDIF
                        NR = NR + 1
                ENDDO
                REWIND(1)
                ALLOCATE(b(NR))
                DO J=1,NR
                        temp(1)=0
                        READ(1,*) temp(1)
                        b(j) = int(temp(1))
                END DO
                CLOSE(1)

        END SUBROUTINE distmesh
	END MODULE mesh
