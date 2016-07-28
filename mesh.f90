MODULE mesh
	implicit none

	CONTAINS
	
	SUBROUTINE poissolve(src_loc,x,src_val) !src_val is used in other codes, it is an array of f(x)*triangle_are @ centroids of triangles
	implicit none
        integer(kind=8)::N,NT,NV,NB,i,j,k,q !nx=elements in x, ny=elements in y, N=total number of elements
        integer, dimension(:,:),allocatable::t
        integer, dimension(:),allocatable::b,row1,col1,row2,col2,ia,ja
        real(kind=8), dimension(:,:),allocatable::p,src_loc !array of points, array of triangle vertices, and big L finite element array
        real(kind=8), dimension(:),allocatable::fu,val1,val2,arr,src_val,x !array of vertices along edge, array of <integral>g_i*f
        real(kind=8), dimension(3,3)::A !We will use this array in the process of finding equations of planes
        real::det,temp !determinants of matrices, values to insert in sparse matrix
        call distmesh(p,t,b) !Builds p,t, and b arrays for later use. 
        NT = size(t(:,1))
        NV = size(p(:,1))
        NB = size(b)
        N = NV !Total number of elements will be the number of vertices
        allocate(fu(N))
        allocate(val1(NT*9),col1(NT*9),row1(NT*9),val2(NT*9),col2(NT*9),row2(NT*9),x(N),src_val(NT),src_loc(NT,2))! Allocate sizes of arrays
        q = 1
        do i = 1,NT !Loop to build the row and col vectors
                do j = 0,8
                        row1(q+j) = t(i,j/3+1)
                        col1(q+j) = t(i,modulo(j,3)+1)
                end do
                q = q+9
        end do

        do i = 1,N
                fu(i)=0
		src_val=0
        end do

        q = 1
        do i = 1,NT !loop through triangles to add entries into the vals vector, and build the f vector
                do j =1,3 !We build the array with x1,y1,x2,y2,...
                        A(j,1) = (p(t(i,j),1))
                        A(j,2) = (p(t(i,j),2))
                        A(j,3) = 1.0
                end do

                det = abs(threedet(A)) !This computes the determinant of that matrix, which is used to compute the area of the triangle in question
                call threeinv(A)


                do j=0,8
                        val1(q+j) = (-1)*(A(1,j/3+1)*A(1,modulo(j,3)+1)+A(2,j/3+1)*A(2,modulo(j,3)+1))*det*0.5
                end do
		q = q+9

                temp = foo((p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0),(p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0)) !2d midpoint rule

                fu(t(i,1)) = fu(t(i,1)) + det*0.5*temp/3.0 !Here, we add the result of the convolution of the basis function with the right hand side of the Poisson equation, which gives us the right hand side vector in the finite element equation.
                fu(t(i,2)) = fu(t(i,2)) + det*0.5*temp/3.0
                fu(t(i,3)) = fu(t(i,3)) + det*0.5*temp/3.0
		
		src_loc(i,1) = p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1)/(3.0)
		src_loc(i,2) = p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2)/(3.0)
		src_val(i) = temp*det*0.5

        end do

        do i=1,NB !This loops through the b array and sets the corresponding row of L to all zeros except for the L(b(i),b(i)) spot. It also sets the f(b(i) cell to zero. This allows for correct evaluation of the edges.
                fu(b(i)) = 0
                do j = 1,NT*9
                        if (row1(j) == b(i)) then
                                if (col1(j) == b(i)) then
                                        val1(j) = 1.0
                                else
                                        val1(j) = 0.0
                                end if
                        end if
                end do
        end do

        j=1
        do i=1,size(col1) !Reprocess information in row1 and col1, to compress them. Redundancies are set to zero, to be ignored later.
                if (col1(i) /= 0 .and. row1(i) /= 0) then
                        col2(j) = col1(i)
                        row2(j) = row1(i)
                        val2(j) = val1(i)

                        do k=i+1,size(col1)
                                if (col1(k)==col2(j) .and. row1(k)==row2(j)) then
                                        col1(k)=0
                                        row1(k)=0
                                        val2(j) = val2(j) + val1(k)
                                end if
                        end do
                        j = j+1
                end if
        end do

        k=0
        do i = 1,size(col2)  !find size of nonzero parts of the row and column arrays
                if(col2(i) /= 0) then
			k=k+1
                end if
        end do

        allocate(ia(k),ja(k),arr(k)) !the points which were set to 0 in col2 are a waste of space and irrelevant to the calculation, so we need our final arrays not to include them. So they should be this size

        j=1
        do i = 1,size(col2) !makes arr and ja, the value and column arrays that will be interacting with PARDISO. They mustn't have the redundant spaces that were set to zero in col2
                if(col2(i)/=0) then
			arr(j) = val2(i)
                        ja(j) = col2(i)
                        ia(j) = row2(i)
                        j=j+1
                end if
        end do
        do i = 1,N !We want X to have an estimate of the solution to begin with... I start with 0, for lack of a better idea
                x(i) = 0
        end do
        temp = 1e-6
        i = NB/4
        q = i
        call mgmres_st(N,size(ia),ia,ja,arr,x,fu,i,q,temp,temp)

	END SUBROUTINE poissolve



	FUNCTION exact(x,y) !known solution to Poisson equation (for debugging purposes)
		implicit none
		real,parameter::pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986
		real(kind=8),intent(in):: x,y
		real(kind=8) exact

		exact = x*x+y*y-1 !zero on the unit circle
		!exact = x*(1-x)*y*(1-y)
		!exact = (-1)*0.5*(1.0/pi)*(1.0/pi)*(sin(pi*x)*sin(pi*y))
		END FUNCTION exact	
	FUNCTION foo(x,y) !right side of the equation
		implicit none
		real,parameter::pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986
		real(kind=8),intent(in):: x,y
		real(kind=8):: foo
		
		foo = 4.0 
		!foo = 2*(x*x-x+y*y-y)
		!foo = sin(pi*x)*sin(pi*y)
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
		real(kind=8), dimension(3,3):: m
		real(kind=8):: threedet

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
			real(kind=8),dimension(3,3)::m
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
