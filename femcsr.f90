PROGRAM fem
	USE mesh
	implicit none
	integer::nx,ny,N,NT,NV,NB,i,j,k,q,r,c,mtype,solver,error,maxfct,mnum,phase,nrhs,msglvl !nx=elements in x, ny=elements in y, N=total number of elements
	integer, dimension(:,:),allocatable::t
	integer, dimension(:),allocatable::b,row1,col1,row2,col2,x,ia,ja,perm,redun
	real, dimension(:,:),allocatable::p !array of points, array of triangle vertices, and big L finite element array
	real, dimension(:),allocatable::f,fu,val1,val2,arr !array of vertices along edge, array of <integral>g_i*f:wq:
	integer, dimension(64):: pt,iparm,dparm
	real, dimension(3,3)::A !We will use this array in the process of finding equations of planes
	real::det,temp !determinants of matrices, values to insert in sparse matrix
	
	mtype = 11
	solver=10
	call pardisoinit(pt,mtype,solver,iparm,dparm,error) !Initialise PARDISO
	iparm(3) = 1
	
	write(*,*) 'Nx = ' !How many elements along x?
	read(*,*) nx
	write(*,*) 'Ny = ' !How many elements along y?
	read(*,*) ny
	N = nx*ny !Total number of elements will be obviously nx*ny
	allocate(f(N),fu(N))
	
	
	call buildmesh(nx,ny,p,t,b) !Builds p,t, and b arrays for later use. 
	NT = size(t(:,1))
	NV = size(p(:,1))
	NB = size(b) 
	allocate(val1(NT*9),col1(NT*9),row1(NT*9),val2(NT*9),col2(NT*9),row2(NT*9),x(N),perm(N))! Allocate sizes of arrays
	q = 1
	do i = 1,NT !Loop to build the row and col vectors
		do j = 0,8
			row1(q+j) = t(i,j/3+1)
			col1(q+j) = t(i,modulo(j,3)+1)
		end do
		q = q+9
	end do
	do i = 1,N
		f(i)=0
	end do	
	q = 1
	do i = 1,NT !loop through triangles to add entries into the vals vector, and build the f vector
		do j =1,3 !We build the matrix with x1,y1,x2,y2,...
			A(j,1) = (p(t(i,j),1))
			A(j,2) = (p(t(i,j),2))
			A(j,3) = 1.0
		end do
		
		det = abs(threedet(A)) !This computes the determinant of that matrix, which is used to compute the area of the triangle in question
		call threeinv(A)
		do j=0,8
			val1(q+j) = (A(1,j/3+1)*A(1,modulo(j,3)+1)+A(2,j/3+1)*A(2,modulo(j,3)+1))*det*0.5
		end do
		q = q+9

		f(t(i,1)) = f(t(i,1)) + det*0.5*foo(p(t(i,1),1),p(t(i,1),2))*0.5 !Here, we add the result of the convolution of the basis function with the right hand side of the Poisson equation, which gives us the right hand side vector in the finite element equation.
		f(t(i,2)) = f(t(i,2)) + det*0.5*foo(p(t(i,2),1),p(t(i,2),2))*0.5
		f(t(i,3)) = f(t(i,3)) + det*0.5*foo(p(t(i,3),1),p(t(i,3),2))*0.5
	end do

	do i=1,NB !This loops through the b array and sets the corresponding row of L to all zeros except for the L(b(i),b(i)) spot. It also sets the f(b(i) cell to zero. This allows for correct evaluation of the edges.
		f(b(i)) = 0
		do j = 1,NT*9
			if (row1(j) == b(i)) then
				if (col1(j) == b(i)) then
					val1(j) = 1
				else
					val1(j) = 0
				end if
			end if
		end do
	end do
	
	j=1
	do i=1,size(col1) !Reprocess information in row1 and col1, to compress them from arrays of size 9NT to arrays of size N.
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
	
	
	!Now PARDISO comes in. Col2 is equivalent to the JA array, val2 is the A array. We need to rewrite row2 to be like IA, since it is not currently in compressed sparse row format. (or maybe not, this version does no such rearrangement. For extra compression, go to fempardiso.f90)

	k=0
	do i = 1,9*NT  !find size of nonzero parts of the row and column arrays, to optimise space usage
		if(col2(i) /= 0) then
				k=k+1
		end if
	end do
	
	allocate(ja(k),arr(k),ia(k)) !this creates arrays of the right size
	
	j=1
	do i = 1,size(col2)
		if(col2(i)/=0) then
			arr(j) = val2(i)
			ja(j) = col2(i)
			ia(j) = row2(i)
			j=j+1
		end if
	end do

	open(1,file='ia.dat',status='new') 
		  do i = 1,size(ia)
			write(1,*) ia(i)
  			end do 
		  close(1)
		  
	open(2,file='ja.dat',status='new') 
		  do i = 1,size(ja)
			write(2,*) ja(i)
		  end do
		  close(2)
	
	open(3,file='fu.dat',status='new') 
		  do i = 1,N
			write(3,*) f(i)
		  end do
		  close(3)

	open(4,file='arr.dat',status='new') 
		  do i = 1,size(arr)
			write(4,*) arr(i)
		  end do
		  close(4)
	open(3,file='p.dat',status='new') 
		  do i = 1,NV
			write(3,*) p(i,1),p(i,2)
		  end do
		  close(5)
	open(6,file='t.dat',status='new') 
		  do i = 1,NT
			write(6,*) t(i,1),t(i,2),t(i,3)
		  end do
		  close(6)
	
	END PROGRAM fem
