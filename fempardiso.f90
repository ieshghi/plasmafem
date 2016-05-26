PROGRAM fem
	USE mesh
	implicit none
	integer::nx,ny,N,NT,NV,NB,i,j,k,q,r,c,mtype,solver,error,maxfct,mnum,phase,nrhs,msglvl !nx=elements in x, ny=elements in y, N=total number of elements
	integer, dimension(:,:),allocatable::t
	integer, dimension(:),allocatable::b,row1,col1,row2,col2,x,ia,ja,perm,redun
	real, dimension(:,:),allocatable::p !array of points, array of triangle vertices, and big L finite element array
	real, dimension(:),allocatable::f,val1,val2,arr !array of vertices along edge, array of <integral>g_i*f:wq:
	integer, dimension(64):: pt,iparm
	real, dimension(64):: dparm
	real, dimension(3,3)::A !We will use this array in the process of finding equations of planes
	real::det,temp !determinants of matrices, values to insert in sparse matrix
	
	call pardisoinit(pt,mtype,solver,iparm,dparm,error) !Initialise PARDISO
	iparm(3) = 1
	iparm(1)=0
	mtype = 11
	solver=0
	
	write(*,*) 'Nx = ' !How many elements along x?
	read(*,*) nx
	write(*,*) 'Ny = ' !How many elements along y?
	read(*,*) ny
	N = nx*ny !Total number of elements will be obviously nx*ny
	allocate(f(N))
	
	
	call buildmesh(nx,ny,p,t,b) !Builds p,t, and b arrays for later use. 
	NT = size(t(:,1))
	NV = size(p(:,1))
	NB = size(b) 
	allocate(val1(NT*9),col1(NT*9),row1(NT*9),val2(NT*9),col2(NT*9),row2(NT*9),x(N),ia(N+1),perm(N))! Allocate sizes of arrays
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
		do j =1,3 !We build the array with x1,y1,x2,y2,...
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
	
	!Now PARDISO comes in. Col2 is equivalent to the JA array, val2 is the A array. We need to rewrite row2 to be like IA, since it is not currently in compressed sparse row format.

	k=0
	do i = 1,9*NT  !find size of nonzero parts of the row and column arrays, to optimise space usage
		if(col2(i) /= 0) then
				k=k+1
		end if
	end do
	
	allocate(ja(k),arr(k)) !the points which were set to 0 in col2 are a waste of space and irrelevant to the calculation, so we need our final arrays not to include them. So they should be this size
	
	do i = 1,size(row2)-1 !this loop puts row2 into an ordered format, so ia can later be in compressed sparse row (CSR) format.
		do k = i+1,size(row2)-1
				if (row2(i) > row2(k) .and. row2(i)/=0 .and. row2(k) /=0) then
				j = row2(i)
				row2(i) = row2(k)
				row2(k) = j
				j = col2(i)
				col2(i) = col2(k)
				col2(k) = j
				temp = val2(i)
				val2(i) = val2(k)
				val2(k) = temp
			end if
		end do
	end do
	

	ia(1)=1
	j=2
	do i = 2,size(row2) !copies over the information from row2, into the correct CSR format for PARDISO
		if (row2(i)>row2(i-1)) then
			ia(j) = i
			j=j+1
		end if
	end do
	j=1
	ia(N+1) = ia(N)+1
		
	do i = 1,size(col2) !makes arr and ja, the arrays that will be interacting with PARDISO. They mustn't have the redundant spaces that were set to zero in col2
		if(col2(i)/=0) then
			arr(j) = val2(i)
			ja(j) = col2(i)
			j=j+1
		end if
	end do
	do i = 1,size(ia)-1 !makes ja strictly increasing in any given row
		do j = ia(i),ia(i+1)-1
			do q = j,ia(i+1)-1
				if (ja(j)>ja(q))then
						k = ja(j)
						ja(j) = ja(q)
						ja(q) = k
				end if
			end do
		end do
	end do
	phase = 11
	msglvl = 1
	nrhs = 1	
	CALL PARDISO_CHKMATRIX(MTYPE,N,ARR,IA,JA,ERROR)
	CALL PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, N, ARR, IA, JA, PERM, NRHS, IPARM, MSGLVL, f, X, ERROR, DPARM)
		
	
	phase = 22

	CALL PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, N, ARR, IA, JA, PERM, NRHS, IPARM, MSGLVL, f, X, ERROR, DPARM)
	
	iparm(8) = 1
	phase = 33
	
!	CALL PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, N, ARR, IA, JA, PERM, NRHS, IPARM, MSGLVL, f, X, ERROR, DPARM)
	
!	phase = -22
		
!	CALL PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, N, ARR, IA, JA, PERM, NRHS, IPARM, MSGLVL, f, X, ERROR, DPARM)

!	phase = -1
	
	!CALL PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, N, ARR, IA, JA,PERM, NRHS, IPARM, MSGLVL, f, X, ERROR, DPARM)

	END PROGRAM fem
