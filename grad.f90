PROGRAM fem
	USE meshgs
	implicit none
	integer(kind=8)::nx,ny,N,NT,NV,NB,i,j,k,q !nx=elements in x, ny=elements in y, N=total number of elements
	integer, dimension(:,:),allocatable::t
	integer, dimension(:),allocatable::b,row1,col1,row2,col2,ia,ja
	real(kind=8), dimension(:,:),allocatable::p !array of points, array of triangle vertices, and big L finite element array
	real(kind=8), dimension(:),allocatable::fu,val1,val2,x,x_exact,arr !array of vertices along edge, array of <integral>g_i*f
	real(kind=8), dimension(3,3)::A !We will use this array in the process of finding equations of planes
	real(kind=8)::det,temp,centx,centy !determinants of matrices, values to insert in sparse matrix
	
!	write(*,*) 'Nx = ' !How many elements along x?
!	read(*,*) nx
!!	write(*,*) 'Ny = ' !How many elements along y?
!	read(*,*) ny
!	N = nx*ny
!	allocate(fu(N))	
!	call buildmesh(nx,ny,p,t,b) !Builds p,t, and b arrays for later use. 

	call distmesh(p,t,b)
	
	NT = size(t(:,1))
	NV = size(p(:,1))
	NB = size(b)

	N = NV !Total number of elements will be the number of vertices
	allocate(fu(N))
	allocate(val1(NT*9),col1(NT*9),row1(NT*9),val2(NT*9),col2(NT*9),row2(NT*9),x(N),x_exact(N))! Allocate sizes of arrays
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
	end do
	q = 1
	do i = 1,NT !loop through triangles to add entries into the vals vector, thus building the left hand side matrix
		do j =1,3 !We build the array with x1,y1,x2,y2,...
			A(j,1) = (p(t(i,j),1))
			A(j,2) = (p(t(i,j),2))
			A(j,3) = 1.0
		end do
			
		det = abs(threedet(A)) !This computes the determinant of that matrix, which is used to compute the area of the triangle in question
		call threeinv(A)
		centx = (p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0) !x-coord of centroid of triangle
		centy = (p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0) !y-coord of centroid of triangle
		do j=0,8
			val1(q+j) = (-1)*(A(1,j/3+1)*A(1,modulo(j,3)+1)+A(2,j/3+1)*A(2,modulo(j,3)+1))*det*0.5/(centx) !evaluating the integral of the basis functions, multiplied by 1/x (grad-schafranov)
		end do
		q = q+9
!	end do

!	do i = 1,NT !build right hand side, looping over all triangles
		centx = (p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0) !x-coord of centroid of triangle
		centy = (p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0) !y-coord of centroid of triangle
		temp = foo(centx,centy)/(centx) !2d midpoint rule, multiplied by 1/x (grad-schafranov)
		
		fu(t(i,1)) = fu(t(i,1)) + det*0.5*temp/3.0 !Here, we add the result of the convolution of the basis function with the right hand side of the Poisson equation, which gives us the right hand side vector in the finite element equation.
		fu(t(i,2)) = fu(t(i,2)) + det*0.5*temp/3.0
		fu(t(i,3)) = fu(t(i,3)) + det*0.5*temp/3.0
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
	
	do i = 1,N
		x_exact(i) = exact(p(i,1),p(i,2))
	end do
	
	open(1,file='./files/conv.dat',access='append')
		write(1,*) maxval(abs(x-x_exact)),sqrt(float(NT))
	close(1)
	
	open(1,file='./files/p.dat') !write stuff to file before memory release, so that it can be plotted in python
		  do i = 1,NV
			write(1,*) p(i,1),p(i,2)
		  end do
		  close(1)
	open(2,file='./files/t.dat') 
		  do i = 1,NT
			write(2,*) t(i,1),t(i,2),t(i,3)
		  end do
		  close(2)
	open(3,file='./files/x.dat') 
		  do i = 1,N
			write(3,*) x(i)
			end do
		  close(3)
	open(3,file='./files/exact.dat') 
		  do i = 1,N
			write(3,*) x_exact(i)
			end do
		  close(3)
	open(3,file='./files/fu.dat') 
		  do i = 1,N
			write(3,*) fu(i)
			end do
		  close(3)
	open(3,file='./files/ia.dat') 
		  do i = 1,size(ia)
				write(3,*) ia(i)
		  end do
		  close(3)
	open(3,file='./files/ja.dat') 
		  do i = 1,size(ja)
				write(3,*) ja(i)
		  end do
		  close(3)
	open(3,file='./files/arr.dat') 
		  do i = 1,size(arr)
				write(3,*) arr(i)
		  end do
		  close(3)

	END PROGRAM fem
