MODULE mesh
  implicit none

  CONTAINS

  FUNCTION exact(x,y,d1,d2,d3) !known solution to Poisson equation (for debugging purposes)
  implicit none
  real *8::pi
  real *8,intent(in):: x,y,d1,d2,d3
  real *8 exact
  pi = 4*atan(1.0)

  !exact = 1.0/8.0*x**4+d1+d2*x**2+d3*(x**4-4*(x**2)*(y**2))
  exact = (1.0/8.0*x**4+d1+d2*x**2+d3*(x**4-4*(x**2)*(y**2)))/sqrt(x)
  !exact = 1-(x*x+y*y) 
  !exact = x*(1-x)*y*(1-y)
  !exact = (-1)*0.5*(1.0/pi)*(1.0/pi)*(sin(pi*x)*sin(pi*y))
  END FUNCTION exact  
  FUNCTION foo(x,y,d1,d2,d3) !right side of the equation
  implicit none
  real,parameter::pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986
  real *8,intent(in):: x,y,d1,d2,d3
  real *8:: foo
  
  !foo = 3.0/2.0*x**2+2*d2+d3*(12*x**2-4*2*y**2-4*2*x**2)
  foo = (24*d1+(x**2)*(24*d2+(x**2)*(24*d3+35)-96*d3*y**2))/(32*x**(5.0d0/2.0d0))
  !foo = -4.0 
  !foo = 2*(x*x-x+y*y-y)
  !foo = sin(pi*x)*sin(pi*y)
  END FUNCTION foo

  function fx(x,y,d1,d2,d3) !x derivative of rhs
  implicit none
  real *8::x,y,d1,d2,d3,fx

  !fx =3*x+d3*(24*x-4*4*x)
  fx = 3*((x**2)*((-8)*d2+(x**2)*(24*d3+35)+32*d3*y**2)-40*d1)/(64*x**(7.0d0/2.0d0))
  endfunction fx

  function fxx(x,y,d1,d2,d3) !x derivative of rhs
  implicit none
  real *8::x,y,d1,d2,d3,fxx

  !fx =3*x+d3*(24*x-4*4*x)
  fxx = 3*(280*d1 + (x**2)*(24*d2+(x**2)*(24**d3+35)-96*d3*(y**2)))/(128*x**(9.0d0/2.0d0))
  endfunction fxx

  function fxy(x,y,d1,d2,d3) !x derivative of rhs
  implicit none
  real *8::x,y,d1,d2,d3,fxy

  fxy = 3*y*d3/(x**(3.0d0/2.0d0))
  endfunction fxy

  function fy(x,y,d1,d2,d3) !y derivative of rhs
  implicit none
  real *8::x,y,d1,d2,d3,fy

  !fy =(-1)*d3*4*4*y
  fy = (-6)*d3*y/sqrt(x)
  endfunction fy

  function exactx(x,y,d1,d2,d3) !x derivative of solution 
  implicit none
  real *8::x,y,d1,d2,d3,exactx
  
!  exactx = 1.0/2.0*x**3+d2*2*x+d3*(4*x**3-8*x*y**2)
  exactx = ((x**2)*(24*d2+7*(8*d3+1)*(x**2)-96*d3*y**2)-8*d1)/(16*x**(3.0d0/2.0d0))
  endfunction exactx

  function exacty(x,y,d1,d2,d3) !y derivative of solution
  implicit none
  real *8::x,y,d1,d2,d3,exacty

!  exacty = (-1)*d3*8*y*x**2
  exacty =  (-8)*d3*y*x**(3.0d0/2.0d0)
  endfunction exacty

  function exactxx(x,y,d1,d2,d3) !x derivative of solution 
  implicit none
  real *8::x,y,d1,d2,d3,exactxx
  
!  exactx = 1.0/2.0*x**3+d2*2*x+d3*(4*x**3-8*x*y**2)
  exactxx = (24*d1+(x**2)*(24*d2+(x**2)*35*(8*d3+1)-96*d3*(y**2)))/(32*(x**(5.0d0/2.0d0)))
  endfunction exactxx

  function exactxy(x,y,d1,d2,d3) !x derivative of solution 
  implicit none
  real *8::x,y,d1,d2,d3,exactxy
  
!  exactx = 1.0/2.0*x**3+d2*2*x+d3*(4*x**3-8*x*y**2)
  exactxy = (-12)*d3*sqrt(x)*y
  endfunction exactxy

  function exactyy(x,y,d1,d2,d3) !x derivative of solution 
  implicit none
  real *8::x,y,d1,d2,d3,exactyy
  
!  exactx = 1.0/2.0*x**3+d2*2*x+d3*(4*x**3-8*x*y**2)
  exactyy = (-8)*d3*x**(3.0d0/2.0d0)
  endfunction exactyy

  SUBROUTINE firstder (d1,d2,d3,x,p,t,b,u,xory) !solves given the boundary u 
  implicit none
  integer::xory !x or y derivative? x=0,y=1
  integer *8::N,NT,NV,NB,i,j,k,q !nx=elements in x, ny=elements in y, N=total number of elements
  integer, dimension(:,:)::t
  integer, dimension(:)::b
  integer, dimension(:),allocatable::row1,col1,row2,col2,ia,ja
  real *8, dimension(:,:)::p!array of points, array of triangle vertices, and big L finite element array
  real *8, dimension(:),allocatable::fu,val1,val2,arr,x,err !array of vertices along edge, array of <integral>g_i*f
  real *8, dimension(3,3)::A !We will use this array in the process of finding equations of planes
  real *8::det,temp !determinants of matrices, values to insert in sparse matrix
  real *8::d1,d2,d3
  real *8,dimension(:)::u
    
  NT = size(t(:,1))
    NV = size(p(:,1))
    NB = size(b)
    N = NV !Total number of elements will be the number of vertices
  allocate(fu(N))
    allocate(val1(NT*9),col1(NT*9),row1(NT*9),val2(NT*9),col2(NT*9),row2(NT*9),x(N))! Allocate sizes of arrays
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
  
  if(xory==0)then
          temp = fx((p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0),(p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0),d1,d2,d3) !2d midpoint rule
  elseif(xory==1)then
          temp = fy((p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0),(p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0),d1,d2,d3) !2d midpoint rule
  elseif(xory==2)then
          temp = fxx((p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0),(p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0),d1,d2,d3) !2d midpoint rule
  elseif(xory==3)then
          temp = fxy((p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0),(p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0),d1,d2,d3) !2d midpoint rule
  else
    write(*,*) "Wrong X/Y selection value"
  endif

        fu(t(i,1)) = fu(t(i,1)) + det*0.5*temp/3.0 !Here, we add the result of the convolution of the basis function with the right hand side of the Poisson equation, which gives us the right hand side vector in the finite element equation.
        fu(t(i,2)) = fu(t(i,2)) + det*0.5*temp/3.0
        fu(t(i,3)) = fu(t(i,3)) + det*0.5*temp/3.0
  
    end do

  do i = 1,size(col2)
  col2(i)=0
  row2(i)=0
  val2(i)=0
  enddo
    
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
    do i=1,NB !This loops through the b array and sets the corresponding row of L to all zeros except for the L(b(i),b(i)) spot. It also sets the f(b(i) cell to zero. This allows for correct evaluation of the edges.
        fu(b(i)) = u(i)
        do j = 1,NT*9
            if (row2(j) == b(i)) then
                if (col2(j) == b(i)) then
                    val2(j) = 1.0
                else
                    val2(j) = 0.0
                end if
            end if
        end do
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
    temp = 1e-16
    i = NB/4
    q = i
    call mgmres_st(N,size(ia),ia,ja,arr,x,fu,i,q,temp,temp)
  endsubroutine firstder

  SUBROUTINE poissolve(d1,d2,d3,src_loc,x,src_val,areas,bound,order) !src_val is used in other codes, it is an array of f(x)*triangle_are @ centroids of triangles
  implicit none
  integer::order
  integer *8::N,NT,NV,NB,i,j,k,q,nedge !nx=elements in x, ny=elements in y, N=total number of elements
  integer, dimension(3)::tri
  integer, dimension(:,:),allocatable::t
  integer, dimension(:),allocatable::b,row1,col1,row2,col2,ia,ja
  real *8, dimension(:,:),allocatable::p,src_loc !array of points, array of triangle vertices, and big L finite element array
  real *8, dimension(:),allocatable::fu,val1,val2,arr,src_val,x,areas !array of vertices along edge, array of <integral>g_i*f
  real *8, dimension(:)::bound
  real *8, dimension(3,3)::A !We will use this array in the process of finding equations of planes
  real *8::d1,d2,d3,eps,del,kap,det,temp !determinants of matrices, values to insert in sparse matrix
  real *8, dimension(2)::centroid, pt1,pt2,pt3
  logical, dimension(:),allocatable::edgetri
  call distmesh(p,t,b,eps,del,kap) !Builds p,t, and b arrays for later use. 
  NT = size(t(:,1))
  NV = size(p(:,1))
  NB = size(b)
  N = NV !Total number of elements will be the number of vertices
  allocate(fu(N))
  allocate(val1(NT*9),col1(NT*9),row1(NT*9),val2(NT*9),col2(NT*9),row2(NT*9),x(N),edgetri(NT))! Allocate sizes of arrays
  q = 1
  do i = 1,NT !Loop to build the row and col vectors
    do j = 0,8
      row1(q+j) = t(i,j/3+1)
      col1(q+j) = t(i,modulo(j,3)+1)
    end do
    q = q+9
  end do

  nedge = 0
  do i=1,nt
    tri = t(i,:)
    k=0
    do j=1,nb
      if(t(i,1)==b(j) .or. t(i,2)==b(j) .or. t(i,3)==b(j))then
        k = 1
      endif
    enddo
    if(k==1)then
      edgetri(i)=.true.
      nedge = nedge+1 
    else
      edgetri(i)=.false.
    endif
  enddo

  allocate(src_val(NT+(2)*nedge),src_loc(2,NT+(2)*nedge),areas(NT+2*nedge))
  do i = 1,N
    fu(i)=0
    src_val(i)=0
  end do

  q = 1
  k = 1
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
    if (order==1)then 
        temp = foo((p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0),(p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0),d1,d2,d3) !2d midpoint rule
    elseif(order==2)then
        temp = fx((p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0),(p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0),d1,d2,d3) !2d midpoint rule
    endif
    fu(t(i,1)) = fu(t(i,1)) + det*0.5*temp/3.0 !Here, we add the result of the convolution of the basis function with the right hand side of the Poisson equation, which gives us the right hand side vector in the finite element equation.
    fu(t(i,2)) = fu(t(i,2)) + det*0.5*temp/3.0
    fu(t(i,3)) = fu(t(i,3)) + det*0.5*temp/3.0
  
    if(edgetri(i))then
      centroid(1) = (p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0)
      centroid(2) = (p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0)
      pt1(1) = (p(t(i,1),1)+p(t(i,2),1)+centroid(1))/3
      pt1(2) = (p(t(i,1),2)+p(t(i,2),2)+centroid(2))/3
      pt2(1) = (p(t(i,2),1)+p(t(i,3),1)+centroid(1))/3
      pt2(2) = (p(t(i,2),2)+p(t(i,3),2)+centroid(2))/3
      pt3(1) = (p(t(i,1),1)+p(t(i,3),1)+centroid(1))/3
      pt3(2) = (p(t(i,1),2)+p(t(i,3),2)+centroid(2))/3
      
      src_loc(1,k) = pt1(1)
      src_loc(2,k) = pt1(2)
      if(order==1)then
          src_val(k) = foo(pt1(1),pt1(2),d1,d2,d3)*det*0.5/3
      elseif(order==2)then
          src_val(k) = fx(pt1(1),pt1(2),d1,d2,d3)*det*0.5/3
      endif
      areas(k) = det*0.5/3
      src_loc(1,k+1) = pt2(1)
      src_loc(2,k+1) = pt2(2)
      if(order==1)then
        src_val(k+1) = foo(pt2(1),pt2(2),d1,d2,d3)*det*0.5/3
      elseif(order==2)then
        src_val(k+1) = fx(pt2(1),pt2(2),d1,d2,d3)*det*0.5/3
      endif
      areas(k+1) = det*0.5/3
      src_loc(1,k+2) = pt3(1)
      src_loc(2,k+2) = pt3(2)
      if(order==1)then
        src_val(k+2) = foo(pt3(1),pt3(2),d1,d2,d3)*det*0.5/3
      elseif(order==2)then
        src_val(k+2) = fx(pt3(1),pt3(2),d1,d2,d3)*det*0.5/3
      endif
      areas(k+2) = det*0.5/3
      k = k+3    
    else
      src_loc(1,k) = (p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0)
      src_loc(2,k) = (p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0)
      src_val(k) = temp*det*0.5
      areas(k) = det*0.5
      k = k+1
    endif
  end do
  do i = 1,size(col2)
    col2(i)=0
    row2(i)=0
    val2(i)=0
  enddo

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
  do i=1,NB !This loops through the b array and sets the corresponding row of L to all zeros except for the L(b(i),b(i)) spot. It also sets the f(b(i) cell to zero. This allows for correct evaluation of the edges.
    fu(b(i)) = bound(i)
    do j = 1,NT*9
      if (row2(j) == b(i)) then
        if (col2(j) == b(i)) then
          val2(j) = 1.0
        else
          val2(j) = 0.0
        end if
      end if
    end do
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
  temp = 1e-14
  i = NB/4
  q = i
  call mgmres_st(N,size(ia),ia,ja,arr,x,fu,i,q,temp,temp)
  END SUBROUTINE poissolve
 
  FUNCTION  linspace(a,b,n) !equivalent of python linspace
    implicit none
    real *8, intent(in):: a,b !start and endpoint
    integer, intent(in):: n !number of elements
    integer:: i ! loop variable
    real *8:: dx, linspace(n)
    dx = (b-a)/(n-1) !spacing between x's
    do i = 1,n
        linspace(i)=a+(i-1)*dx ! fill the output array
    end do
  END FUNCTION linspace

  FUNCTION  linspace2(a,b,n) !equivalent of python linspace
    implicit none
    real *8, intent(in):: a,b !start and endpoint
    integer, intent(in):: n !number of elements
    integer:: i ! loop variable
    real *8:: dx, linspace2(n)
    dx = (b-a)/(n) !spacing between x's
    do i = 1,n
        linspace2(i)=a+(i-1)*dx ! fill the output array
    end do
  endfunction linspace2


  FUNCTION threedet(m) !Find the determinant of a 3x3 matrix
  implicit none
  real *8, dimension(3,3):: m
  real *8:: threedet

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
  real *8,dimension(:,:),allocatable::p
  integer,dimension(:),allocatable::b
  integer::nx,ny,i,j !i,j are loop variables
  real *8,dimension(nx)::x !arrays of x and y positions
  real *8,dimension(ny)::y
  real *8::w,z
  !w=1
  !z=2
  x = linspace(0.0d0,1.0d0,nx)
  y = linspace(0.0d0,1.0d0,ny)

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
    real *8,dimension(3,3)::m
    real *8::det,a,b,c,d,e,f,g,h,i
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
  

  SUBROUTINE distmesh(p,t,b,eps,del,kap)
        implicit none
        INTEGER *8, PARAMETER :: maxrecs = 1000000
        INTEGER *8 :: J, NR, ios
        CHARACTER(LEN=100) :: inputfile
        CHARACTER(LEN=1) :: junk
        integer,dimension(:,:),allocatable::t
        real *8,dimension(:,:),allocatable::p
        real *8,dimension(3)::temp
        integer,dimension(:),allocatable::b
  real *8::eps,del,kap
        NR = 0
  open(unit=1,file='infiles/params.txt')
  read(1,*) eps, del, kap
  close(1)

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
