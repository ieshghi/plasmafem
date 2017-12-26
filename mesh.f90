module mesh
  use functions
  implicit none
  contains
  
  subroutine weirdder (c,x,p,t,b,u,xory,uat,uxat) !Solves grad-shafranov with no oversampling, and given a boundary u. Also given
!          solutions uat (solution to first equation), and uxat (can contain first derivative in x or y)
  use functions
  implicit none
  integer::xory !Which solve are we doing? 
  integer *8::n,nt,nv,nb,i,j,k,q 
  integer, dimension(:,:)::t
  integer, dimension(:)::b
  integer, dimension(:),allocatable::row1,col1,row2,col2,ia,ja
  real *8, dimension(:),optional::uxat
  real *8, dimension(:,:)::p
  real *8, dimension(:),allocatable::fu,val1,val2,arr,x,err 
  real *8, dimension(3,3)::a
  real *8::det,temp
  real *8::cx,cy,c
  real *8,dimension(:)::u,uat

  nt = size(t(:,1))
  nv = size(p(:,1))
  nb = size(b)
  n = nv !total number of elements will be the number of vertices
  allocate(fu(n))
  allocate(val1(nt*9),col1(nt*9),row1(nt*9),val2(nt*9),col2(nt*9),row2(nt*9),x(n))! allocate sizes of arrays
  q = 1
  do i = 1,nt !loop to build the row and col vectors
    do j = 0,8
      row1(q+j) = t(i,j/3+1)
      col1(q+j) = t(i,modulo(j,3)+1)
    end do
  q = q+9
  end do
  do i = 1,n
    fu(i)=0 !initialise right hand side vector to 0
  end do
  q = 1
  do i = 1,nt !loop through triangles to add entries into the stiffness matrix and right hand side
    do j =1,3 !Trick to get the integral of dot product of two plane basis functions. First build this array A
      a(j,1) = (p(t(i,j),1))
      a(j,2) = (p(t(i,j),2))
      a(j,3) = 1.0
    end do

    det = abs(threedet(a)) !this computes the determinant of that matrix, which is used to compute the area of the triangle in question
    call threeinv(a) !We invert the array
    do j=0,8
      val1(q+j) = (-1)*(a(1,j/3+1)*a(1,modulo(j,3)+1)+a(2,j/3+1)*a(2,modulo(j,3)+1))*det*0.5 !This somehow just gives the answer...
      val1(q+j) = val1(q+j) + (-1)*1./9*stiff((p(t(i,1),1)+p(t(i,2),1)+&
        p(t(i,3),1))/3,(p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/3,c)*det*0.5 !And this gives us the correction to the stiffness matrix,
!which we evaluate with a coarse midpoint rule. 
!for the "weird poisson equation"
    end do
    q = q+9 !We have a bizarre looping structure. Each triangle contributes to nine spots in the matrix
    cx = (p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0) !coordinates of the centroid of the triangle in question
    cy = (p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0)
    if(xory==0)then !xory can take four values, each of which corresponds to a different derivative. 
      temp = gsrhsx((uat(t(i,1))+uat(t(i,2))+uat(t(i,3)))/3,cx,cy,c) !2d midpoint rule for the right hand side. 
    elseif(xory==1)then
      temp = gsrhsy((uat(t(i,1))+uat(t(i,2))+uat(t(i,3)))/3,cx,cy,c)
    elseif(xory==2)then
      temp = gsrhsxx((uat(t(i,1))+uat(t(i,2)) + uat(t(i,3)))/3,(uxat(t(i,1))+uxat(t(i,2))+uxat(t(i,3)))/3,cx,cy,c)
    elseif(xory==3)then
      temp = gsrhsxy((uat(t(i,1))+uat(t(i,2)) + uat(t(i,3)))/3,(uxat(t(i,1))+uxat(t(i,2))+uxat(t(i,3)))/3,cx,cy,c)!uxat seems like it
      !      contains an x-derivative of the function, but in this case we are feeding in a y-derivative. I didn't want to declare
      !      two redundant arrays.
    else
      write(*,*) "wrong x/y selection value"
    endif

    fu(t(i,1)) = fu(t(i,1)) + det*0.5*temp/3.0 !here, we add the result of the convolution of the basis function with the right hand side of the poisson equation, which gives us the right hand side vector in the finite element equation.
    fu(t(i,2)) = fu(t(i,2)) + det*0.5*temp/3.0
    fu(t(i,3)) = fu(t(i,3)) + det*0.5*temp/3.0
  end do

  do i = 1,size(col2)
    col2(i)=0
    row2(i)=0
    val2(i)=0
  enddo
    
  j=1
  do i=1,size(col1) !reprocess information in row1 and col1, to compress them. redundancies are set to zero, to be ignored later.
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
  do i=1,nb !this loops through the b array and sets the corresponding row of l to all zeros except for the l(b(i),b(i)) spot. it
!  also sets the f(b(i) cell to whatever our boundary conditions say it has to be. this allows for correct evaluation of the edges.
    fu(b(i)) = u(i)
    do j = 1,nt*9
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
  allocate(ia(k),ja(k),arr(k)) !the points which were set to 0 in col2 are a waste of space and irrelevant to the calculation, so we need our final arrays not to include them. so they should be this size

  j=1
  do i = 1,size(col2) !makes arr and ja, the value and column arrays that will be interacting with pardiso. they mustn't have the redundant spaces that were set to zero in col2
    if(col2(i)/=0) then
      arr(j) = val2(i)
      ja(j) = col2(i)
      ia(j) = row2(i)
      j=j+1
    end if
  end do
  do i = 1,n !we want x to have an estimate of the solution to begin with... i start with 0, for lack of a better idea
    x(i) = 0
  end do
  temp = 1e-16
  i = nb/4
  q = i*4
  call mgmres_st(n,size(ia),ia,ja,arr,x,fu,i,q,temp,temp)
  endsubroutine weirdder

  subroutine gssolve(c,p,t,b,src_loc,x,src_val,areas,bound,order,guess) !Different G-S solver, which oversamples 3x the boundary
             !triangles, and spits out a srcloc array (centroids of triangles), and srcval (right hand side at those points)
  use functions
  implicit none
  integer::order
  integer *8::n,nt,nv,nb,i,j,k,q,nedge
  integer, dimension(3)::tri
  integer, dimension(:,:)::t
  integer, dimension(:)::b
  real *8, dimension(:,:)::p 
  integer, dimension(:),allocatable::row1,col1,row2,col2,ia,ja
  real *8, dimension(:,:),allocatable::src_loc
  real *8, dimension(:),allocatable::fu,val1,val2,arr,src_val,x,areas
  real *8, dimension(:)::bound,guess
  real *8, dimension(3,3)::a
  real *8::cgs,cx,cy,det,temp,c
  real *8, dimension(2)::centroid, pt1,pt2,pt3
  logical, dimension(:),allocatable::edgetri
  
  nt = size(t(:,1))
  nv = size(p(:,1))
  nb = size(b)
  n = nv !total number of elements will be the number of vertices
  allocate(fu(n))
  allocate(val1(nt*9),col1(nt*9),row1(nt*9),val2(nt*9),col2(nt*9),row2(nt*9),x(n),edgetri(nt))! allocate sizes of arrays
  q = 1
  do i = 1,nt !loop to build the row and col vectors
    do j = 0,8
      row1(q+j) = t(i,j/3+1)
      col1(q+j) = t(i,modulo(j,3)+1)
    end do
    q = q+9
  end do

  nedge = 0 !we want to count the number of edge triangles, to oversample them
  do i=1,nt
    tri = t(i,:)
    k=0
    do j=1,nb
      if(t(i,1)==b(j) .or. t(i,2)==b(j) .or. t(i,3)==b(j))then
        k = 1
      endif
    enddo
    if(k==1)then
      edgetri(i)=.true. !edgetri tells us which triangles are on the edge
      nedge = nedge+1 !nedge tells us how many triangles are on the edge
    else
      edgetri(i)=.false.
    endif
  enddo

  allocate(src_val(nt+(2)*nedge),src_loc(2,nt+(2)*nedge),areas(nt+2*nedge)) !we give our arrays space for oversampling
  do i = 1,n
    fu(i)=0
    src_val(i)=0 !initialise arrays to 0
  end do

  q = 1
  k = 1
  do i = 1,nt !loop through triangles to add entries into the vals vector, and build the f vector.
    do j =1,3 !use the same trick as in weirdder to compute the dot product of the basis functions
      a(j,1) = (p(t(i,j),1))
      a(j,2) = (p(t(i,j),2))
      a(j,3) = 1.0
    end do

    det = abs(threedet(a)) !this computes the determinant of that matrix, which is used to compute the area of the triangle in question
    call threeinv(a) !invert the array

    if(order==1)then !Order=1 is the initial solve, which is done by iteration, so the stiffness matrix needs not be modified
    do j=0,8
      val1(q+j) = (-1)*(a(1,j/3+1)*a(1,modulo(j,3)+1)+a(2,j/3+1)*a(2,modulo(j,3)+1))*det*0.5
    end do
    elseif(order==2)then !order=2 is the first derivative solve. This is done directly, so we must modify the stiffness matrix
    do j=0,8
      val1(q+j) = (-1)*(a(1,j/3+1)*a(1,modulo(j,3)+1)+a(2,j/3+1)*a(2,modulo(j,3)+1))*det*0.5 !initial matrix
      val1(q+j) = val1(q+j) + (-1)*1./9*stiff(p(t(i,1),1)+p(t(i,2),1)+&
        p(t(i,3),1),p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2),c)*det*0.5 !Modified piece
    end do
    endif
    q = q+9
    
    cx = (p(t(i,1),1)+p(t(i,2),1)+p(t(i,3),1))/(3.0) !centroid position
    cy = (p(t(i,1),2)+p(t(i,2),2)+p(t(i,3),2))/(3.0)
    cgs = (guess(t(i,1))+guess(t(i,2))+guess(t(i,3)))/3 !value of our guess, or the previous solution, @ the centroid
    if (order==1)then 
      temp = gsrhs(cgs,cx,cy,c) !Evaluate right hand side at the centroid
    elseif(order==2)then
      temp = gsrhsx(cgs,cx,cy,c)
    endif
    fu(t(i,1)) = fu(t(i,1)) + det*0.5*temp/3.0 !here, we add the result of the convolution of the basis function with the right hand side of the poisson equation, which gives us the right hand side vector in the finite element equation.
    fu(t(i,2)) = fu(t(i,2)) + det*0.5*temp/3.0
    fu(t(i,3)) = fu(t(i,3)) + det*0.5*temp/3.0
  
    if(edgetri(i))then !this is where we oversample. We generate three triangles per triangle, and midpoint rule each of them.
      centroid(1) = cx 
      centroid(2) = cy
      pt1(1) = (p(t(i,1),1)+p(t(i,2),1)+centroid(1))/3 !pt1 is the centroid of the first sub-triangle, pt2 etc.
      pt1(2) = (p(t(i,1),2)+p(t(i,2),2)+centroid(2))/3
      pt2(1) = (p(t(i,2),1)+p(t(i,3),1)+centroid(1))/3
      pt2(2) = (p(t(i,2),2)+p(t(i,3),2)+centroid(2))/3
      pt3(1) = (p(t(i,1),1)+p(t(i,3),1)+centroid(1))/3
      pt3(2) = (p(t(i,1),2)+p(t(i,3),2)+centroid(2))/3
      
      src_loc(1,k) = pt1(1) !location of the triangle is its centroid
      src_loc(2,k) = pt1(2)
      if(order==1)then 
          src_val(k) = gsrhs((cgs + guess(t(i,1)) + guess(t(i,2)))/3,pt1(1),pt1(2),c)*det*0.5/3 !midpoint rule the RHS at the centroid
      elseif(order==2)then
          src_val(k) = gsrhsx((cgs + guess(t(i,1)) + guess(t(i,2)))/3,pt1(1),pt1(2),c)*det*0.5/3
      endif
      areas(k) = det*0.5/3 !We also need the areas of the triangles, for error estimation. 
      src_loc(1,k+1) = pt2(1)
      src_loc(2,k+1) = pt2(2)
      if(order==1)then
          src_val(k+1) = gsrhs((cgs + guess(t(i,2)) + guess(t(i,3)))/3,pt2(1),pt2(2),c)*det*0.5/3
      elseif(order==2)then
          src_val(k+1) = gsrhsx((cgs + guess(t(i,2)) + guess(t(i,3)))/3,pt2(1),pt2(2),c)*det*0.5/3
      endif
      areas(k+1) = det*0.5/3
      src_loc(1,k+2) = pt3(1)
      src_loc(2,k+2) = pt3(2)
      if(order==1)then
          src_val(k+2) = gsrhs((cgs + guess(t(i,1)) + guess(t(i,3)))/3,pt3(1),pt3(2),c)*det*0.5/3
      elseif(order==2)then
          src_val(k+2) = gsrhsx((cgs + guess(t(i,1)) + guess(t(i,3)))/3,pt3(1),pt3(2),c)*det*0.5/3
      endif
      areas(k+2) = det*0.5/3
      k = k+3    
    else !If the triangle is not on the edge, we do the normal thing. No need to oversample.
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
  do i=1,size(col1) !reprocess information in row1 and col1, to compress them. redundancies are set to zero, to be ignored later.
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
  do i=1,nb !this loops through the b array and sets the corresponding row of l to all zeros except for the l(b(i),b(i)) spot. it also sets the f(b(i) cell to zero. this allows for correct evaluation of the edges.
    fu(b(i)) = bound(i)
    do j = 1,nt*9
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
  allocate(ia(k),ja(k),arr(k)) !the points which were set to 0 in col2 are a waste of space and irrelevant to the calculation, so we need our final arrays not to include them. so they should be this size

  j=1
  do i = 1,size(col2) !makes arr and ja, the value and column arrays that will be interacting with pardiso. they mustn't have the redundant spaces that were set to zero in col2
    if(col2(i)/=0) then
      arr(j) = val2(i)
      ja(j) = col2(i)
      ia(j) = row2(i)
      j=j+1
    end if
  end do
  do i = 1,n !we want x to have an estimate of the solution to begin with... i start with 0, for lack of a better idea
    x(i) = 0
  end do
  temp = 1e-16
  i = nb/4
  q = i*4
  call mgmres_st(n,size(ia),ia,ja,arr,x,fu,i,q,temp,temp)
  end subroutine gssolve
 
  function  linspace2(a,b,n) !equivalent of python linspace, except it doesn't include the endpoint
    use functions
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


  function threedet(m) !find the determinant of a 3x3 matrix
  use functions
  implicit none
  real *8, dimension(3,3):: m
  real *8:: threedet

  threedet =  m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + m(1,3)*m(2,1)*m(3,2)
  threedet = threedet - m(1,3)*m(2,2)*m(3,1) - m(1,1)*m(2,3)*m(3,2) - m(1,2)*m(2,1)*m(3,3)
  end function threedet  
  
  subroutine threeinv(m) !inverts 3x3 matrix. m is the matrix we want to invert
    real *8,dimension(3,3)::m
    real *8::det,a,b,c,d,e,f,g,h,i
    det = threedet(m)

    a = m(2,2)*m(3,3)-m(2,3)*m(3,2)
    b = (-1)*(m(2,1)*m(3,3)-m(2,3)*m(3,1))
    c = m(2,1)*m(3,2)-m(2,2)*m(3,1)
    d = (-1)*(m(1,2)*m(3,3)-m(1,3)*m(3,2))
    e = m(1,1)*m(3,3)-m(1,3)*m(3,1)
    f = (-1)*(m(1,1)*m(3,2)-m(1,2)*m(3,1))
    g = m(1,2)*m(2,3)-m(1,3)*m(2,2)
    h = (-1)*(m(1,1)*m(2,3)-m(1,3)*m(2,1))
    i = m(1,1)*m(2,2)-m(1,2)*m(2,1)

    m(1,1) = a
    m(1,2) = d
    m(1,3) = g
    m(2,1) = b
    m(2,2) = e
    m(2,3) = h
    m(3,1) = c
    m(3,2) = f
    m(3,3) = i
    m = m*(1/det)
  end subroutine threeinv
  
  subroutine distmesh(p,t,b,d1,d2,d3,d4,c,gam)
  use functions
  implicit none
  integer *8, parameter :: maxrecs = 1000000
  integer *8 :: j, nr, ios
  character(len=100) :: inputfile
  character(len=1) :: junk
  integer,dimension(:,:),allocatable::t
  real *8,dimension(:,:),allocatable::p
  real *8,dimension(3)::temp
  integer,dimension(:),allocatable::b
  real *8::d1,d2,d3,d4,c,gam
  nr = 0
  open(unit=1,file='infiles/params.txt')
  read(1,*) d1,d2,d3,d4,c,gam
  close(1)

    open(unit=1,file='infiles/p.txt')
        do j=1,maxrecs
            read(1,*,iostat=ios) junk
            if (ios /= 0) exit
            if (j == maxrecs) then
                stop
            endif
            nr = nr + 1
        enddo
        rewind(1)
        allocate(p(nr,2))
        do j=1,nr
            read(1,*) p(j,1),p(j,2)
        end do
        close(1)
        nr = 0
    open(unit=1,file='infiles/t.txt')
        do j=1,maxrecs
            read(1,*,iostat=ios) junk
            if (ios /= 0) exit
            if (j == maxrecs) then
                stop
            endif
            nr = nr + 1
        enddo
        rewind(1)
        allocate(t(nr,3))
        do j=1,nr
            read(1,*) temp(1),temp(2),temp(3)
    t(j,1) = int(temp(1))
            t(j,2) = int(temp(2))
            t(j,3) = int(temp(3))
        end do
        close(1)
        nr = 0
    open(unit=1,file='infiles/b.txt')
        do j=1,maxrecs
            read(1,*,iostat=ios) junk
            if (ios /= 0) exit
            if (j == maxrecs) then
                stop
            endif
            nr = nr + 1
        enddo
        rewind(1)
        allocate(b(nr))
        do j=1,nr
            temp(1)=0
            read(1,*) temp(1)
            b(j) = int(temp(1))
        end do
        close(1)

    end subroutine distmesh
  end module mesh
