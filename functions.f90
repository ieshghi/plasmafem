module functions
implicit none
contains

  function gsrhs(uat,x,y) !C=1
  implicit none
  real *8,intent(in)::x,y,uat
  real *8::gsrhs

  gsrhs = 3.0d0/4.0d0*uat/(x**2)+x**(3.0d0/2.0d0)
  endfunction gsrhs

  function gsrhsx(uat,x,y)
  implicit none
  real *8,intent(in)::x,y,uat
  real *8::gsrhsx

  gsrhsx = 3.0d0/2.0d0*x**(1.0d0/2.0d0) - 1.5d0*uat/(x**3)! + 3./4*uxat/(x**2)
  endfunction gsrhsx
  
  function stiff(x,y)
  implicit none
  real *8,intent(in)::x,y
  real *8::stiff

  stiff = 3./4/(x**2)

  endfunction stiff

  function gsrhsy(uat,x,y)
  implicit none
  real *8,intent(in)::x,y,uat
  real *8::gsrhsy

  gsrhsy = 0! 3./4*uyat/(x**2)
  endfunction gsrhsy
  
  function gsrhsxx(uat,uxat,x,y)
  implicit none
  real *8,intent(in)::x,y,uxat,uat
  real *8::gsrhsxx

  gsrhsxx =  0 - 3*uxat/(x**3) + 9./2*uat/(x**4) + 3./4*1./(x**0.5d0) !3./4*uxxat/(x**2) 
  endfunction gsrhsxx
  
  function gsrhsyy(uat,x,y)
  implicit none
  real *8,intent(in)::x,y,uat
  real *8::gsrhsyy

  gsrhsyy = 0 !3./4*uyyat/(x**2)
  endfunction gsrhsyy
  
  function gsrhsxy(uat,uyat,x,y)
  implicit none
  real *8,intent(in)::x,y,uat,uyat
  real *8::gsrhsxy

  gsrhsxy = 0 - 3./2*uyat/(x**3)!3./4*uxyat/(x**2) 
  endfunction gsrhsxy

  function exact(x,y,d1,d2,d3) !known solution to poisson equation (for debugging purposes)
  implicit none
  real *8::pi
  real *8,intent(in):: x,y,d1,d2,d3
  real *8:: exact
  pi = 4*atan(1.0)

  !exact = 1.0/8.0*x**4+d1+d2*x**2+d3*(x**4-4*(x**2)*(y**2))
  exact = (1.0/8.0*x**4+d1+d2*x**2+d3*(x**4-4*(x**2)*(y**2)))/sqrt(x)
  !exact = 1-(x*x+y*y)
  !exact = x*(1-x)*y*(1-y)
  !exact = (-1)*0.5*(1.0/pi)*(1.0/pi)*(sin(pi*x)*sin(pi*y))
  end function exact
  function foo(x,y,d1,d2,d3) !right side of the equation
  implicit none
  real,parameter::pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986
  real *8,intent(in):: x,y,d1,d2,d3
  real *8:: foo

  !foo = 3.0/2.0*x**2+2*d2+d3*(12*x**2-4*2*y**2-4*2*x**2)
  foo = (24*d1+(x**2)*(24*d2+(x**2)*(24*d3+35)-96*d3*y**2))/(32*x**(5.0d0/2.0d0))
  !foo = -4.0
  !foo = 2*(x*x-x+y*y-y)
  !foo = sin(pi*x)*sin(pi*y)
  end function foo

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

  !exactx = 1.0/2.0*x**3+d2*2*x+d3*(4*x**3-8*x*y**2)
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
  exactyy = (-8)*d3*x**(3.0d0/2.0d0)

  endfunction exactyy

endmodule functions
