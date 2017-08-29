module functions
implicit none
contains

  function gsrhs(uat,x,y,c) !C=1, right hand side of Grad Shafranov equation after switching to u=psi/sqrt(x)
  implicit none
  real *8,intent(in)::x,y,c,uat
  real *8::gsrhs

  gsrhs = 3.0d0/4.0d0*uat/(x**2)+c*x**(3.0d0/2.0d0)
  endfunction gsrhs

  function gsrhsx(uat,x,y,c) !Right hand side of x derivative
  implicit none
  real *8,intent(in)::x,y,c,uat
  real *8::gsrhsx

  gsrhsx = 3.0d0/2.0d0*c*x**(1.0d0/2.0d0) - 1.5d0*uat/(x**3)! + 3./4*uxat/(x**2)
  endfunction gsrhsx
  
  function stiff(x,y,c) !modification to the stiffness matrix (true for all derivatives)
  implicit none
  real *8,intent(in)::x,y,c
  real *8::stiff

  stiff = 3.0d0/4.0d0/(x**2)

  endfunction stiff

  function gsrhsy(uat,x,y,c) !y derivative of RHS
  implicit none
  real *8,intent(in)::x,y,c,uat
  real *8::gsrhsy

  gsrhsy = 0! 3./4*uyat/(x**2)
  endfunction gsrhsy
  
  function gsrhsxx(uat,uxat,x,y,c) !xx derivative of RHS
  implicit none
  real *8,intent(in)::x,y,c,uxat,uat
  real *8::gsrhsxx

  gsrhsxx =  0 - 3.0d0*uxat/(x**3) + 9.0d0/2.0d0*uat/(x**4) + 3.0d0/4.0d0*1.0d0*c/(x**0.5d0) !3./4*uxxat/(x**2) 
  endfunction gsrhsxx
  
  function gsrhsyy(uat,x,y,c) !yy derivative of RHS
  implicit none
  real *8,intent(in)::x,y,c,uat
  real *8::gsrhsyy

  gsrhsyy = 0 !3./4*uyyat/(x**2)
  endfunction gsrhsyy
  
  function gsrhsxy(uat,uyat,x,y,c) !xy derivative of RHS
  implicit none
  real *8,intent(in)::x,y,c,uat,uyat
  real *8::gsrhsxy

  gsrhsxy = 0 - 3.0d0/2.0d0*uyat/(x**3)!3./4*uxyat/(x**2) 
  endfunction gsrhsxy

  function exact(x,y,c,d1,d2,d3) !known solution to GS equation with 0 boundary conditions on a tokamak with parameters d1,d2,d3
  implicit none
  real *8::pi
  real *8,intent(in):: x,y,c,d1,d2,d3
  real *8:: exact
  pi = 4*atan(1.0)

  exact = (c/8.0*x**4+d1+d2*x**2+d3*(x**4-4*(x**2)*(y**2)))/sqrt(x)
  end function exact
  
  function exactx(x,y,c,d1,d2,d3) !x derivative of solution
  implicit none
  real *8::x,y,c,d1,d2,d3,exactx

  exactx = ((x**2)*(24*d2+7*(8*d3+c)*(x**2)-96*d3*y**2)-8*d1)/(16*x**(3.0d0/2.0d0))
  endfunction exactx

  function exacty(x,y,c,d1,d2,d3) !y derivative of solution
  implicit none
  real *8::x,y,c,d1,d2,d3,exacty

  exacty =  (-8)*d3*y*x**(3.0d0/2.0d0)
  endfunction exacty

  function exactxx(x,y,c,d1,d2,d3) !xx derivative of solution
  implicit none
  real *8::x,y,c,d1,d2,d3,exactxx

  exactxx = (24*d1+(x**2)*(24*d2+(x**2)*35*(8*d3+c)-96*d3*(y**2)))/(32*(x**(5.0d0/2.0d0)))
  endfunction exactxx

  function exactxy(x,y,c,d1,d2,d3) !xy derivative of solution
  implicit none
  real *8::x,y,c,d1,d2,d3,exactxy

  exactxy = (-12)*d3*sqrt(x)*y
  endfunction exactxy

  function exactyy(x,y,c,d1,d2,d3) !yy derivative of solution
  implicit none
  real *8::x,y,c,d1,d2,d3,exactyy
  exactyy = (-8)*d3*x**(3.0d0/2.0d0)

  endfunction exactyy

endmodule functions
