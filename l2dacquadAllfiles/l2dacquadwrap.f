      subroutine l2dacquadwrap(src_loc,src_val,targ_loc,targ_norm,
     1           n,m,ntj,iside,pot)
c
c     This subroutine computes the fields due to layer potentials
c     on a closed curve Gamma governed by the Laplace equation in 2D. 
c     Let D denote the interior of Gamma and E its exterior.
c   
c     The double layer potential D(mu) is singular and 
c     the limiting processes
c            x -> x_0   (x \in D, x_0 \in Gamma) 
c            x -> x_0   (x \in E, x_0 \in Gamma) 
c     yield
c     -1/2 mu + D^*(mu)  from the interior and 
c     +1/2 mu + D^*(mu)  from the exterior, where D^*
c     denotes the principal value integral.
c
c     The normal derivative of S(sigma) denoted S'(sigma)
c     takes the value
c     +1/2 sigma + S'^*(sigma)  from the interior and 
c     -1/2 sigma + S'^*(sigma)  from the exterior, where S'^*
c     denotes the principal value integral.
c
c     The normal derivative of D(mu) denoted D'(mu)
c     is hypersingular and has no jump (the limiting finite part
c     integral is the same from either side).
c
c     This subroutine computes the principal value integrals: 
c
c     pot = S(sigma)(x) + D^*(mu)(x) for x \in Gamma
c     potn = S'^*(sigma)(x) + D'(mu)(x) for x \in Gamma
c
c     It computes GRAD as the one-sided limit from either
c     the interior or exterior. Thus, 
c
c     if iside = 1
c
c     grad = lim      \nabla (S(sigma) + D(mu))   .
c            x -> x_0   
c            x \in D, x_0 \in Gamma) 
c
c     If iside = -1
c
c     grad = lim      \nabla (S(sigma) + D(mu)) 
c            x -> x_0   
c            x \in E, x_0 \in Gamma) 
c
c     Method:
c     A variant of the fast multipole algorithm is used to compute a
c     local expansion of the solution at a set of 
c     targets near the boundary, from which the boundary values
c     are computed. 
c     If iside = 1, the targets are located in the interior of the 
c     closed boundary curve defined by (xs,ys,) which is assumed to be 
c     positively oriented with (rnx,rny) the outward normal to the 
c     closed curve.
c     If iside = -1, the targets are located in the exterior of the 
c     closed boundary curve defined by (xs,ys,).
c
c     INPUT:
c 
c     src_loc    location of the points where the source density is measured
c                2 by n matrix, the first row corresponding to x values, the
c                second row to y values
c     src_val    value of the source density at the src_loc locations
c                vector with length n
c                note that the density includes the Jacobian of the transformation
c                from the computational domain to the square
c     targ_loc   location of the points where the potential is evaluated
c                2 by m matrix, the first row corresponding to x values, the
c                second row to y values
c     targ_norm  components of the outward  normal vector at the m target locations
c                2 by m matrix, the first row corresponding to x components,
c                the second row to y components
c     n          number of source points
c     m          number of target points
c     ntj   order of local expansion
c     iside 1 means interior, -1 means exterior
c
c     OUTPUT:
c
c     pot   is the principal value of the potential on the boundary.
c     potn  is the principal value/finite part of the normal
c            derivative of the potential on the boundary.
c     grad  is the limiting value of the gradient of the potential 
c            at the corresponding points 
c            (from the side determined by iside).
c
      implicit real *8 (a-h,o-z)
      integer iside
      integer n
      integer m
      integer ntj
      integer ier
      integer iprec
      real *8 src_loc(2,n)
      real *8 src_val(n)
      real *8 targ_loc(2,m)
      real *8 targ_norm(2,m)
      real *8 centers(2,m)
      real *8 rnx(m),rny(m)
      real *8 source(2,n)
      real *8 sigma(n)
      complex *16 dipvec(n)
      complex *16 pot(m)
      complex *16 grad(2,m)
      complex *16 ztarg,hess(3)
      real *8 scj(m)
      complex *16 jexps(ntj+1,m)
c
c
      do j=1,n
         sigma(j)=0.0d0
      enddo   
c
      pi = 4.0d0*datan(1.0d0)
      rnx=targ_norm(1,:)
      rny=targ_norm(2,:)
c      
c
      do j = 1,m
         if (j.eq.1) then 
            h2 = (targ_loc(1,2)-targ_loc(1,m))**2
     1          +(targ_loc(2,2)-targ_loc(2,m))**2 
         else if (j.eq.m) then 
            h2 = (targ_loc(1,1)-targ_loc(1,m-1))**2 
     1          +(targ_loc(2,1)-targ_loc(2,m-1))**2 
         else 
            h2 = (targ_loc(1,j+1)-targ_loc(1,j-1))**2 
     1          +(targ_loc(2,j+1)-targ_loc(2,j-1))**2 
         endif 
         h2 = sqrt(h2)
         write(33,*) h2
c
c        h2 is approx 2*local boundary spacing.
c
c        Now locate QBX expansion centers
         centers(1,j) = targ_loc(1,j) - iside*3*rnx(j)*h2
         centers(2,j) = targ_loc(2,j) - iside*3*rny(j)*h2
      enddo
c
c
c     x and y components of the gradient of the Green's function
c     that is part of the source density
c
c
      do i=1,n
         dipvec(i)=src_val(i)*(dcmplx(1.0d0,0.0d0))
      enddo
c
c
c     Initialize jexps, the Taylor expansions calculated 
c     by the FMM routine
c
c
      do i = 1,m
         do j = 1,ntj+1
            jexps(j,i) = 0.0d0
         enddo
      enddo
c
c
      iprec = 4
c
      call cfmm2dacquads(ier,iprec,n,src_loc,
     1        0,sigma,1,dipvec,
     1        m,centers,jexps,ntj,scj)
c
c     now evaluate local expansions at desired boundary points
c     using nearest (undersampled) expansion center.
c
c      WRITE(*,*) jexps
      ifgrad = 0
      ifhess = 0
c
      do j = 1,m
         ztarg = dcmplx(targ_loc(1,j),targ_loc(2,j))
         rscale = scj(j)
         call l2dtaeval(rscale,centers(1,j),jexps(1,j),
     1        ntj,ztarg,pot(j),ifgrad,grad(1,j),ifhess,hess)
c
         pot(j) = -pot(j)/(2*pi)
      enddo
c
c      WRITE(*,*) pot
c
      return
      end
