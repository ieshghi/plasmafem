      subroutine l2dacquadwrapl(xs,ys,dsdt,rnx,rny,h,ns,
     1           ifcharge,sigma,ifdipole,mu,nover,nunder,
     2           ntj,iside,pot,potn,grad)
c
c     This subroutine computes the fields due to layer potentials
c     on a closed curve Gamma governed by the Laplace equation in 2D. 
c     Let D denote the interior of Gamma and E its exterior.
c   
c     The single layer potential S(sigma) is weakly singular and 
c     the limiting processes
c            x -> x_0   (x \in D, x_0 \in Gamma) 
c            x -> x_0   (x \in E, x_0 \in Gamma) 
c     are well-defined.

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
c     xs    x-component of equispaced boundary points with respect to 
c           some parametrization t,
c           assumed to be positively ordered (counterclockwise).
c     ys    y-component of equispaced boundary points with respect to 
c           some parametrization t,
c           assumed to be positively ordered (counterclockwise).
c     dsdt  derivative of arclength with respect to t 
c           at the corresponding points
c     rnx   x-component of outward normal
c           at the corresponding points
c     rny   y-component of outward normal
c           at the corresponding points
c     h     (equal) spacing in the parametrization t
c     ns    number of points
c     sigma single layer density at the corresponding points
c     mu    double layer density at the corresponding points
c     nover oversampling factor
c     nover undersampling factor
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
      real *8 xs(ns)
      real *8 ys(ns)
      real *8 dsdt(ns)
      real *8 rnx(ns),rny(ns)
      complex *16 sigma(ns)
      complex *16 mu(ns)
      complex *16 pot(ns), gg, zk, grad(2,ns)
      complex *16 potn(ns)
      complex *16 ztarg,hess(3)
      real *8, allocatable :: centers(:,:)
      real *8, allocatable :: xsover(:)
      real *8, allocatable :: ysover(:)
      real *8, allocatable :: source(:,:)
      real *8, allocatable :: dsdtover(:)
      real *8, allocatable :: rnxover(:)
      real *8, allocatable :: rnyover(:)
      real *8, allocatable :: scj(:)
      complex *16, allocatable :: sigmaover(:)
      complex *16, allocatable :: muover(:)
      complex *16, allocatable :: jexps(:,:)
      complex *16, allocatable :: dipvec(:)
c
      pi = 4.0d0*datan(1.0d0)
      nsover = ns*nover
      nt = ns/nunder
c
      allocate(centers(2,nt))
      allocate(scj(nt))
      allocate(jexps(0:ntj,nt))
c
      allocate(xsover(nsover))
      allocate(ysover(nsover))
      allocate(source(2,nsover))
      allocate(dsdtover(nsover))
      allocate(rnxover(nsover))
      allocate(rnyover(nsover))
      allocate(sigmaover(nsover))
      allocate(muover(nsover))
      allocate(dipvec(nsover))
c
c     create local expansion center points (the nt target points)
c     Assumes boundary is positively oriented.
c
ccc      write(6,*) ' nunder =',nunder
ccc      write(13,*) ' nunder =',nunder
ccc      write(6,*) ' nt =',nt
ccc      write(13,*) ' nt =',nt
      do jj = 1,nt
         j = nunder*jj-(nunder-1)
ccc         j = nunder*jj-(nunder/2)
         if (j.eq.1) then 
            h2 = (xs(2)-xs(ns))**2 + (ys(2)-ys(ns))**2 
         else if (j.eq.ns) then 
            h2 = (xs(1)-xs(ns-1))**2 + (ys(1)-ys(ns-1))**2 
         else 
            h2 = (xs(j+1)-xs(j-1))**2 + (ys(j+1)-ys(j-1))**2 
         endif 
         h2 = sqrt(h2)
         write(33,*) h2
c
c        h2 is approx 2*local boundary spacing.
c
         centers(1,jj) = xs(j) - iside*3*rnx(j)*h2
         centers(2,jj) = ys(j) - iside*3*rny(j)*h2
      enddo
c
      ifprint=1
      if (ifprint.eq.1) then
         do i = 1,ns
            write(11,*) xs(i),ys(i)
         enddo
         do i = 1,nt
            write(12,*) centers(1,i),centers(2,i)
         enddo
      endif
c
c     oversample curve and densities
c
      call finter(xs,xsover,ns,nsover)
      call finter(ys,ysover,ns,nsover)
      call finter(dsdt,dsdtover,ns,nsover)
      call finter(rnx,rnxover,ns,nsover)
      call finter(rny,rnyover,ns,nsover)
c
      if (ifcharge.eq.1) call finterc(sigma,sigmaover,ns,nsover)
      if (ifdipole.eq.1) call finterc(mu,muover,ns,nsover)
c
      ifprint=1
      if (ifprint.eq.1) then
         do i = 1,nsover
            write(21,*) xsover(i),ysover(i)
         enddo
      endif
c
c     scale densities to include (smooth) quadrature weights.
c
      do i = 1,nsover
         if (ifcharge.eq.1) 
     1      sigmaover(i) = sigmaover(i)*dsdtover(i)*h/nover
         if (ifdipole.eq.1) muover(i) = muover(i)*dsdtover(i)*h/nover
      enddo
c
      rscale = 1.0d0
      do i = 1,nsover
         source(1,i) = xsover(i)
         source(2,i) = ysover(i)
         if (ifdipole.eq.1) 
     1      dipvec(i) = -muover(i)*(dcmplx(rnxover(i),rnyover(i)))
      enddo
      if (ifprint.eq.1) then
         do i = 1,nsover
            write(24,*) dipvec(i)
         enddo
      endif
      do i = 1,nt
         do j = 0,ntj
            jexps(j,i) = 0.0d0
         enddo
      enddo
c
      islow = 0
      if (islow.eq.1) then
      do i = 1,nt
         scj(i) = rscale
         if (ifcharge.eq.1) then
            call l2dformta_add(ier,rscale,source,sigmaover,nsover,
     1        centers(1,i),ntj,jexps(0,i))
         endif
         if (ifdipole.eq.1) then
         call l2dformta_dp_add(ier,rscale,source,dipvec,
     1        nsover,centers(1,i),ntj,jexps(0,i))
         endif
      enddo
      endif
c
      iprec = 4
      if (islow .eq. 0) then
ccc      call cfmm2dacquad(ier,iprec,nsover,source,
ccc     1        ifcharge,sigmaover,ifdipole,dipvec,
ccc     1        nt,centers,jexps,ntj)
      call cfmm2dacquads(ier,iprec,nsover,source,
     1        ifcharge,sigmaover,ifdipole,dipvec,
     1        nt,centers,jexps,ntj,scj)
c
      endif
c
c     now evaluate local expansions at desired boundary points
c     using nearest (undersampled) expansion center.
c
      ifgrad = 1
      ifhess = 0
c
ccc      write(6,*) ' nunder =',nunder
ccc      write(13,*) ' nunder =',nunder
ccc      write(6,*) ' nt =',nt
ccc      write(13,*) ' nt =',nt
ccc      t0 = second()
      do j = 1,ns
         jj = (j+nunder-1)/nunder
ccc         jj = mod(j+nunder/2,ns)/nunder
         if (jj.eq.0) jj = nt
ccc         write(6,*)' j,jj are ',j,jj
         ztarg = dcmplx(xs(j),ys(j))
         rscale = scj(jj)
ccc         write(6,*)' jj rscale are ',jj,rscale
ccc         write(13,*)' jj rscale are ',jj,rscale
         call l2dtaeval(rscale,centers(1,jj),jexps(0,jj),
     1        ntj,ztarg,pot(j),ifgrad,grad(1,j),ifhess,hess)
c
         pot(j) = -pot(j)/(2*pi)
c
c     subtract diagonal limiting term
c
         if (ifdipole.eq.1) pot(j) = pot(j) + iside*mu(j)/2.0d0
c
         potn(j) = grad(1,j)*rnx(j) + grad(2,j)*rny(j)
         potn(j) = -potn(j)/(2*pi)
c
c     subtract diagonal limiting term
c
         if (ifcharge.eq.1) potn(j) = potn(j) - iside*sigma(j)/2.0d0
c
         grad(1,j) = -grad(1,j)/(2*pi)
         grad(2,j) = -grad(2,j)/(2*pi)
      enddo
ccc      t1 = second()
ccc      call prin2(' time for evals is *',t1-t0,1)
      return
      end
