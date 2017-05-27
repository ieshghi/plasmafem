c
c     Test acquad for SLP and DLP 
c     in two dimensions for the Laplace equation 
c     
c     We define the field in a domain R as that induced by an 
c     exterior point source and the field in the exterior domain
c     E as that induced by an interior point source. 
c     The boundary of R, denoted Gamma, is assumed to be positively
c     oriented and the normal is assumed to be outward from R.
c     Green's theorem states that
c
c            u = S[du/dn] - D[u]   inside
c            u = -S[du/dn] + D[u]   outside
c
c     assuming G = -log(r)/(2\pi).
c     
c     The sign change in Green's identity follows because
c     of our convention that the normal vector points into 
c     the exterior.
c
c     We check the formula above at off surface points. We then
c     take the limit from the inside/outside
c     and, denoting the PV of D by D^*, we have
c
c     (1/2) u = S[du/dn] - D^*[u]     inside
c     (1/2) u = -S[du/dn] + D^*[u]  outside.
c
c     We check the error in these formula.
c     We then check the limit of the normal derivative of 
c     Green's identities:
c
c     (1/2) du/dn = S'^*[du/dn] - D'[u]   inside
c     (1/2) du/dn = -S'^*[du/dn] + D'[u]   outside
c
      program testgreen
      implicit real *8 (a-h,o-z)
      parameter (nmax = 1000000)
      complex *16 uin(nmax)
      complex *16 dudnin(nmax)
      complex *16 uout(nmax)
      complex *16 dudnout(nmax)
      complex *16 sigma(nmax)
      complex *16 rmu(nmax)
      dimension x(nmax)
      dimension y(nmax)
      dimension rnx(nmax)
      dimension rny(nmax)
      dimension dsdt(nmax)
      complex *16 dlp(nmax)
      complex *16 dlph(nmax)
      complex *16 slp(nmax)
      complex *16 slph(nmax)
      complex *16 dlpout(nmax)
      complex *16 dlphout(nmax)
      complex *16 slpout(nmax)
      complex *16 slphout(nmax)
      complex *16 potnsin(nmax)
      complex *16 potndin(nmax)
      complex *16 potnsout(nmax)
      complex *16 potndout(nmax)
      complex *16 gradin(2,nmax)
      complex *16 gradout(2,nmax)
      complex *16 gradsin(2,nmax)
      complex *16 gradsout(2,nmax)
      complex *16 graddin(2,nmax)
      complex *16 graddout(2,nmax)
      complex *16 dummy(nmax)
      complex *16 eye,z,h1,h0,zk,u,ux,uy
      complex *16 ssum,dsum
      external glapgrad
      external glap
c
      call prini(6,13)
      pi = 4.0d0*datan(1.0d0)
      eye = dcmplx(0.0d0,1.0d0)
c
c     create smooth geometry, define number of discretization points
c     and step size h in parametrization.
c
      a = 1.2d0
      xa = 1.2d0
      xb = 0.5d0
      nosc = 5
      b = 1.3d0
      n = 500
c
      h = 2.0d0*pi/n
c
c
c     define points in exterior and interior of R.
c
      xout = 5.75d0
      yout = 0.2d0
      xin = 0.75d0
      yin = 0.2d0
c
      write(15,*) xout, yout
      write(15,*) xin, yin
c
      do i = 1,n
         theta = 2.0d0*pi*(i-1)/n
c
c        alex's domain R
c
         wave = 1.0d0 + 0.3d0*dsin(nosc*theta)
         dwavedt = nosc*0.3d0*dcos(nosc*theta)
         wave = 1.0d0 
         dwavedt = 0.0d0
c
         x(i) = dcos(theta)*wave
         y(i) = dsin(theta)*wave
         dx = -sin(theta)*wave + cos(theta)*dwavedt
         dy = cos(theta)*wave + sin(theta)*dwavedt
         rnorm = dsqrt(dx*dx + dy*dy)
         rnx(i) = dy/rnorm
         rny(i) = -dx/rnorm
         dsdt(i) = rnorm
c
c        create dudn and u for interior test (dudnin, uin)
c
         rr2 = (x(i)- xout)**2 + (y(i) - yout)**2
	       rr = dsqrt(rr2)
         u = dlog(rr)/(2*pi)
         h1 = 1.0d0/rr/(2*pi)
         ux = h1*(x(i)-xout)/rr
         uy = h1*(y(i)-yout)/rr
         dudnin(i) = ux*rnx(i) + uy*rny(i)
         uin(i) = u
         gradin(1,i) = ux
         gradin(2,i) = uy
c
c        create dudn and u for exterior test (dudnout, uout)
c
         rr2 = (x(i)- xin)**2 + (y(i) - yin)**2
	 rr = dsqrt(rr2)
         u = dlog(rr)/(2*pi)
         h1 = 1.0d0/rr/(2*pi)
         ux = h1*(x(i)-xin)/rr
         uy = h1*(y(i)-yin)/rr
         gradout(1,i) = ux
         gradout(2,i) = uy
         dudnout(i) = ux*rnx(i) + uy*rny(i)
         uout(i) = u
      enddo
c
c     Check Green's formula at interior and exterior points.
c     Note: when evaluating the doulbe layer kernel, we want
c     ux, uy below to be the derivatives of u
c     with respect to the source (x(i),y(i)) not the target,
c     hence the extra minus sign.
c
c
      ssumin = 0.0d0
      dsumin = 0.0d0
      ssumout = 0.0d0
      dsumout = 0.0d0
      do i = 1,n
         w = dsdt(i)
         rr2 = (xin -x(i))**2 + (yin - y(i))**2
	 rr = dsqrt(rr2)
         u = -dlog(rr)/(2*pi)
         h1 = -1.0d0/rr/(2*pi)
         ux = -h1*(xin-x(i))/rr
         uy = -h1*(yin-y(i))/rr
         ssumin = ssumin + w*dudnin(i)*u
         dsumin = dsumin + w*uin(i)*(ux*rnx(i) + uy*rny(i))
c
         rr2 = (xout -x(i))**2 + (yout - y(i))**2
	 rr = dsqrt(rr2)
         u = -dlog(rr)/(2*pi)
         h1 = -1.0d0/rr/(2*pi)
         ux = -h1*(xout-x(i))/rr
         uy = -h1*(yout-y(i))/rr
         ssumout = ssumout + w*dudnout(i)*u
         dsumout = dsumout + w*uout(i)*(ux*rnx(i) + uy*rny(i))
      enddo
      ssumin = ssumin*h
      dsumin = dsumin*h
      ssumout = ssumout*h
      dsumout = dsumout*h
      rr2 = (xin -xout)**2 + (yin - yout)**2
      rr = dsqrt(rr2)
      h0 = dlog(rr)/(2*pi)
      call prin2( ' uexact is *',h0,2)
      call prin2( ' ssumin is *',ssumin,1)
      call prin2( ' dsumin is *',dsumin,1)
      call prin2( ' ssumout is *',ssumout,1)
      call prin2( ' dsumout is *',dsumout,1)
      call prin2( ' interior S-D is *',ssumin-dsumin,1)
      call prin2( ' exterior -(S-D) is *',-(ssumout-dsumout),1)
c
c     compute single layer potential at discretization pts.
c
      nover = 1
      nunder = 1
      ntj = 6
      norder = 8
      ifcharge = 1
      ifdipole = 0
      iin = 1
      iout = -1
c
c     interior SLP calls
c
      call l2dacquadwrap(x,y,dsdt,rnx,rny,h,n,ifcharge,dudnin,
     1     ifdipole,uin,nover,nunder,ntj,iin,slp,potnsin,gradsin)
      call prin2(' x is *',x,n) 
      call prin2(' y is *',y,n) 
      call prin2(' dsdt is *',dsdt,n) 
      call prin2(' rnx is *',rnx,n) 
      call prin2(' rny is *',rny,n) 
      call prin2(' dudnin is *',dudnin,2*n) 
      call prin2(' uin is *',uin,2*n) 
      call prin2(' h is *',h,1) 
      call prinf(' ifcharge is *',ifcharge,1) 
      call prinf(' ifdipole is *',ifdipole,1) 
      call prinf(' nover is *',nover,1) 
      call prinf(' nunder is *',nunder,1) 
      call prinf(' ntj is *',ntj,1) 
      call prinf(' iin is *',iin,1) 
cccc      call l2dqbxwrap(x,y,dsdt,rnx,rny,h,n,ifcharge,dudnin,
cccc     1     ifdipole,uin,nover,nunder,ntj,iin,slp,potnsin,gradsin)
ccc      call slpquaddirectc(ier,norder,slph,x,y,dsdt,dudnin,h,n,
ccc     1     glap,p1,zk,i1,i2)
      call prin2(' slp is *',slp,20) 
      call prin2(' gradsin is *',gradsin,20) 
c
c     exterior SLP calls
c
      call prin2(' x = *',x,n)
      call prin2(' y = *',y,n)
      call prin2(' rnx = *',rnx,n)
      call prin2(' rny = *',rny,n)
      call prin2(' h = *',h,1)
      call prin2(' dudnin = *',dudnout,1)
      call prinf(' nover = *',nover,1)
      call prinf(' nunder = *',nunder,1)
ccc      call l2dacquadwrap(x,y,dsdt,rnx,rny,h,n,ifcharge,dudnout,
ccc     1    ifdipole,uout,nover,nunder,ntj,iout,slpout,potnsout,gradsout)
ccc      call l2dqbxwrap(x,y,dsdt,rnx,rny,h,n,ifcharge,dudnout,
ccc     1     ifdipole,uout,nover,nunder,ntj,iout,slpout,potnsout,gradsout)
ccc      call slpquaddirectc(ier,norder,slphout,x,y,dsdt,dudnout,h,n,
ccc     1     glap,p1,zk,i1,i2)
c
cc      call prin2(' gradsout is *',gradsout,2*2*n) 
ccc      call prin2(' slp is *',slp,2*n) 
ccc      call prin2(' slph is *',slph,2*n) 
c
      ifcharge = 0
      ifdipole = 1
c
c     interior DLP calls
c
      call l2dacquadwrap(x,y,dsdt,rnx,rny,h,n,ifcharge,dudnin,
     1     ifdipole,uin,nover,nunder,ntj,iin,dlp,potndin,graddin)
cccc      call l2dqbxwrap(x,y,dsdt,rnx,rny,h,n,ifcharge,dudnin,
cccc     1     ifdipole,uin,nover,nunder,ntj,iin,dlp,potndin,graddin)
ccc      call dlpquaddirectc(ier,norder,dlph,x,y,rnx,rny,dsdt,uin,h,n,
ccc     1     glapgrad,p1,zk,i1,i2)
c
c     exterior DLP calls
c
ccc      call l2dacquadwrap(x,y,dsdt,rnx,rny,h,n,ifcharge,dudnout,
ccc     1     ifdipole,uout,nover,nunder,ntj,iout,dlpout,potndout,graddout)
ccc      call l2dqbxwrap(x,y,dsdt,rnx,rny,h,n,ifcharge,dudnout,
ccc     1    ifdipole,uout,nover,nunder,ntj,iout,dlpout,potndout,graddout)
ccc      call dlpquaddirectc(ier,norder,dlphout,x,y,rnx,rny,dsdt,uout,h,n,
ccc     1     glapgrad,p1,zk,i1,i2)
c
cc      call prin2(' graddout is *',graddout,2*2*n) 
cc      call prin2(' gradout is *',gradout,2*2*n) 
ccc      call prin2(' dlp is *',dlp,2*n) 
ccc      call prin2(' dlph is *',dlph,2*n) 
c
c
      err = 0.0d0
      err2 = 0.0d0
      err3 = 0.0d0
      err4 = 0.0d0
      rnorm = 0.0d0
      errah = 0.0d0
      do i = 1,n
         rnorm = rnorm + abs(uin(i))**2
         err =  err+abs(real(slp(i)-dlp(i)-uin(i)/2.0d0))**2
         err2 = err2+abs(real(slph(i)-dlph(i)-uin(i)/2.0d0))**2
         err3 = err3+abs(real(-slpout(i)+dlpout(i)-uout(i)/2.0d0))**2
         err4 = err4+abs(real(-slphout(i)+dlphout(i)-uout(i)/2.0d0))**2
      enddo
      call prin2(' green id error acquad int= *',dsqrt(err/rnorm),1)
      call prin2(' green id error alpert int = *',dsqrt(err2/rnorm),1)
      call prin2(' green id error acquad ext = *',dsqrt(err3/rnorm),1)
      call prin2(' green id error alpert ext = *',dsqrt(err4/rnorm),1)
      rnormg = 0.0d0
      errg = 0.0d0
      do i = 1,n
         rnormg = rnormg + abs(gradin(1,i))**2 + abs(gradin(2,i))**2
         errg = errg+abs(real(gradsin(1,i)-graddin(1,i)-gradin(1,i)))**2
         errg = errg+abs(real(gradsin(2,i)-graddin(2,i)-gradin(2,i)))**2
      enddo
      call prin2(' grad error acquad int= *',dsqrt(errg/rnormg),1)
      rnormg = 0.0d0
      errg = 0.0d0
      do i = 1,n
         rnormg = rnormg + abs(gradout(1,i))**2 + abs(gradout(2,i))**2
         errg = errg+abs(real(-gradsout(1,i)+
     1                   graddout(1,i)-gradout(1,i)))**2
         errg = errg+abs(real(-gradsout(2,i)+
     1                   graddout(2,i)-gradout(2,i)))**2
      enddo
      call prin2(' grad error acquad ext= *',dsqrt(errg/rnormg),1)
c
c     check dudn formula from interior
c
c
      err = 0.0d0
      rnorm = 0.0d0
      do i = 1,n
         rnorm = rnorm + abs(dudnin(i))**2
         err =  err+abs(real(potnsin(i)-potndin(i)-dudnin(i)/2.0d0))**2
      enddo
      call prin2(' dudn error int = *',dsqrt(err/rnorm),1)
c
c      check dudn formula from exterior
c
      err = 0.0d0
      rnorm = 0.0d0
      do i = 1,n
         rnorm = rnorm + abs(dudnout(i))**2
         err =  err+abs(real(-potnsout(i)+
     1               potndout(i)-dudnout(i)/2.0d0))**2
      enddo
      call prin2(' dudn error ext = *',dsqrt(err/rnorm),1)

      stop
      end
c
c
      subroutine glap(xt,yt,x,y,p1,zk,i1,i2,gg)
      implicit real *8 (a-h,o-z)
      complex *16 gg,h0,h1,zk,z,eye,zs
c
c     compute Green's function -log(r)/(2pi) for Laplace equation.
c
      pi = 4.0d0*datan(1.0d0)
c
      rr2 = (x-xt)**2 + (y-yt)**2
      rr = dsqrt(rr2)
      gg = -dlog(rr)/(2*pi)
      return
      end
c
      subroutine glapgrad(xt,yt,x,y,p1,zk,i1,i2,gx,gy)
      implicit real *8 (a-h,o-z)
      complex *16 eye,gx,gy,h0,h1,zk,z,zs
c
c     compute gradient of Green's function for Laplace equation.
c
      pi = 4.0d0*datan(1.0d0)
c
      rr2 = (x-xt)**2 + (y-yt)**2
      rr = dsqrt(rr2)
      h1 = -1.0d0/rr/(2*pi)
      gx = -h1*(xt-x)/rr
      gy = -h1*(yt-y)/rr
      return
      end

