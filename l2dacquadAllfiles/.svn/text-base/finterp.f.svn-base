c
      subroutine finter(fcoarse,ffine,ncoarse,nfine)
c
C     Fourier interpolation subroutine for real-valued
C     periodic function.
C
C     INPUT:
C
C     fcoarse   = array of function values on coarse grid
C     ncoarse   = number of coarse grid points
C                IMPORTANT: ncoarse is assumed to be even
C     nfine    = number of fine grid points 
C                IMPORTANT: nfine/ncoarse must be an integer
C
C     OUTPUT:
C
C     ffine = array of function values on fine grid.
C             The first entry positions of the two arrays are 
C             chosen to coincide
C
      implicit real *8 (a-h,o-z)
      integer *4 ncoarse,nfine,lw,ipt,ip2
      real *8 fcoarse(ncoarse),ffine(nfine)
      complex *16, allocatable :: work(:)
c
c
      ip2 = 2*(ncoarse) + 10
      ipt = ip2 + 2*(nfine) + 10
      itot = ipt + nfine 
      allocate(work(itot))
c
      call zffti(ncoarse,work(1))
      call zffti(nfine,work(ip2))
c
c---- Write fcoarse data into workspace and compute forward transform
c
      do i = 1,ncoarse
	 work(ipt+i-1) = fcoarse(i)
      enddo
      call zfftf(ncoarse,work(ipt),work(1))
c
c     Embed in finer transform array in standard FFT
c     format, with zero padding as appropriate.
c
      do i = 0,ncoarse/2
	 work(ipt+i) = work(ipt+i)/ncoarse
      enddo
c
      do i = 1,ncoarse/2-1
	 work(ipt+nfine-i) = work(ipt+ncoarse-i)/ncoarse
      enddo
      do i = ncoarse/2 +1,nfine-ncoarse/2
	 work(ipt+i) = 0
      enddo
c
c     Compute inverse transform and write out fine data.
c
      call zfftb(nfine,work(ipt),work(ip2))
      do i = 1,nfine
	 ffine(i) = dreal(work(ipt+i-1))
      enddo
      return
      end
