C*****7***************************************************************72

	  subroutine poissonSolver(m,n,mm1,nm1,idimy,
     &                         zeta,psi,a,b,c,wdim,dx,dy)

C*****7***************************************************************72
C Given the values of the vorticity at each gridpoint, this routine
C will solve the poisson equation:
C               delSquared(psi) = zeta
C to get the value of psi at each gridpoint.  It accomplished this by
C calling a routine pois which is the actual solver.  This particular
C routine simply sets up the appropriate variables to enable pois to 
C work correctly.  It then assigns the appropriate values to the 
C vorticity based on the output from pois.  The routines called from
C this routine and from pois were all obtained from the www site of
C the Scientific Computation Division of UCAR (the FISHPAK package).
C They have been modified somewhat, but the core of the routines remain
C*****7***************************************************************72
C Variables used by poissonSolver (and by pois):
C
C*****7***************************************************************72
C** Variable declarations

      implicit none
      integer i,j,n,m,wdim,iError
      integer nm1,mm1,idimy
      double precision a(mm1),b(mm1),c(mm1),zeta(m,n),psi(m,n)
      double precision ry(idimy,nm1),dx,dy
      double precision w(wdim)
      
C*****7***************************************************************72
C** Constant declarations

      integer iflg,nperod,mperod
      parameter (nperod=0,mperod=0)
      
      double precision dxsq,dysq
      
      iflg = 0
      dysq = dy*dy
      dxsq = dx*dx
      
C*****7***************************************************************72
C** Assign appropriate values to the ry array

      do 810 j=1,nm1
        do 820 i=1,idimy
          ry(i,j) = dysq*zeta(i,j)
 820    continue
 810  continue

C*****7***************************************************************72
C** With necessary values assigned, call the solver

      write(*,*) 'pois'
      call pois(iflg,nperod,nm1,mperod,mm1,a,b,c,idimy,ry,w)
      write(*,*) 'exit pois'

C*****7***************************************************************72
C** Based on values of the ry array on exit from pois, assign
C** the solution values to psi

      do 830 j=1,nm1
        do 840 i=1,idimy
          psi(i,j) = ry(i,j)
 840    continue
 830  continue

C*****7***************************************************************72
C** Now take care of the boundary conditions on PSI for a doubly
C** periodic domain

      do 850 j=1,n
        psi(m,j) = psi(1,j)
 850  continue

      do 860 i=1,m
        psi(i,n) = psi(i,1)
 860  continue

C*****7***************************************************************72

      return
      end