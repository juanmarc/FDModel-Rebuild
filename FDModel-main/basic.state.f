C*****7***************************************************************72

	  subroutine basicState(n,m,nm1,mm1,radius,theta,zeta,psi,
     &                      a,b,c,wdim,dx,dy,idimy)

C*****7***************************************************************72
C  This routine will establish the initial vorticity profile for the
C  run.  Once the profile is computed, it will then call poissonSolver
C  to invert the vorticity profile and get the streamfunction profile
C*****7***************************************************************72
C  Variables used by basicState:
C
C*****7***************************************************************72
C** Variable declarations

      implicit none
      
      integer i,j,m,n,mm1,nm1,k,wdim,idimy
      
      double precision radius(m,n),theta(m,n),zeta(m,n),psi(m,n)
      double precision a(mm1),b(mm1),c(mm1),dx,dy,integral,correction
      
      double precision sumCos,r,fn1,fn2,fn3,s1,s2,s3,area
      
C*****7***************************************************************72
C** Constant declarations.  r1,r2,r3 and r4 are used for the ring
C** configuration.  zeta1,zeta2,zeta2 and zetaZero are also used
C** in the ring configuration

      double precision r1,r2,r3,r4
      parameter(r1=15.0D3,r2=22.5D3,r3=25.0D3,r4=32.5D3)
      
      double precision zero,one,two,three,four,eight,oneFourth
      double precision zeta1,zeta2,zeta3,zetaZero,pi
      parameter(zero=0.0D0,one=1.0D0,two=2.0D0,three=3.0D0)
      parameter(four=4.0D0,eight=8.0D0,onefourth=2.5D-1)
      parameter(zeta1=4.1825D-4,zeta2=7.000808D-3)
      parameter(zeta3=1.0D-5,zetaZero=zero)
      
C*****7***************************************************************72
C** Calculate the vorticity field at each gridpoint.
C** Use the ring vorticity profile found in Schubert et al.

      do 200, j=1,n
        do 210 i=1,m
          
          sumCos = zero
          
          do 220 k=1,8
            sumCos = sumCos + dcos(k*theta(i,j))
 220      continue
 
          r=radius(i,j)
          if (r .le. r1) then
            zeta(i,j) = zetaZero + zeta1 + zero
          elseif ((r .gt. r1) .and. (r .le. r2)) then
            fn1 = (r-r1)/(r2-r1)
            fn2 = (r2-r)/(r2-r1)
            fn3 = (r2-r+1)/(r2-r1)
            s1 = one - three*(fn1*fn1) + two*(fn1*fn1*fn1)
            s2 = one - three*(fn2*fn2) + two*(fn2*fn2*fn2)
            s3 = one - three*(fn3*fn3) + two*(fn3*fn3*fn3)
            zeta(i,j) = zetaZero+zeta1*s1+zeta2*s2+zeta3*sumCos*s3
          elseif ((r .gt. r2) .and. (r .le. r3)) then
            zeta(i,j) = zetaZero + zeta2 + zeta3*sumCos
          elseif ((r .gt. r3) .and. (r .le. r4)) then
            fn2 = (r-r3)/(r4-r3)
            fn3 = (r-r3)/(r4-r3)
            s2 = one - three*(fn2*fn2) + two*(fn2*fn2*fn2)
            s3 = one - three*(fn3*fn3) + two*(fn3*fn3*fn3)
            zeta(i,j) = zetaZero + zeta2*s2 + zeta3*sumCos*s3
          else
            zeta(i,j) = zetaZero + zero + zero
          endif
                  
 210    continue
 200  continue

C*****7***************************************************************72
C** Equate the boundaries so that the field is doubly periodic

      do 230 j=1,n
        zeta(m,j) = zeta(1,j)
 230  continue

      do 240 i=1,m
        zeta(i,n) = zeta(i,1)
 240  continue
 
C*****7***************************************************************72
C** Ensure the net circulation is zero.  First, compute the area of 
C** the domain

      area = (mm1*dx)*(nm1*dy)

C*****7***************************************************************72
C** Then compute the area integrated vorticity divided by the total
C** area to get the circulation correction factor.  The integral is
C** computed via trapezoidal quadrature

      integral = zero
      do 250 j=1,nm1
        do 260 i=1,mm1
          integral = integral + oneFourth*dx*dy*
     &               (zeta(i,j)+zeta(i,j+1)+zeta(i+1,j)+zeta(i+1,j+1))
 260    continue
 250  continue
      correction = integral/area
      
C*****7***************************************************************72
C** Then use the correction factor to adjust the vorticity to have
C** zero net circulation

      do 270 j=1,n
        do 280 i=1,m
          zeta(i,j) = zeta(i,j) - correction
 280    continue
 270  continue
 
C*****7***************************************************************72
C** Now invert the vorticity field to get the streamfunction field

      call poissonSolver(m,n,mm1,nm1,idimy,
     &                   zeta,psi,a,b,c,wdim,dx,dy)

C*****7***************************************************************72

      return
      end