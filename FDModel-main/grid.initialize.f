C*****7***************************************************************72

	  subroutine gridInitialize(n,m,dx,dy,radius,theta)

C*****7***************************************************************72
C This routine will set up the grid to be used in the run.  It will
C compute the radius and angle at each grid point.  The angle is in
C reference to the positive x axis.  That is, due east is 0, due north
C is pi/2, due west is pi, and due south is 3pi/2.  Also note that
C the center of the grid is given a radius of 0, and all other radii
C are simply the radial distance from the center.
C
C !!!!NOTE:
C     For this routine to work correctly, n and m must be ODD.
C*****7***************************************************************72
C Variables used by gridInitialize
C
C*****7***************************************************************72
C** Variable declarations

      implicit none
      include 'fd.common'
      integer xmdpt,ymdpt,i,j
      double precision x,y
      
C*****7***************************************************************72
C** Constant declarations

      double precision pi,one,two,three,four,zero
      parameter (zero=0.0D0,one=1.0D0,two=2.0D0,three=3.0D0,four=4.0D0)
      
      pi = four*datan(one)
      
C*****7***************************************************************72
C** Determine the midpoints in x and y

      xmdpt = int((m+1)/2)
      ymdpt = int((n+1)/2)
      
C*****7***************************************************************72
C** Calculate the radius and angle at each point.  First determine
C** the y distance in relation to the center and then determine the
C** x distance in relation to the center.  From these values, get the 
C** radius and angle.

      do 100 j=1,n
        y = (one*(j-ymdpt))*dy
        do 110 i=1,m
          x = (one*(i-xmdpt))*dx
          radius(i,j) = dsqrt(one*(x*x + y*y))
          
C** since atan will return a value in the interval (-pi/2,pi/2), we 
C** must adjust the result so that if falls in the interval [0,2pi]          
        
          if ((x .gt. zero) .and. (y .ge. zero)) then
            theta(i,j) = datan(one*y/x)
          elseif ((x .gt. zero) .and. (y .le. zero)) then
            theta(i,j) = two*pi + datan(one*y/x)
          elseif (x .lt. zero) then
            theta(i,j) = pi + datan(one*y/x)
          elseif ((x .eq. zero) .and. (y .gt. zero)) then
            theta(i,j) = pi/two
          elseif ((x .eq. zero) .and. (y .lt. zero)) then
            theta(i,j) = three*pi/two
          else
            theta(i,j) = zero
          endif
          
 110    continue 
 100  continue      

C*****7***************************************************************72
      return
      end