      subroutine pois (iflg,nperod,n,mperod,m,a,b,c,idimy,y,ierror,w)
c
c
c***********************************************************************
c
c          version  2  october 1976  including errata  october 1976
c
c         documentation for this program is given in
c
c        efficient fortran subprograms for the solution of
c            elliptic partial differential equations
c
c                              by
c
c          paul swarztrauber   and  roland sweet
c
c             technical note tn/ia-109   july 1975
c
c       national center for atmospheric research  boulder,colorado 80307
c
c        which is sponsored by the national science foundation
c
c***********************************************************************
c
c
c
c     subroutine pois solves the linear system of equations
c
c          a(i)*x(i-1,j) + b(i)*x(i,j) + c(i)*x(i+1,j)
c
c          + x(i,j-1) - 2.*x(i,j) + x(i,j+1) = y(i,j)
c
c               for i = 1,2,...,m  and  j = 1,2,...,n.
c
c     the indices i+1 and i-1 are evaluated modulo m, i.e.,
c     x(0,j) = x(m,j) and x(m+1,j) = x(1,j), and x(i,0) may be equal to
c     0, x(i,2), or x(i,n) and x(i,n+1) may be equal to 0, x(i,n-1), or
c     x(i,1) depending on an input parameter.
c
c
c     * * * * * * * * * *     on input     * * * * * * * * * *
c
c     iflg
c       = 0  on initial entry to pois or if n and nperod are changed
c            from previous call.
c       = 1  if n and nperod are unchanged from previous call to pois.
c
c       note0  a call with iflg = 1 is about 1 percent faster than a
c              call with iflg = 0  .
c
c     nperod
c       indicates the values which x(i,0) and x(i,n+1) are assumed to
c       have.
c
c       = 0  if x(i,0) = x(i,n) and x(i,n+1) = x(i,1).
c       = 1  if x(i,0) = x(i,n+1) = 0  .
c       = 2  if x(i,0) = 0 and x(i,n+1) = x(i,n-1).
c       = 3  if x(i,0) = x(i,2) and x(i,n+1) = x(i,n-1).
c
c     n
c       the number of unknowns in the j-direction.  n is dependent on
c       nperod and must have the following form0
c
c            nperod
c              = 0 or 2, then n = (2**p)(3**q)(5**r)
c              = 1, then n = (2**p)(3**q)(5**r)-1
c              = 3, then n = (2**p)(3**q)(5**r)+1
c
c            where p, q, and r may be any non-negative integers.  n must
c            be greater than 2  .
c
c     mperod
c       = 0 if a(1) and c(m) are not zero
c       = 1 if a(1) = c(m) = 0
c
c     m
c       the number of unknowns in the i-direction.  m may be any integer
c       greater than 1  .
c
c     a,b,c
c       one-dimensional arrays of length m which specify the
c       coefficients in the linear equations given above.
c
c     idimy
c       the row (or first) dimension of the two-dimensional array y as
c       it appears in the program calling pois.  this parameter is used
c       to specify the variable dimension of y.  idimy must be at least
c       m.
c
c     y
c       a two-dimensional array which specifies the values of the right
c       side of the linear system of equations given above.  y must be
c       dimensioned at least m*n.
c
c     w
c       a one-dimensional array-which must be provided by the user for
c       work space.  the length of w must be at least 6n+5m.
c
c
c     * * * * * * * * * *     on output     * * * * * * * * * *
c
c     y
c       contains the solution x.
c
c     ierror
c       an error flag which indicates invalid input parameters  except
c       for number zero, a solution is not attempted.
c
c       = 0  no error.
c       = 1  m .le. 2  .
c       = 2  n .le. 2 or not of the right form.
c       = 3  idimy .lt. m.
c       = 4  nperod .lt. 0 or nperod .gt. 3  .
c
c     w
c       contains intermediate values that must not be destroyed if
c       pois will be called again with iflag = 1  .
c
c
c
      external        trid       ,tridp
      dimension       y(idimy,1)
      dimension       w(1)       ,b(1)       ,a(1)       ,c(1)
      ierror = 0
      if (m .le. 2) ierror = 1
c     i = ncheck(nperod-2*((2*nperod)/3)+n)
c     if (n.le.2 .or. i.gt.1) ierror = 2
      if (idimy .lt. m) ierror = 3
      if (nperod.lt.0 .or. nperod.gt.3) ierror = 4
      if (ierror .ne. 0) go to 105
      iwdim1 = 6*n+1
      iwdim2 = iwdim1+m
      iwdim3 = iwdim2+m
      iwdim4 = iwdim3+m
      iwdim5 = iwdim4+m
      do 101 i=1,m
         a(i) = -a(i)
         c(i) = -c(i)
         k = iwdim5+i-1
         w(k) = -b(i)+2.
  101 continue
      if (mperod .eq. 0) go to 102
      call poisgn (nperod,n,iflg,m,a,w(iwdim5),c,idimy,y,w(1),
     1             w(iwdim1),w(iwdim2),w(iwdim3),w(iwdim4),trid)
      go to 103
  102 call poisgn (nperod,n,iflg,m,a,w(iwdim5),c,idimy,y,w(1),
     1             w(iwdim1),w(iwdim2),w(iwdim3),w(iwdim4),tridp)
  103 do 104 i=1,m
         a(i) = -a(i)
         c(i) = -c(i)
  104 continue
  105 return
      end
