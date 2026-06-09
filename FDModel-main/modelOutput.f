C*****7***************************************************************72

      subroutine modelOutput(n,m,zeta,psi)

C*****7***************************************************************72
C This routine will print out the values of the vorticity array and
C the streamfunction array each time it is called
C*****7***************************************************************72
C Variables used by modelOutput:
C
C*****7***************************************************************72
C** Variable declarations:

      implicit none

      integer n,m
      double precision zeta(m,n),psi(m,n)

C*****7***************************************************************72
C** Write out the arrays

      write(25,*) zeta
      write(25,*) psi

C*****7***************************************************************72

      return
      end