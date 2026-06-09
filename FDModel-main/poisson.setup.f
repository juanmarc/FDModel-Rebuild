C*****7***************************************************************72

      subroutine poissonSetup(mm1,a,b,c,dx,dy)

C*****7***************************************************************72
C This routine will set up the arrays representing the coefficients
C in the arrays A, B, and C determined by the following representation
C of the Laplaian of the function PSI in finite difference form:
C
C         A(I)*PSI(I-1,J) + B(I)*PSI(I,J) + C(I)*PSI(I+1,J)
C         + PSI(I,J-1) -2*PSI(I,J) + PSI(I,J+1) = ZETA(I,J)
C
C Where the normal finite difference representation of the Laplacian
C of a function has been multiplied by (dy)^2.  The reason for this is
C done by its own subroutine is because the values of A, B, and C will
C not change through the course of the run.  So, to save computer time,
C we define it once and just pass it to the function needing it.
C
C*****7***************************************************************72
C Variables used by poissionSetup:
C
C*****7***************************************************************72
C** Variable declarations

      implicit none
      
      integer mm1,i
      double precision a(mm1),b(mm1),c(mm1)
      double precision dx,dy,dxsq,dysq,temp
      
C*****7***************************************************************72
C** Constant declarations

      double precision two
      parameter (two = 2.0D0)
      
C*****7***************************************************************72
C** Determine dxsq, dysq, and dysq/dxsq

      dxsq = dx*dx
      dysq = dy*dy
      temp = dysq/dxsq
      
C*****7***************************************************************72
C** Compute the contents of the arrays

      do 800 i=1,mm1
        a(i) = temp
        b(i) = -two*temp
        c(i) = temp
 800  continue
      
C*****7***************************************************************72

      return
      end