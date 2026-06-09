C*****7***************************************************************72

      program fd_model

C*****7***************************************************************72
C  This is the main code for the model.  All actions are initiated from
C  the code sequence.  For documentation on the variables used, see the 
C  file fd.common.f
C*****7***************************************************************72
C Variables not covered in fd.common
C             K = the do loop index to iterate through each output
C                 increment
C  hr, min, sec = the clock time (in model time, not real time)
C                 commensurate with the time step
C*****7***************************************************************72

      implicit none
      include 'fd.common'
      integer k,hr,min
      real sec

C*****7***************************************************************72
C**  Read in the specifiable model parameters.  Set up the coefficient
C**  arrays needed for the poisson solver.  And, set up the grid.
C**  Establish the initial vorticity field

      write(*,*) 'call setup'
      call setup(restart,dt,lastGood,nu,state,outputInc,n,m,nm1,
     &            mm1,idimy,maxTimeCount,order,scndOrder,frthOrder,
     & 			 dx,dy,g,f0,beta,omega)
	  write(*,*) 'poissonSetup'
      call poissonSetup(mm1,a,b,c,dx,dy)
	  write(*,*) 'gridInitialize'
      call gridInitialize(n,m,dx,dy,radius,theta)

C*****7***************************************************************72
C**  If the run is a restart, then open the appropriate files as 
C**  'old' files, read-in the 'initial" vorticity field from the
C**  output file, and compute the time count based on the outputInc and
C**  the lastGood hour (restartInitialize does these last two parts)
 
      if (restart) then
         write(*,*) 'restart'
C        open(unit=25,file='fd.data',status='old')
C        call restartInitialize(n,m,timeCount,outputInc,zeta,psi,
C     &                         lastGood,dt)

C*****7***************************************************************72
C**  Otherwise, open the files as 'new' files, define the 
C**  initial vorticity field, print out the initial fields, and set the
C**  timeCount equal to 0.

      else
        open(unit=25,file='fd.data',status='unknown')
        write(*,*) 'basicState'
        call basicState(n,m,nm1,mm1,radius,theta,zeta,psi,
     &                  a,b,c,wdim,dx,dy,idimy)
        write(*,*) 'modelOutput'
        call modelOutput(n,m,zeta,psi)
      endif
      
      
      write(*,*) 'done'

      close (unit=25)
      
C*****7***************************************************************72

      stop
      end program fd_model
