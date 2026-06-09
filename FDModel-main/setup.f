C*****7***************************************************************72

	  subroutine setup(restart,dt,lastGood,nu,state,outputInc,n,m,nm1,
     & 				   mm1,idimy,maxTimeCount,order,scndOrder,frthOrder,
     & 				   dx,dy,g,f0,beta,omega)

C*****7***************************************************************72
C  This routine will read in the appropriate run parameters from the 
C  file 'fd.parameters' and thenprint them out to the file indicated
C  by the variable outFileParameters.  In doing so, the size of the 
C  domain, the type of vortex, and whether or not the run is a restart
C  is established
C*****7***************************************************************72
C** Declarations

      implicit none
      
      double precision nu,dx,dy,dt,g,f0,beta,omega
      
      integer maxTimeCount,outputInc,month,day,yr,state
      integer order,n,m,nm1,mm1,idimy,openStatus
      
      real lastGood
      
      logical restart,scndOrder,frthOrder
      
      character*24 outFileParameters
      character*24 runName
      character*3 restartReply
      character*3 outFileStatus
      
C*****7***************************************************************72
C** Read in the parameters

      open(unit=1,file='fd.parameters')
      
      read(1,*) outFileParameters
      read(1,*) runName
      read(1,*) month, day, yr
      read(1,*) dx
      read(1,*) dy
      read(1,*) n
      read(1,*) m
      read(1,*) dt
      read(1,*) g
      read(1,*) f0
      read(1,*) beta
      read(1,*) omega
      read(1,*) nu
      read(1,*) restartReply
      read(1,*) lastGood
      read(1,*) maxTimeCount
      read(1,*) outputInc
      read(1,*) order
      read(1,*) state
      
      close(unit=1)
      
C*****7***************************************************************72
C** Determine some additional parameters

      nm1 = n-1
      mm1 = m-1
      idimy = m-1

C*****7***************************************************************72
C** Based on the contents of restartReply, determine whether or not
C** this run is a restart of a previous run.  restartReply must be
C** of the form 'YES', 'yes', 'NO', or 'no'

      if ((restartReply .eq. 'NO') .or. (restartReply .eq. 'no')) then
        restart = .false.
        outFileStatus = 'NEW'
      else
        restart = .true.
        outFileStatus = 'OLD'
      endif
      
C*****7***************************************************************72
C** Based on the value of order, determine if this is a 2nd or a
C** 4th order run

      if (order .eq. 0) then
        scndOrder = .true.
        frthOrder = .false.
      else
        scndOrder = .false.
        frthOrder = .true.
      endif
      
C*****7***************************************************************72
C** Write out the appropriate parameters to the file indicated by
C** the variable outFileParameters

      open(unit=50,file=outFileParameters, status='unknown') 
      write(50,*) 'NAME OF RUN:'
      write(50,*) runName
      write(50,*) 'DATE OF RUN (mo day yr):'
      write(50,220) month,' ',day,' ',yr
      write(50,*) 'HORIZ. GRID SPACING (X DIRECTION) IN METERS (dx):'
      write(50,250) dx
      write(50,*) 'HORIZ. GRID SPACING (Y DIRECTION) IN METERS (dy):'
      write(50,250) dy
      write(50,*) 'NUMBER OF HORIZONTAL GRID POINTS (X DIRECTION):'
      write(50,200) n
      write(50,*) 'NUMBER OF HORIZONTAL GRID POINTS (Y DIRECTION):'
      write(50,200) m
      write(50,*) 'TIME STEP IN SECONDS (dt):'
      write(50,250) dt
      write(50,*) 'GRAVITIONAL ACCELERATION (m/s^2) (g):'
      write(50,250) g
      write(50,*) 'CORIOLIS PARAMETER (1/s) (f0):'
      write(50,290) f0
      write(50,*) 'BETA (beta):'
      write(50,290) beta
      write(50,*) 'VISCOSITY (m^2/s) (nu):'
      write(50,280) nu
      write(50,*) 'RESTART OF PREVIOUS RUN? (Y/N):'
      write(50,*) restartReply
      write(50,*) 'IF RESTART, LAST HOUR OF GOOD DATA:'
      write(50,250) lastGood
      write(50,*) 'NUMBER OF TIME COUNTS (maxTimeCount):'
      write(50,200) maxTimeCount
      write(50,*) 'ITERATIONS PER OUTPUT INCREMENT (outputInc):'
      write(50,200) outputInc
      write(50,*) 'NUMBER OF OUTPUT INCREMENTS:'
      write(50,200) maxTimeCount/outputInc
      write(50,*) 'OUTPUT INCREMENT IN HOURS:'
      write(50,250) outputInc*dt/3600
      write(50,*) '2ND ORDER (0) OR 4TH ORDER (1) (order):'
      write(50,200) order
      write(50,*) 'RING(0), MONOPOLE(1), OR CIRCUMOLAR(2) (state):'
      write(50,200) state
      
 250  format(f13.2)
 280  format(f13.4)
 290  format(f13.8)
 200  format(i13)
 220  format(i3,a1,i3,a1,i3)
 
 	  close(unit=50)
 	  
C*****7***************************************************************72

      return
      end