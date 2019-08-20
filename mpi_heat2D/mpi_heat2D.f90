! ****************************************************************************
! FILE: mpi_heat2D.f
! OTHER FILES: draw_heat.c mpi_heat2D.h
! DESCRIPTION:  
!   HEAT2D Example - Parallelized Fortran Version
!   This example is based on a simplified two-dimensional heat
!   equation domain decomposition.  The initial temperature is computed to be
!   high in the middle of the domain and zero at the boundaries.  The
!   boundaries are held at zero throughout the simulation.  During the
!   time-stepping, an array containing two domains is used; these domains
!   alternate between old data and new data.
! 
!   In this parallelized version, the grid is decomposed by the master
!   process and then distributed by cols to workers processes.  At each
!   time step, workers processes must exchange border data with neighbors,
!   because a grid point's current temperature depends upon it's previous
!   time step value plus the values of the neighboring grid points.  Upon
!   completion of all time steps, the worker processes return their results
!   to the master process.
! 
!   Two data files are produced: an initial data set and a final data set.
!   An X graphic of these two states displays after all calculations have
!   completed.
! 
! AUTHOR: Blaise Barney - adapted from D. Turner's serial C version. Converted
!   to MPI: George L. Gusciora (1/95)
! LAST REVISED: 06/12/13 Blaise Barney
! ****************************************************************************
! Explanation of constants and variables
!   NXPROB                         =  x dimension of problem grid 
!   NYPROB                         =  y dimension of problem grid
!   STEPS                          =  number of time steps
!   MAXWORKER                      =  maximum number of workers tasks
!   MINWORKER                      =  minimum number of workers tasks
!   BEGIN, LTAG, RTAG, DONE        =  message tags
!   NONE                           =  indicates no neighbor
!   CX, CY                         =  used in heat equation
!   u                              =  array for grids
!   taskid,MASTER                  =  taskids 
!   numworkers                     =  number of workers processes
!   numtasks                       =  number of tasks
!   avecol,cols,offset,extra       =  for sending cols of data
!   dest, source                   =  to - from for message send-receive
!   left,right                     =  neighbor tasks
!   msgtype                        =  for message types
!   rc,start,end                   =  misc
!   i,ix,iy,iz,it                  =  loop variables
! ---------------------------------------------------------------------------

      program heat2D
        include 'mpif.h'
        include 'mpi_heat2D.h'

!     Routine for creating the X graph of the wave
        external draw_heat

        integer STEPS,MAXWORKER,MINWORKER,BEGIN,LTAG,RTAG,DONE,NONE,MASTER
        parameter(STEPS=50)
        parameter(MAXWORKER=8)
        parameter(MINWORKER=3)
        parameter(BEGIN=1)
        parameter(LTAG=2)
        parameter(RTAG=3)
        parameter(DONE=4)
        parameter(NONE=0)
        parameter(MASTER=0)
        integer taskid,numtasks,numworkers,avecol,cols,offset,extra, &
          dest,source,left,right,msgtype, &
          rc,start,end,i,ix,iy,iz,it,ierr
        integer status(MPI_STATUS_SIZE)

!     First, find out my taskid and how many tasks are running */
        call MPI_INIT( ierr )
        call MPI_COMM_RANK( MPI_COMM_WORLD, taskid, ierr )
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numtasks, ierr )
        numworkers = numtasks-1

        if (taskid .eq. MASTER) then
!     ****************************** master code *****************************
!     Check if numworkers is within range - quit if not
        if ((numworkers.lt.MINWORKER).or.(numworkers.gt.MAXWORKER)) then
          print *,'MP_PROCS needs to be between', MINWORKER+1,'and', &
            MAXWORKER+1, 'for this exercise'
          print *,'Quitting...'
          call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
        end if

        print *, 'Starting mpi_heat2D with',numworkers,'worker tasks.'

!     Initialize grid 
        print *,'Grid size: X=',NXPROB,' Y=',NYPROB,' Time steps=',STEPS
        print *,'Initializing grid and writing initial.dat file...'
        call inidat
        call prtdat(1)

!     Distribute work to workers.  Must first figure out how many cols to
!     send and what to do with extra cols.  
        avecol=NYPROB/numworkers
        extra=mod(NYPROB,numworkers)
        offset=1
        do i=1, numworkers
          if (i .le. extra) then
            cols=avecol+1
          else
            cols=avecol 
          end if
!       Tell each worker which other workers are its neighbors, since
!       they must exchange data with each other later. 
          if (i .eq. 1) then
            left=NONE
          else  
            left=i-1
          end if
          if (i .eq. numworkers) then
            right=NONE
          else
            right=i+1
          end if 
!       Now send startup information to each worker
          dest = i
        call MPI_SEND( offset, 1, MPI_INTEGER, dest, BEGIN, &
          MPI_COMM_WORLD, ierr )
        call MPI_SEND( cols, 1, MPI_INTEGER, dest, BEGIN, &
          MPI_COMM_WORLD, ierr )
        call MPI_SEND( left, 1, MPI_INTEGER, dest, BEGIN, &
          MPI_COMM_WORLD, ierr )
        call MPI_SEND( right, 1, MPI_INTEGER, dest, BEGIN, &
          MPI_COMM_WORLD, ierr )
        call MPI_SEND( u(1,offset,1), cols * NXPROB, MPI_REAL, dest, &
          BEGIN, MPI_COMM_WORLD, ierr )
          print *,'Sent to=',dest,'offset=',offset,'cols=',cols, &
            'left=',left,'right=',right
          offset = offset + cols
        end do
            
!     Now wait for results from all workers tasks 
        do i=1, numworkers
          source = i
          msgtype = DONE
          call MPI_RECV( offset, 1, MPI_INTEGER, source, &
            msgtype, MPI_COMM_WORLD, status, ierr )
          call MPI_RECV( cols, 1, MPI_INTEGER, source, &
            msgtype, MPI_COMM_WORLD, status, ierr )
          call MPI_RECV( u(1,offset,1), cols * NXPROB, MPI_REAL, &
            source,msgtype,MPI_COMM_WORLD,status,ierr)
        end do

!     Print and show results 
        print *,'Creating final.dat file and generating graph...'
        call prtdat(2)
        call draw_heat()


!     End of master code
        call MPI_FINALIZE(ierr)
        end if

        if (taskid .ne. MASTER)  then
!     ****************************** worker code ******************************
!     Initialize everything - including the borders - to zero 
        do ix=1,NXPROB
          do iy=1,NYPROB
            do iz=1,2
              u(ix,iy,iz) = 0.0
            end do
          end do
        end do

!     Now receive my offset, cols, neighbors and grid partition from master 
        source = MASTER
        msgtype = BEGIN
        call MPI_RECV( offset, 1, MPI_INTEGER, source, &
          msgtype, MPI_COMM_WORLD, status, ierr )
        call MPI_RECV( cols, 1, MPI_INTEGER, source, &
          msgtype, MPI_COMM_WORLD, status, ierr )
        call MPI_RECV( left, 1, MPI_INTEGER, source, &
          msgtype, MPI_COMM_WORLD, status, ierr )
        call MPI_RECV( right, 1, MPI_INTEGER, source, &
          msgtype, MPI_COMM_WORLD, status, ierr )
        call MPI_RECV( u(1,offset,1),cols*NXPROB,MPI_REAL, source, &
          msgtype, MPI_COMM_WORLD, status, ierr )

!     Determine border elements.  Need to consider first and last columns.
!     Obviously, col 1 can't exchange with col 1-1.  Likewise, the last
!     col can't exchange with last+1. 
        start=offset
        end=offset+cols-1
        if (offset .eq. 1) then
          start=2
        end if
        if (offset + cols .gt. NYPROB) then
          end=end-1
        end if
        print *,'task=',taskid,'  start=',start,'  end=',end 

!     Begin doing STEPS iterations.  Must communicate border cols with
!     neighbors.  If I have the first or last grid col, then I only need to 
!     communicate with one neighbor. 
        iz=1
        do it=1, STEPS
          if (left .ne. NONE) then
            call MPI_SEND( u(1,offset,iz), NXPROB, MPI_REAL, left, &
              RTAG, MPI_COMM_WORLD, ierr )

            source = left
            msgtype = LTAG
            call MPI_RECV( u(1,offset-1,iz),NXPROB,MPI_REAL, source, &
              msgtype, MPI_COMM_WORLD, status, ierr )
          end if

          if (right .ne. NONE) then
            call MPI_SEND(u(1,offset+cols-1,iz),NXPROB,MPI_REAL, &              
              right,LTAG,MPI_COMM_WORLD,ierr)

            source = right
            msgtype = RTAG
            call MPI_RECV(u(1,offset+cols,iz),NXPROB,MPI_REAL,source, &
              msgtype, MPI_COMM_WORLD, status, ierr )
          end if
        
!       Now call update to update the value of grid points
          call update(start,end,u(1,1,iz),u(1,1,3-iz))
          
          iz=3-iz
        end do

!     Send my portion of final results back to master 
        call MPI_SEND( offset, 1, MPI_INTEGER, MASTER, DONE, &
          MPI_COMM_WORLD, ierr )
        call MPI_SEND( cols, 1, MPI_INTEGER, MASTER, DONE, &
          MPI_COMM_WORLD, ierr )
        call MPI_SEND( u(1,offset,iz),cols*NXPROB,MPI_REAL,MASTER, &
          DONE, MPI_COMM_WORLD, ierr )

!     End of worker code
        call MPI_FINALIZE(ierr)
        end if

        end

!****************************************************************************
        subroutine update (start, end, u1, u2)
!****************************************************************************
        include 'mpi_heat2D.h'
        integer start, end, ix, iy
        real*4 u1, u2
        dimension u1(NXPROB,NYPROB),u2(NXPROB,NYPROB)

        do iy=start, end
          do ix=2, NXPROB-1
              u2(ix,iy) = u1(ix,iy) &
                + CX * ( u1(ix+1,iy) + u1(ix-1,iy) - 2.0 * u1(ix,iy)) &
                + CY * ( u1(ix,iy+1) + u1(ix,iy-1) - 2.0 * u1(ix,iy))
          end do
        end do
      end

!*****************************************************************************
        subroutine inidat
!*****************************************************************************
        include 'mpi_heat2D.h'
        integer ix,iy

        do ix=0,NXPROB-1
          do iy=0,NYPROB-1
            u(ix+1,iy+1,1) = float(ix*(NXPROB-ix-1) * iy*(NYPROB-iy-1))
          end do
        end do
        end


!**************************************************************************
        subroutine prtdat(i)
!**************************************************************************
        include 'mpi_heat2D.h'
        integer i,ix, iy
        character*11 fname

        if (i .eq. 1) then
          fname ='initial.dat'
        else if (i .eq. 2) then
          fname = 'final.dat'
        end if

        open(21, file=fname, form='formatted')
        do ix=1,NXPROB
          do iy=1,NYPROB
            write(21,'(f8.1,1x,$)')u(ix,iy,1)
          end do
          write(21,'(1x)')
        end do
        close(21)
        end