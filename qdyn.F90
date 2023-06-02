program qdyn
  use FFTW3
  use mod_init
  use mod_propag

! -------------------------------------------------------------------!
!                               / Qdyn \                             !
! Program for N-dimensional numerical quantum propagation on a grid. !
! Authors: Jiri Janos, Jiri Suchan, Petr Slavicek (2023)             !
! -------------------------------------------------------------------!

   implicit none
   integer     :: n
   real(DP)    :: time

!--Initialization--!

  call init()

!--Propagation mode--!


do n=1, nstep
  time = n*dt
  if (rank .eq. 1) then
      call propag_1d(wfx,wfp,ngrid,theta_v1,kin_p1,dt) 
      call normalize1d(wfx,ngrid,dx)
  end if

  if (rank .eq. 2) then
      call propag_2d(wf2x,wf2p,ngrid,theta_v2,kin_p2,dt)
      call normalize2d(wf2x,ngrid,dx)
  endif

  if (rank .eq. 3) then
      call propag_3d(wf3x,wf3p,ngrid,theta_v3,kin_p3,dt)
      call normalize3d(wf3x,ngrid,dx)
  endif
 
!Printing
      if (modulo(time,dtwrite) .eq. 0 ) then
        write(*,'(F8.1,a)') time
        if(rank .eq. 1) call printwf1d(wfx,x,v1)
        if(rank .eq. 2) call printwf2d(wf2x,x,y,v2)
        if(rank .eq. 3) call printwf3d(wf3x,x,y,z,v3)
      end if


end do

!--Finalization--!

! this should be also subroutine
write(*,*) "JOB DONE."

  open(666,file='wf.chk', action='WRITE', iostat=iost)
  write(666,*) "#QDYN checkpoint file for reloading WF to program"
  write(666,*) "#Rank:",rank,"Pot:",pot,"Ngrid:",ngrid
  if(rank .eq. 1) write(666,*) wfx 
  if(rank .eq. 2) write(666,*) wf2x
  if(rank .eq. 3) write(666,*) wf3x
  close(666)

end program qdyn

