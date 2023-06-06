program qdyn
  use mod_vars
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

!--Initialization--!

  call read_input()
  call init()


!--Propagation mode--!

do n=1, nstep
  time = n*dt

  select case(rank)
  case(1)
      call propag_1d(wfx,wfp,theta_v1,kin_p1) 
      !jj - I should already project out in the init part
      if (run.eq.1 .and. nstates.eq.2) call project_out_1d(wfxgs,wfx)
      call normalize_1d(wfx)
      call update_energy_1d(wfx)

  case(2)
      call propag_2d(wf2x,wf2p,theta_v2,kin_p2)
      call normalize_2d(wf2x)
      call update_energy_2d(wf2x, energy(1))

  case(3)
      call propag_3d(wf3x,wf3p,theta_v3,kin_p3)
      call normalize_3d(wf3x)
      call update_energy_3d(wf3x, energy(1))
  end select
 
!Printing
  if (modulo(time,dtwrite) .eq. 0 ) then
    write(*,'(F8.1,a,F14.9,a)') time, ' a.u.; E=', energy(1), ' a.u.'
    call printen()
    !TODO: change to case
    !TODO: delete x and v1 from printing
    if(rank .eq. 1) call printwf_1d(wfx,x,v1)
    !TODO: slower writing of wf.. for 2 and 3 dim it is too much data
!    if(rank .eq. 2) call printwf_2d(wf2x,x,y,v2)
!    if(rank .eq. 3) call printwf_3d(wf3x,x,y,z,v3)
  end if


end do

!--Finalization--!

!TODO: this should be also subroutine
write(*,*) "JOB DONE."

call print_chk()

end program qdyn

