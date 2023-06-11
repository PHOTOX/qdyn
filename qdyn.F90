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
   ! TODO: istate variable added
   integer     :: n

!--Initialization--!

  call read_input()
  call init()


!--Propagation mode--!

select case(run)
  ! Real time propagation
  case(0)
    do n=1, nstep
      time = n*dt

      select case(rank)
        case(1)
          call propag_1d(wfx(1,:),wfp,theta_v1,kin_p1) 
          call update_norm()
          call update_energy_1d(wfx(1,:))

        case(2)
          call propag_2d(wf2x,wf2p,theta_v2,kin_p2)
          call update_norm()
          call update_energy_2d(wf2x, energy(1))

        case(3)
          call propag_3d(wf3x,wf3p,theta_v3,kin_p3)
          call update_norm()
          call update_energy_3d(wf3x, energy(1))
      end select

      !print information
      if (modulo(time,dtwrite) .eq. 0 ) then
        write(*,'(F8.1,a,F14.9,a,F9.7)') time, ' a.u.; E=', energy(1), ' a.u.; norm=', norm
        call printen()
        !TODO: change to case
        !TODO: delete x and v1 from printing
        if(rank .eq. 1) call printwf_1d(1,x,v1)
        !TODO: slower writing of wf.. for 2 and 3 dim it is too much data
      !    if(rank .eq. 2) call printwf_2d(wf2x,x,y,v2)
      !    if(rank .eq. 3) call printwf_3d(wf3x,x,y,z,v3)
      end if
    end do

  case(1)
  ! Imaginary time propagation
    do istate=1, nstates
      write(*,'(a,I2)') "* Optimizing state ",istate
      do n=1, nstep
        time = n*dt

        !select dimension to propagate
        select case(rank)
          case(1)
            call propag_1d(wfx(istate,:),wfp,theta_v1,kin_p1) 
            !jj - I should already project out in the init part
!            if (istate.ge.2) then
              do jstate=1,istate-1
              write(*,*) "Projecting out state:", jstate
               call project_out_1d(wfx(jstate,:),wfx(istate,:))
              end do
!            end if
            call normalize_1d(wfx(istate,:))
            call update_energy_1d(wfx(istate,:))

          case(2)
            call propag_2d(wf2x,wf2p,theta_v2,kin_p2)
            call normalize_2d(wf2x)
            call update_energy_2d(wf2x, energy(1))

          case(3)
            call propag_3d(wf3x,wf3p,theta_v3,kin_p3)
            call normalize_3d(wf3x)
            call update_energy_3d(wf3x, energy(1))
        end select

        !print information
        if (modulo(time,dtwrite) .eq. 0 ) then
          write(*,'(F8.1,a,F14.9,a,F14.9,a)') time, ' a.u.; E=', energy(1), ' a.u.; dE=', energy_diff, ' a.u.'
          call printen_state(istate)
          !TODO: change to case
          !TODO: delete x and v1 from printing
          if(rank .eq. 1) call printwf_1d(istate,x,v1)
          !TODO: slower writing of wf.. for 2 and 3 dim it is too much data
        !    if(rank .eq. 2) call printwf_2d(wf2x,x,y,v2)
        !    if(rank .eq. 3) call printwf_3d(wf3x,x,y,z,v3)
        end if
      end do
    end do

end select
 
!--Finalization--!

!TODO: this should be also subroutine
write(*,*) "JOB DONE."

call print_chk()

end program qdyn

