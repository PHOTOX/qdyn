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

  write(*,*) "====== Qdyn ====="

!--Reading input--!
  call read_input()

!--Initialization--!
  call init()

!--Propagation mode--!
write(*,*) 
write(*,*) "### Propagation ###"
write(*,*)

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

          !jj - field function ready for evaluation
          !print *,elmag_field(time)

        case(2)
          call propag_2d(wf2x(1,:,:),wf2p,theta_v2,kin_p2)
          call update_norm()
          call update_energy_2d(wf2x(1,:,:))

        case(3)
          call propag_3d(wf3x,wf3p,theta_v3,kin_p3)
          call update_norm()
          call update_energy_3d(wf3x, energy(1))
      end select

      !print information
      if (modulo(time,dtwrite) .eq. 0 ) then
        write(*,'(F8.1,a,F14.9,a,F9.7)') time, ' a.u.; E=', energy(1), ' a.u.; norm=', norm
        call printen()
        !TODO: delete x and v1 from printing
        if (print_wf) then
          select case(rank)
          case(1)
           call printwf_1d(1,x,v1)
          case(2)
            call printwf_2d(1,x,y,v2)
          case(3)
            call printwf_3d(wf3x,x,y,z,v3)
          end select
        end if
        if(rank .eq. 1) call printwf_1d(1,x,v1)
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
              do jstate=1,istate-1
               call project_out_1d(wfx(jstate,:),wfx(istate,:))
              end do
            call normalize_1d(wfx(istate,:))
            call update_energy_1d(wfx(istate,:))

          case(2)
            call propag_2d(wf2x(istate,:,:),wf2p,theta_v2,kin_p2)
              do jstate=1,istate-1
               call project_out_2d(wf2x(jstate,:,:),wf2x(istate,:,:))
              end do
            call normalize_2d(wf2x(istate,:,:))
            call update_energy_2d(wf2x(istate,:,:))

          case(3)
            call propag_3d(wf3x,wf3p,theta_v3,kin_p3)
            call normalize_3d(wf3x)
            call update_energy_3d(wf3x, energy(1))
        end select

        !print information
        if (modulo(time,dtwrite) .eq. 0 ) then
          write(*,'(F8.1,a,F14.9,a,F14.9,a)') time, ' a.u.; E=', energy(1), ' a.u.; dE=', energy_diff, ' a.u.'
          call printen_state(istate)
          !TODO: delete x and v1 from printing
          if (print_wf) then
            select case(rank)
            case(1)
             call printwf_1d(istate,x,v1)
            case(2)
              call printwf_2d(istate,x,y,v2)
            !jj modify case(3)
            case(3)
              call printwf_3d(wf3x,x,y,z,v3)
            end select
          end if
        end if
      end do
    end do

end select
 
!--Finalization--!

!TODO: this should be also subroutine
write(*,*) "JOB DONE."

call print_chk()

end program qdyn

