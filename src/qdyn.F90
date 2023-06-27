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
          call propag_rt_1d(wfx(1,:)) 
          call update_norm()

        case(2)
          call propag_rt_2d(wf2x(1,:,:))
          call update_norm()

        case(3)
          call propag_rt_3d(wf3x(1,:,:,:))
          call update_norm()
      end select

      !print information
      if ((modulo(time,dtwrite).eq.0).or.(n.eq.nstep)) then
        select case(rank)
        case(1)
          call update_energy_1d(wfx(1,:))
        case(2)
          call update_energy_2d(wf2x(1,:,:))
        case(3)
          call update_energy_3d(wf3x(1,:,:,:))
        end select
        call printen()

        write(*,'(F8.1,a,F14.9,a,F9.7)') time, ' a.u.; E=', energy(1), ' a.u.; norm=', norm
        !TODO: delete x and v1 from printing
        if (print_wf) then
          select case(rank)
          case(1)
           call printwf_1d(1)
          case(2)
            call printwf_2d(1)
          case(3)
            call printwf_3d(1)
          end select
        end if

        if (use_field) call print_field()

      end if
    end do

  case(1)
  ! Imaginary time propagation
    do istate=1, nstates
      write(*,'(a,I2)') "* Optimizing state ",istate
      do n=1, nstep
        time = n*dt

        select case(rank)
          case(1)
            call propag_it_1d(wfx(istate,:)) 
              do jstate=1,istate-1
               call project_out_1d(wfx(jstate,:),wfx(istate,:))
              end do
            call normalize_1d(wfx(istate,:))

          case(2)
            call propag_it_2d(wf2x(istate,:,:))
              do jstate=1,istate-1
               call project_out_2d(wf2x(jstate,:,:),wf2x(istate,:,:))
              end do
            call normalize_2d(wf2x(istate,:,:))

          case(3)
            call propag_it_3d(wf3x(istate,:,:,:))
              do jstate=1,istate-1
               call project_out_3d(wf3x(jstate,:,:,:),wf3x(istate,:,:,:))
              end do
            call normalize_3d(wf3x(istate,:,:,:))
        end select

        !print information
        if ((modulo(time,dtwrite).eq.0).or.(n.eq.nstep)) then
          select case(rank)
          case(1)
            call update_energy_1d(wfx(istate,:))
          case(2)
            call update_energy_2d(wf2x(istate,:,:))
          case(3)
            call update_energy_3d(wf3x(istate,:,:,:))
          end select
          call printen_state(istate)

          write(*,'(F8.1,a,F14.9,a,F14.9,a)') time, ' a.u.; E=', energy(1), ' a.u.; dE=', energy_diff, ' a.u.'
          if (print_wf) then
            select case(rank)
            case(1)
             call printwf_1d(istate)
            case(2)
              call printwf_2d(istate)
            case(3)
              call printwf_3d(istate)
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

