program qdyn
   use mod_vars
   use FFTW3
   use mod_init
   use mod_propag
   use mod_exactfactor

   ! -------------------------------------------------------------------!
   !                               / Qdyn \                             !
   ! Program for N-dimensional numerical quantum propagation on a grid. !
   ! Authors: Jiri Janos, Jiri Suchan, Petr Slavicek (2025)             !
   ! Copyright PHOTOX group                                             !
   ! -------------------------------------------------------------------!

   implicit none
   integer :: n, istate, jstate

   write(*, *) "====== Qdyn ====="

   !--Reading input--!
   call read_input()

   !--Initialization--!
   call init()

   !--Propagation mode--!
   write(*, *)
   write(*, *) "### Propagation ###"
   write(*, *)

   select case(run)
      ! Real time propagation
   case(0)
      do n = 1, nstep
         time = n * dt

         select case(rank)
         case(1)
            call propag_rt_1d()
            if (exact_factor) call ef_savehistory_1d()
         case(2)
            call propag_rt_2d(wf2x(1, :, :))
         case(3)
            call propag_rt_3d(wf3x(1, :, :, :))
         end select

         call update_norm()

         !print information
         if ((modulo(time, dtwrite)==0).or.(n==nstep)) then
            !update and print energies
            select case(rank)
            case(1)
               call update_total_energy_1d()
            case(2)
               call update_energy_2d(wf2x(1, :, :))
            case(3)
               call update_energy_3d(wf3x(1, :, :, :))
            end select
            call printen()

            write(*, '(a,I10,a,F8.1,a,F14.9,a,F9.7)') 'Step: ', n, ' ; ', time, ' a.u.; E=', energy(1), &
                  ' a.u.; norm=', norm

            !if multistate problem, transform wf to adiabatic basis and print populations of ad. and diab. states
            if (nstates>1) then
               call wf_adiab_trans()
               call print_pop()
            end if

            !print wave function
            if (print_wf) then
               do istate = 1, nstates
                  select case(rank)
                  case(1)
                     call printwf_1d(istate)
                     if (nstates>1) call printwf_ad_1d(istate)
                  case(2)
                     call printwf_2d(istate)
                  case(3)
                     call printwf_3d(istate)
                  end select
               end do
            end if

            !print field
            if (field_coupling) call print_field()

         end if

         ! calculate and print exact factorization
         if (exact_factor) then
            !calculate GI exact factorization quantites and print them
            ! this is done during every printing step and then at the end of the simualtion
            if ((modulo(time, dtwrite)==0).or.(n==nstep)) then
               select case(rank)
               case(1)
                  call exact_factor_1d(n)
                  call print_ef_1d()
               end select
            end if

            !calculate GD-TDPES and print it
            !time derivative of wf is necessary for GD-TDPES which can be calculated only two steps later using central difference formula
            !the calcualtion comes two steps after the beginning because this cannot be done during initialization and also for the
            !last step of the dynamics
            if ((n==efhistory / 2).or.(modulo((n - efhistory / 2) * dt, dtwrite)==0).or.(n==nstep)) then
               select case(rank)
               case(1)
                  call exact_factor_gd_tdpes_1d(n)
                  call print_ef_gd_tdpes_1d(n)
               end select
            end if
         end if

      end do

   case(1)
      ! Imaginary time propagation
      do istate = 1, nstates
         write(*, '(a,I2)') "* Optimizing state ", istate
         do n = 1, nstep
            time = n * dt

            select case(rank)
            case(1)
               call propag_it_1d(wfx(istate, :))
               do jstate = 1, istate - 1
                  call project_out_1d(wfx(jstate, :), wfx(istate, :))
               end do
               call normalize_1d(wfx(istate, :))

            case(2)
               call propag_it_2d(wf2x(istate, :, :))
               do jstate = 1, istate - 1
                  call project_out_2d(wf2x(jstate, :, :), wf2x(istate, :, :))
               end do
               call normalize_2d(wf2x(istate, :, :))

            case(3)
               call propag_it_3d(wf3x(istate, :, :, :))
               do jstate = 1, istate - 1
                  call project_out_3d(wf3x(jstate, :, :, :), wf3x(istate, :, :, :))
               end do
               call normalize_3d(wf3x(istate, :, :, :))
            end select

            !print information
            if ((modulo(time, dtwrite)==0).or.(n==nstep)) then
               select case(rank)
               case(1)
                  call update_energy_1d(wfx(istate, :))
               case(2)
                  call update_energy_2d(wf2x(istate, :, :))
               case(3)
                  call update_energy_3d(wf3x(istate, :, :, :))
               end select
               call printen_state(istate)

               write(*, '(a,I10,a,F8.1,a,F14.9,a,F14.9,a)') 'Step: ', n, ' ; ', time, ' a.u.; E=', energy(1), &
                     ' a.u.; dE=', energy_diff, ' a.u.'
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

      call print_chk(nstates) ! this prints last IT wave function which can be than loaded to next calculation

   end select

   !--Finalization--!
   write(*, *) "JOB DONE."

end program qdyn
