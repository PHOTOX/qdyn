module mod_exactfactor

use mod_vars
use mod_utils

  implicit none
  integer,private     :: file_unit !TODO: used only when plotting for states is used, probabaly remove
! TODO: define the time derivative variable here since it won't be used in any other module. Or the question is allocation

  !TODO: dont forget to save gi ef quantities in the history during initialization

CONTAINS

!=== EXACT FACTORIZATION ===!
subroutine exact_factor_1d(step, gaugedep)
  integer, intent(in)                    :: step
  character(len=2), intent(in)           :: gaugedep
  complex(DP), dimension(nstates,xngrid) :: wf_td_1d


  if (gaugedep=='gi') then ! calculation of gauge-independent EF quantities
    write(*,*) '* GI EF, step = ', step


  elseif (gaugedep=='gd') then ! calculation of gauge-dependent EF quantities
    write(*,*) '* GD EF, step = ', step

    ! calculate time derivative of wf
    wf_td_1d = 0.0d0
    if ((step >= efhistory) .and. (step < nstep)) then
      ! fifth-order central difference formula for steps during propagation
      call wf_time_der_1d(wf_td_1d, '5th-order_cent')
    elseif (step < efhistory) then
      ! third-order forward difference formula for initialization
      call wf_time_der_1d(wf_td_1d, '3th-order_forw')
    elseif (step==nstep) then
      ! fifth-order backward difference formula for the last step
      call wf_time_der_1d(wf_td_1d, '5th-order_back')
    end if

  end if

end subroutine

!=== TIME DERIVATIVES ===!
subroutine wf_time_der_1d(wf_td_1d, order)
  character(len=14), intent(in)                         :: order
  complex(DP), dimension(nstates,xngrid), intent(inout) :: wf_td_1d

  write(*,*) " - Calculating time derivative of wf: ", order
  if (order == '5th-order_cent') then
    write(*,*) " - Central formula for time derivative"
  elseif (order == '3th-order_forw') then
    write(*,*) " - Forward formula for time derivative"
  elseif (order == '5th-order_back') then
    write(*,*) " - Backward formula for time derivative"
  end if

end subroutine

!=== SAVE HISTORY ===!
subroutine ef_savehistory_1d()
! save wave function and phase history

  do j=1, nstates ! looping through states
    do i=1, efhistory-1 ! looping through history
      wfx_history(efhistory+1-i,j,:) = wfx_history(efhistory-i,j,:) ! shifting previous wfx
    end do
    wfx_history(1,j,:) = wfx(j,:) ! saving current wfx
  end do

end subroutine

!=== INITIALIZATION ===!
subroutine init_ef
  !TODO: finish initialization of ef
  write(*,*) "---------"
  write(*,*) "Exact factorization initialized"

  !TODO: dont forget to save gi ef quantities in the history during initialization
  write(*,*) "Save GI EF quantities at the beginning"

  select case(rank)
  case(1)
    call exact_factor_1d(0, 'gi')
  end select
  call print_ef('gi')


end subroutine

!=== PRINTING ===!
subroutine print_ef(gaugedep)
  character(len=2), intent(in)  :: gaugedep


  write(*,*) "First I need to create the files in init()!"

  if (gaugedep=='gi') then
    write(*,*) 'write GI'
  elseif (gaugedep=='gd') then
    write(*,*) 'write GD'
  end if

end subroutine

end module
