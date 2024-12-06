module mod_exactfactor

use mod_vars
use mod_utils

  implicit none
  integer,private     :: file_unit !TODO: used only when plotting for states is used, probabaly remove
! TODO: define the time derivative variable here since it won't be used in any other module. Or the question is allocation

CONTAINS

!=== EXACT FACTORIZATION ===!
subroutine exact_factor_1d(step)
  integer, intent(in)       :: step

  write(*,*) "Calculate final EF quantities that are then printed"

  if (step .ge. efhistory) then
    call wf_time_der_1d()
  end if

    
end subroutine

!=== TIME DERIVATIVES ===!
subroutine wf_time_der_1d() !TODO: this is currently braket functions for utils
!  complex(DP)               :: wf_time_der_1d=0.0d0
  integer                   :: i

  write(*,*) "Time derivative would be calculated"



end subroutine

!=== Save history ===!
subroutine ef_savehistory_1d()
! save wave function and phase history
  integer                       :: i, j

  ! saving wf history
  write(*,*) "Saving wf history!"

  do j=1, nstates ! looping through states
    do i=1, efhistory-1 ! looping through history
      wfx_history(i+1,j,:) = wfx_history(i,j,:) ! shifting previous wfx
    end do
    wfx_history(1,j,:) = wfx(j,:) ! saving current wfx
  end do

  write(*,*) "Saving phase history?"
  !TODO: save phase history, I will need to calculate the phase for this


end subroutine

!=== PRINTING ===!
subroutine print_ef()

  write(*,*) "First I need to create the files in init()!"
  !write(101,'(F8.1,5(F16.9))') time, energy, energy_diff, norm
    
end subroutine

end module
