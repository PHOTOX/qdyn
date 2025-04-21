module mod_exactfactor

   use mod_vars
   use mod_utils

   implicit none
   complex(DP), parameter :: imag_unit = cmplx(0.0d0, 1.0d0)
   real(DP), parameter :: ef_zero = 0.00000001D0
   integer, private :: file_unit
   character(len = 50), private :: file_name

CONTAINS

   !=== EXACT FACTORIZATION ===!
   ! calculation of all Exact Factorization quantities but GD-TDPES
   subroutine exact_factor_1d(step)
      integer, intent(in) :: step
      complex(DP), dimension(nstates, xngrid) :: grad_C_ef_1d

      ! nuclear density
      nucdens_1d = calc_nucdens_1d(wfx)

      ! set the gauge
      call set_gauge_1d(wfx, nucdens_1d, phase_1d, grad_phase_1d)

      ! calculate vector potential
      call calc_vecpot_1d(wfx, nucdens_1d, vecpot_1d)

      ! calculate electronc coefficients
      call calc_elcoef_1d(wfx, nucdens_1d, phase_1d, C_ef_1d, grad_C_ef_1d)

      ! calculate gi-tdpes
      call calc_gitdpes_1d(wfx, nucdens_1d, C_ef_1d, grad_C_ef_1d, vecpot_1d, gi_tdpes_1d)

   end subroutine

   ! todo: finish GD TDPES
   ! calculation of GD-TDPES only
   subroutine exact_factor_gd_tdpes_1d(step)
      integer, intent(in) :: step
      integer :: current_time_index ! center of the history arrays to get the actual wf for which we calculate GD-TDPES
      integer :: hindex, state
      character(len = 14) :: deriv_order
      complex(DP), dimension(nstates, xngrid) :: td_wf ! wave function time derivative
      complex(DP), dimension(xngrid) :: sum_over_states ! the summation of derivatives over all states appearing in GD-TDPES
      real(DP), dimension(efhistory, xngrid) :: nucdens_hist, phase_hist
      real(DP), dimension(xngrid) :: td_phase, grad_phase ! nuclear phase time derivative, history of phase gradient is not
      ! necessary but we need some variable to go to the function set_gauge_1d()

      ! assess the history available and the numerical formula used for the derivative and set the current index in the history
      if ((step >= efhistory) .and. (step < nstep)) then
         ! fifth-order central difference formula for steps during propagation
         deriv_order = '5th-order_cent'
         ! center of the history arrays to get the actual wf for which we calculate GD-TDPES,
         ! if efhistory=5 -> current_time_index=3 (center of the history array)
         current_time_index = efhistory - efhistory / 2 ! maybe unnecesarly difficult but works for now
      elseif (step < efhistory) then
         ! third-order forward difference formula for initialization
         deriv_order = '3th-order_forw'
         current_time_index = 1 ! here the current_time_index is the first because we use forward formula
      elseif (step==nstep) then
         ! fifth-order backward difference formula for the last step
         deriv_order = '5th-order_back'
         current_time_index = efhistory ! here the current_time_index is the last because we use backward formula
      end if

      ! calculate nuclear density
      nucdens_1d = calc_nucdens_1d(wfx_history(current_time_index, :, :))

      ! calculate nuclear phase history
      do hindex = 1, efhistory
         nucdens_hist(hindex, :) = calc_nucdens_1d(wfx_history(hindex, :, :))
         call set_gauge_1d(wfx_history(hindex, :, :), nucdens_hist(hindex, :), phase_hist(hindex, :), grad_phase)
      end do

      ! calculate nuclear phase time derivative
      call phase_time_der_1d(td_phase, phase_hist, deriv_order)

      ! calculate time derivative of wf
      call wf_time_der_1d(td_wf, sum_over_states, deriv_order)

      ! construct GD-TDPES
      gd_tdpes_1d = 0.0d0

      do state = 1, nstates
         do i = 1, xngrid
            if (nucdens_1d(i) > ef_zero) then
               !todo: check why the sign is opposite
               !todo: seems there are small problems at the first and last step
               gd_tdpes_1d(i) = gd_tdpes_1d(i) &
                     + (aimag(conjg(wfx_history(current_time_index, state, i)) * td_wf(state, i)) &
                           - abs(wfx_history(current_time_index, state, i))**2 * td_phase(i)) / nucdens_1d(i)
            end if
         end do
      end do

   end subroutine

   !=== SET GAUGE ===!
   subroutine set_gauge_1d(wf, nucdens, S, gradS)
      real(DP), dimension(xngrid), intent(in) :: nucdens
      complex(DP), dimension(nstates, xngrid), intent(in) :: wf
      real(DP), dimension(xngrid), intent(inout) :: S, gradS

      if (ef_gauge == 'S0') then
         S = 0.0d0
         gradS = 0.0d0
      elseif (ef_gauge == 'A0') then
         gradS = vecpotS0_1d(wf, nucdens)
         S = 0.0d0
         S(1) = gradS(1) * dx
         do i = 2, xngrid
            S(i) = S(i - 1) + gradS(i) * dx
         end do
      end if

   end subroutine

   !=== NUCLEAR DENSITY ===!
   function calc_nucdens_1d(wfx)
      complex(DP), dimension(nstates, xngrid), intent(in) :: wfx
      real(DP), dimension(xngrid) :: calc_nucdens_1d

      calc_nucdens_1d = 0.0d0
      do i = 1, nstates
         calc_nucdens_1d = calc_nucdens_1d + abs(wfx(i, :))**2
      end do
   end function

   !=== VECTOR POTENTIAL ===!
   subroutine calc_vecpot_1d(wf, nucdens, A)
      real(DP), dimension(xngrid), intent(in) :: nucdens
      complex(DP), dimension(nstates, xngrid), intent(in) :: wf
      real(DP), dimension(xngrid), intent(inout) :: A

      if (ef_gauge == 'S0') then
         A = vecpotS0_1d(wf, nucdens)
      elseif (ef_gauge == 'A0') then
         A = 0.0d0
      end if
   end subroutine

   function vecpotS0_1d(wf, nucdens)
      real(DP), dimension(xngrid), intent(in) :: nucdens
      complex(DP), dimension(nstates, xngrid), intent(in) :: wf
      complex(DP), dimension(nstates, xngrid) :: tmp_wf ! for FFT
      integer :: istate
      real(DP), dimension(xngrid) :: vecpotS0_1d

      vecpotS0_1d = 0.0d0

      ! calculating the expectation value with i*px operator
      do istate = 1, nstates ! states to be propagated
         ! iFFT -> k
         call dfftw_plan_dft_1d(plan_backward, xngrid, wf(istate, :), wfp, FFTW_BACKWARD, FFTW_ESTIMATE)
         call dfftw_execute_dft(plan_backward, wf(istate, :), wfp)
         call dfftw_destroy_plan(plan_backward)

         wfp = wfp / dsqrt(real(xngrid, kind = DP))

         do i = 1, xngrid
            wfp(i) = - cmplx(0.0D0, 1.0D0) * px(i) * wfp(i) !hbar d\dx psi = - i*px*psi
         end do

         ! FFT -> x
         call dfftw_plan_dft_1d(plan_forward, xngrid, wfp, tmp_wf(istate, :), FFTW_FORWARD, FFTW_ESTIMATE)
         call dfftw_execute_dft(plan_forward, wfp, tmp_wf(istate, :))
         call dfftw_destroy_plan(plan_forward)

         tmp_wf(istate, :) = tmp_wf(istate, :) / dsqrt(real(xngrid, kind = DP))

         vecpotS0_1d = vecpotS0_1d + aimag(conjg(wf(istate, :)) * tmp_wf(istate, :))
      end do

      ! checking for division by zero and calculating vecpot
      do i = 1, xngrid
         if (nucdens(i) > ef_zero) then
            vecpotS0_1d(i) = vecpotS0_1d(i) / nucdens(i)
         else
            vecpotS0_1d(i) = 0.0D0
         end if
      end do

   end function

   !=== ELECTRONIC COEFFICIENTS ===!
   subroutine calc_elcoef_1d(wf, nucdens, S, C, dC)
      real(DP), dimension(xngrid), intent(in) :: nucdens, S
      complex(DP), dimension(nstates, xngrid), intent(in) :: wf
      complex(DP), dimension(nstates, xngrid), intent(inout) :: C, dC
      integer :: istate

      C = 0.0d0
      dC = 0.0d0

      ! Calculating electronic coefficients and the gradient of electronic coefficients.
      ! However, I first calculate C with lower threshold for zero to make it smooth even at the edges which is necessary for
      ! smooth derivative. I calculate the derivative in this extended interval. Then, I enforce the condition of C and dC being
      ! zero out of the interval where nuclear density is bigger than zero_ef.
      do istate = 1, nstates
         do i = 1, xngrid
            if (nucdens(i) > ef_zero * 0.1d0) then ! lower threshold here, C is in extended range
               C(istate, i) = wf(istate, i) / dsqrt(nucdens(i)) * exp(-imag_unit * S(i))
            end if
         end do
         do i = 1, xngrid
            dC(istate, :) = dfdx_1d(C(istate, :), xngrid, dx, i) ! performing derivative in the extended range
         end do
      end do

      ! now enforcing the true zero_ef (the derivative will now be smooth in this range, even at edges)
      do istate = 1, nstates
         do i = 1, xngrid
            if (nucdens(i) < ef_zero) then ! lower threshold here, C is in extended range
               C(istate, i) = 0.0d0
               dC(istate, i) = 0.0d0
            end if
         end do
      end do

   end subroutine

   function dfdx_1d(f, dim, dx, index)
      ! Numeric derivative of function f using fifth-order central formula except for edges where respective smaller-order
      ! formulas are used.
      integer, intent(in) :: dim ! dimension of the array
      complex(DP), intent(in) :: f(dim) ! function for differentiation
      real(DP), intent(in) :: dx ! distance between neighbouring gridpoints
      integer :: index ! index of the point where we want the derivative
      complex(DP), dimension(dim) :: dfdx_1d
      real(DP), parameter :: c1 = 1.0d0 / 12.0d0
      real(DP), parameter :: c2 = 8.0d0 / 12.0d0

      do index = 1, dim
         if(index==1) dfdx_1d(index) = (f(index + 1) - f(index)) / dx
         if(index==2 .or. index==dim - 1) dfdx_1d(index) = 0.50d0 * (f(index + 1) - f(index - 1)) / dx
         if(index>2 .and. index<dim - 1) dfdx_1d(index) = (c1 * (f(index - 2) - f(index + 2)) + &
               c2 * (f(index + 1) - f(index - 1))) / dx
         if(index==dim) dfdx_1d(index) = -(f(index) - f(index - 1)) / dx
      end do

   end function

   !=== GI-TDPES ===!
   subroutine calc_gitdpes_1d(wf, nucdens, C, dC, A, gi_tdpes)
      real(DP), dimension(xngrid), intent(in) :: nucdens, A
      complex(DP), dimension(nstates, xngrid), intent(in) :: wf, C, dC
      real(DP), dimension(5, xngrid), intent(inout) :: gi_tdpes
      real(DP), dimension(xngrid) :: Vint
      integer :: istate, jstate

      gi_tdpes = 0.0d0

      do istate = 1, nstates
         do jstate = 1, nstates
            ! calculate <psi|H_el|psi> term
            gi_tdpes(1, :) = gi_tdpes(1, :) + REAL(CONJG(C(istate, :)) * H1(istate, jstate, :) * C(jstate, :))
            ! calculate <psi|V_int|psi> term
            if (field_coupling) then
               Vint = - dipole_coupling(istate, jstate, :) * elmag_field(time)
            else ! if the field is off, we don't have elmag_field, thus, we need to set up manually
               Vint = 0.0d0
            end if
            gi_tdpes(2, :) = gi_tdpes(2, :) + REAL(CONJG(C(istate, :)) * Vint(:) * C(jstate, :))
         end do
         ! calculate Sum_k [<grad_k psi| grad_k psi> / 2M_k] term
         gi_tdpes(3, :) = gi_tdpes(3, :) + 0.5D0 / mass_x * ABS(dC(istate, :))**2
      end do
      ! calculate Sum_k [A_k^2 / 2M_k] term
      gi_tdpes(4, :) = - 0.5D0 / mass_x * A(:)**2
      ! sum up all the terms in to the total GI-TDPES
      gi_tdpes(5, :) = gi_tdpes(1, :) + gi_tdpes(2, :) + gi_tdpes(3, :) + gi_tdpes(4, :)

   end subroutine

   !=== TIME DERIVATIVES ===!
   subroutine wf_time_der_1d(td_wf, sum_over_states, order)
      character(len = 14), intent(in) :: order
      complex(DP), dimension(nstates, xngrid), intent(inout) :: td_wf
      complex(DP), dimension(xngrid), intent(out) :: sum_over_states ! the summation of derivatives over all states in GD-TDPES

      integer :: state

      td_wf = 0.0d0
      sum_over_states = 0.0D0
      do state = 1, nstates
         ! there could be check if nuc_density_1d>efzero but I will do that once calculating gd_tdpes_1d
         ! checking for division by zero and calculating vecpot
         do i = 1, xngrid
            if (nucdens_1d(i) > ef_zero) then
               if (order == '5th-order_cent') then
                  td_wf(state, i) = ((-1.0D0 / 12.0D0) * (wfx_history(5, state, i) - wfx_history(1, state, i)) + &
                        (8.0D0 / 12.0D0) * (wfx_history(4, state, i) - wfx_history(2, state, i))) / dt
               elseif (order == '3th-order_forw') then
                  ! comes from https://web.media.mit.edu/~crtaylor/calculator.html
                  ! 	f_x = (-3*f[i+0]+4*f[i+1]-1*f[i+2])/(2*1.0*h**1)
                  td_wf(state, i) = (-3.0d0 * wfx_history(1, state, i) + 4.0d0 * wfx_history(2, state, i)&
                        - 1.0d0 * wfx_history(3, state, i)) / (2.0d0 * dt)
               elseif (order == '5th-order_back') then
                  ! comes from https://web.media.mit.edu/~crtaylor/calculator.html
                  ! f_x = (3*f[i-4]-16*f[i-3]+36*f[i-2]-48*f[i-1]+25*f[i+0])/(12*1.0*h**1)
                  td_wf(state, i) = (3.0d0 * wfx_history(1, state, i) - 16.0d0 * wfx_history(2, state, i) &
                        + 36.0d0 * wfx_history(3, state, i) - 48.0d0 * wfx_history(4, state, i) &
                        + 25.0d0 * wfx_history(5, state, i)) / (12.0d0 * dt)
               end if
               sum_over_states(i) = sum_over_states(i) + 0.5D0 * (conjg(td_wf(state, i)) * wfx_history(3, state, i) + &
                     conjg(wfx_history(3, state, i)) * td_wf(state, i)) / nucdens_1d(i)
            end if ! else it's zero which is achieve byt defining the variables as zero above
         end do
      end do

   end subroutine

   subroutine phase_time_der_1d(td_phase, phase_hist, order)
      real(DP), dimension(xngrid), intent(inout) :: td_phase
      real(DP), dimension(efhistory, xngrid), intent(in) :: phase_hist
      character(len = 14), intent(in) :: order

      td_phase = 0.0D0
      do i = 1, xngrid
         if (nucdens_1d(i) > ef_zero) then
            if (order == '5th-order_cent') then
               td_phase(i) = ((-1.0D0 / 12.0D0) * (phase_hist(5, i) - phase_hist(1, i)) + &
                     (8.0D0 / 12.0D0) * (phase_hist(4, i) - phase_hist(2, i))) / dt
            elseif (order == '3th-order_forw') then
               ! comes from https://web.media.mit.edu/~crtaylor/calculator.html
               ! 	f_x = (-3*f[i+0]+4*f[i+1]-1*f[i+2])/(2*1.0*h**1)
               td_phase(i) = (-3.0d0 * phase_hist(1, i) + 4.0d0 * phase_hist(2, i) - 1.0d0 * phase_hist(3, i)) / (2.0d0 * dt)
            elseif (order == '5th-order_back') then
               ! comes from https://web.media.mit.edu/~crtaylor/calculator.html
               ! f_x = (3*f[i-4]-16*f[i-3]+36*f[i-2]-48*f[i-1]+25*f[i+0])/(12*1.0*h**1)
               td_phase(i) = (3.0d0 * phase_hist(1, i) - 16.0d0 * phase_hist(2, i) + 36.0d0 * phase_hist(3, i)&
                     - 48.0d0 * phase_hist(4, i) + 25.0d0 * phase_hist(5, i)) / (12.0d0 * dt)
            end if
         end if
      end do
   end subroutine


   !=== SAVE HISTORY ===!
   subroutine ef_savehistory_1d()
      ! save wave function and phase history

      do j = 1, nstates ! looping through states
         do i = 1, efhistory - 1 ! looping through history
            wfx_history(efhistory + 1 - i, j, :) = wfx_history(efhistory - i, j, :) ! shifting previous wfx
         end do
         wfx_history(1, j, :) = wfx(j, :) ! saving current wfx
      end do

   end subroutine

   !=== INITIALIZATION ===!
   subroutine init_ef()
      write(*, *) "---------"
      write(*, *) "Exact factorization initialization"

      if(rank == 1) then
         allocate(wfx_history(efhistory, nstates, xngrid), nucdens_1d(xngrid), phase_1d(xngrid), grad_phase_1d(xngrid))
         allocate(vecpot_1d(xngrid), C_ef_1d(nstates, xngrid), gi_tdpes_1d(5, xngrid), gd_tdpes_1d(xngrid))
      end if

      ! nuclear density file initialization
      file_unit = 500
      file_name = 'nuclear_density.dat'
      open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
      if (iost /= 0) then
         write(*, *) "Error opening file: ", file_name
         stop 1
      end if
      if(rank == 1) then
         write(file_unit, *) "# x   rho_nuc"
      elseif(rank == 2) then
         write(file_unit, *) "# x   y   rho_nuc"
      elseif(rank == 3) then
         write(file_unit, *) "# x   y   z   rho_nuc"
      end if
      close(file_unit)
      open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)
      write(*, *)"Nuclear density file:           " // file_name

      ! nuclear phase and its gradient file initialization
      file_unit = 501
      file_name = 'nuclear_phase.dat'
      open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
      if (iost /= 0) then
         write(*, *) "Error opening file: ", file_name
         stop 1
      end if
      if(rank == 1) then
         write(file_unit, *) "# x   S   grad_S"
      elseif(rank == 2) then
         write(file_unit, *) "# x   y   S   grad_S"
      elseif(rank == 3) then
         write(file_unit, *) "# x   y   z   S   grad_S"
      end if
      close(file_unit)
      open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)
      write(*, *)"Nuclear phase file:             " // file_name

      ! vector potential file initialization
      file_unit = 502
      file_name = 'tdvp.dat'
      open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
      if (iost /= 0) then
         write(*, *) "Error opening file: ", file_name
         stop 1
      end if
      if(rank == 1) then
         write(file_unit, *) "# x   TDVP"
      elseif(rank == 2) then
         write(file_unit, *) "# x   y   TDVP"
      elseif(rank == 3) then
         write(file_unit, *) "# x   y   z   TDVP"
      end if
      close(file_unit)
      open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)
      write(*, *)"TD vector potential file:       " // file_name

      ! gauge-independent TDPES file initialization
      file_unit = 503
      file_name = 'gi-tdpes.dat'
      open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
      if (iost /= 0) then
         write(*, *) "Error opening file: ", file_name
         stop 1
      end if
      if(rank == 1) then
         write(file_unit, *) "#           x                     <psi|H_el|psi>             <psi|V_int|psi>   &
               Sum_k[<grad_k psi|grad_k psi>/2M_k]  Sum_k [A_k^2/2M_k]            total GI-TDPES"
      elseif(rank == 2) then
         write(file_unit, *) "# x   y   ;   <psi|H_el|psi>   ;   <psi|V_int|psi>   ;   Sum_k [<grad_k psi| grad_k psi> / 2M_k]   &
               ;   Sum_k [A_k^2 / 2M_k]   ;   total GI-TDPES"
      elseif(rank == 3) then
         write(file_unit, *) "# x   y   z   ;   <psi|H_el|psi>   ;   <psi|V_int|psi>   ;   &
               Sum_k [<grad_k psi| grad_k psi> / 2M_k]   ;   Sum_k [A_k^2 / 2M_k]   ;   total GI-TDPES"
      end if
      close(file_unit)
      open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)
      write(*, *)"Gauge-indep. TDPES file:        " // file_name

      ! gauge-dependent TDPES file initialization
      file_unit = 504
      file_name = 'gd-tdpes.dat'
      open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
      if (iost /= 0) then
         write(*, *) "Error opening file: ", file_name
         stop 1
      end if
      if(rank == 1) then
         write(file_unit, *) "# x   GD-TDPES"
      elseif(rank == 2) then
         write(file_unit, *) "# x   y   GD-TDPES"
      elseif(rank == 3) then
         write(file_unit, *) "# x   y   z   GD-TDPES"
      end if
      close(file_unit)
      open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)
      write(*, *)"Gauge-dep. TDPES file:          " // file_name

      ! electronic coefficients file initialization
      file_unit = 505
      file_name = 'el_coefficients_ef.dat'
      open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
      if (iost /= 0) then
         write(*, *) "Error opening file: ", file_name
         stop 1
      end if
      if(rank == 1) then
         write(file_unit, *) "# x   C_el [Re(C_el1), Im(C_el1), Re(C_el2), Im(C_el2), ...]"
      elseif(rank == 2) then
         write(file_unit, *) "# x   y   C_el [Re(C_el1), Im(C_el1), Re(C_el2), Im(C_el2), ...]"
      elseif(rank == 3) then
         write(file_unit, *) "# x   y   z   C_el [Re(C_el1), Im(C_el1), Re(C_el2), Im(C_el2), ...]"
      end if
      close(file_unit)
      open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)
      write(*, *)"El. coefficients file:          " // file_name

      select case(rank)
      case(1)
         call exact_factor_1d(0)
         call print_ef_1d()
      end select

   end subroutine

   !=== PRINTING ===!
   ! Print everything but GD-TDPES
   subroutine print_ef_1d()

      ! print nuclear density
      file_unit = 500
      write(file_unit, '(A,F10.3,A)') " # time ", time, " a.u."
      do i = 1, size(x)
         write(file_unit, *) x(i), nucdens_1d(i)
      end do

      ! print nuclear phase and its gradient
      file_unit = 501
      write(file_unit, '(A,F10.3,A)') " # time ", time, " a.u."
      do i = 1, size(x)
         write(file_unit, *) x(i), phase_1d(i), grad_phase_1d(i)
      end do

      ! print TDVP
      file_unit = 502
      write(file_unit, '(A,F10.3,A)') " # time ", time, " a.u."
      do i = 1, size(x)
         write(file_unit, *) x(i), vecpot_1d(i)
      end do

      ! print el. coefficients
      file_unit = 505
      write(file_unit, '(A,F10.3,A)') " # time ", time, " a.u."
      do i = 1, size(x)
         write(file_unit, '(F21.16, 30(E20.10, E20.10))') x(i), C_ef_1d(:, i)
      end do

      ! print GI-TDPES
      file_unit = 503
      write(file_unit, '(A,F10.3,A)') " # time ", time, " a.u."
      do i = 1, size(x)
         write(file_unit, '(F23.16,5E28.12)') x(i), gi_tdpes_1d(1, i), gi_tdpes_1d(2, i), gi_tdpes_1d(3, i), &
               gi_tdpes_1d(4, i), gi_tdpes_1d(5, i)
      end do

   end subroutine

   ! Print everything GD-TDPES only
   subroutine print_ef_gd_tdpes_1d(step)
      integer, intent(in) :: step
      real(DP) :: time

      if (step/=nstep) then ! time for GD-TDPES when we are calculating the time-derivative two steps later
         time = (step - efhistory / 2) * dt
      else ! for the last step, the time is the same
         time = step * dt
      end if

      ! print GD-TDPES
      file_unit = 504
      write(file_unit, '(A,F10.3,A)') " # time ", time, " a.u."
      do i = 1, size(x)
         write(file_unit, '(F23.16,E28.12)') x(i), gd_tdpes_1d(i)
      end do

   end subroutine

end module
