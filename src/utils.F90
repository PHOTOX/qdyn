module mod_utils

   use mod_vars
   use FFTW3
   use fparser, ONLY : initf, parsef, evalf, EvalErrType, EvalErrMsg

   implicit none
   integer, private :: file_unit

CONTAINS

   !=== SCALAR PRODUCT ===!
   ! <bra|ket> (bra wavefunction is converted to conjugate in this function)
   function braket_1d(bra, ket)

      complex(DP), intent(in) :: bra(:), ket(:)
      complex(DP) :: braket_1d
      integer :: i

      braket_1d = 0.0d0

      do i = 1, xngrid
         braket_1d = braket_1d + dx * conjg(bra(i)) * ket(i)
      end do

   end function braket_1d

   function braket_2d(bra, ket)

      complex(DP), intent(in) :: bra(:, :), ket(:, :)
      real(DP) :: braket_2d
      integer :: i, j

      braket_2d = 0.0d0

      do i = 1, xngrid
         do j = 1, yngrid
            braket_2d = braket_2d + dx * dy * conjg(bra(i, j)) * ket(i, j)
         end do
      end do

   end function braket_2d

   function braket_3d(bra, ket)

      complex(DP), intent(in) :: bra(:, :, :), ket(:, :, :)
      real(DP) :: braket_3d
      integer :: i, j, k

      braket_3d = 0.0d0

      do i = 1, xngrid
         do j = 1, yngrid
            do k = 1, zngrid
               braket_3d = braket_3d + dx * dy * dz * conjg(bra(i, j, k)) * ket(i, j, k)
            end do
         end do
      end do

   end function braket_3d

   !=== PROJECT OUT STATE ===!
   ! |psi> = sum c_i |phi_i>
   ! |psi_new> = |psi> - <phi_i|psi>|phi_i>
   ! Projecting out the i-th state wavefunction phi_i from the total wavefunction
   subroutine project_out_1d(phi_i, wfx)

      complex(DP), intent(inout) :: wfx(:)
      complex(DP), intent(in) :: phi_i(:)
      complex(DP) :: rot(xngrid)
      real(DP) :: c_i

      c_i = braket_1d(phi_i, wfx)
      wfx = wfx - c_i * phi_i

      ! WARNING: numerical instability causes optimization to lower states rotated by 90 degrees in the imaginary plane.
      !It is purely numerical and can be removed by projecting out 90 degrees rotated wf.
      if (project_rot) then
         rot = dcmplx(0.0d0, 1.0d0) * phi_i
         c_i = braket_1d(rot, wfx)
         wfx = wfx - c_i * rot
      end if

   end subroutine

   subroutine project_out_2d(phi_i, wf2x)

      complex(DP), intent(inout) :: wf2x(:, :)
      complex(DP), intent(in) :: phi_i(:, :)
      complex(DP) :: rot(xngrid, yngrid)
      real(DP) :: c_i

      c_i = braket_2d(phi_i, wf2x)
      wf2x = wf2x - c_i * phi_i

      ! WARNING: numerical instability causes optimization to lower states rotated by 90 degrees in the imaginary plane.
      !It is purely numerical and can be removed by projecting out 90 degrees rotated wf.
      if (project_rot) then
         rot = dcmplx(0.0d0, 1.0d0) * phi_i
         c_i = braket_2d(rot, wf2x)
         wf2x = wf2x - c_i * rot
      end if

   end subroutine

   subroutine project_out_3d(phi_i, wf3x)

      complex(DP), intent(inout) :: wf3x(:, :, :)
      complex(DP), intent(in) :: phi_i(:, :, :)
      complex(DP) :: rot(xngrid, yngrid, zngrid)
      real(DP) :: c_i

      c_i = braket_3d(phi_i, wf3x)
      wf3x = wf3x - c_i * phi_i

      ! WARNING: numerical instability causes optimization to lower states rotated by 90 degrees in the imaginary plane.
      !It is purely numerical and can be removed by projecting out 90 degrees rotated wf.
      if (project_rot) then
         rot = dcmplx(0.0d0, 1.0d0) * phi_i
         c_i = braket_3d(rot, wf3x)
         wf3x = wf3x - c_i * rot
      end if

   end subroutine

   !=== NORMALIZATION ===!
   subroutine update_norm()
      implicit none
      integer :: istate

      norm = 0.0d0
      select case(rank)
      case(1)
         do istate = 1, nstates
            norm = norm + braket_1d(wfx(istate, :), wfx(istate, :))
         end do
      case(2)
         norm = braket_2d(wf2x(1, :, :), wf2x(1, :, :))
      case(3)
         norm = braket_3d(wf3x(1, :, :, :), wf3x(1, :, :, :))
      end select

      !TODO: here should be norm check with some clever threshold
      ! currenlty very simple patch
      if (dabs(dsqrt(norm) - 1.0d0) >= norm_thresh .and. run == 0) then
         write(*, '(a,F10.8,a)') "WARNING! Norm (", norm, ") exceeded threshold"
         write(*, *) "Renormalization!"
         select case(rank)
         case(1)
            wfx(:, :) = wfx(:, :) / dsqrt(norm)
            norm = 0.0d0
            do istate = 1, nstates
               norm = norm + braket_1d(wfx(istate, :), wfx(istate, :))
            end do
         case(2)
            call normalize_2d(wf2x(1, :, :))
            norm = braket_2d(wf2x(1, :, :), wf2x(1, :, :))
         case(3)
            call normalize_3d(wf3x(1, :, :, :))
            norm = braket_3d(wf3x(1, :, :, :), wf3x(1, :, :, :))
         end select

      end if

   end subroutine update_norm

   subroutine normalize_1d(wf)
      ! currently normalizing only first state
      complex(DP), intent(inout) :: wf(:)
      real(DP) :: norm

      norm = braket_1d(wf, wf)

      wf = wf / dsqrt(norm)

   end subroutine normalize_1d

   subroutine normalize_2d(wf)

      complex(DP), intent(inout) :: wf(:, :)
      real(DP) :: norm

      norm = braket_2d(wf, wf)

      wf = wf / dsqrt(norm)

   end subroutine normalize_2d

   subroutine normalize_3d(wf)

      complex(DP), intent(inout) :: wf(:, :, :)
      real(DP) :: norm

      norm = braket_3d(wf, wf)

      wf = wf / dsqrt(norm)

   end subroutine normalize_3d

   !=== POPULATIONS ===!
   subroutine update_pop()
      implicit none
      integer :: istate

      diab_pop = 0.0d0
      ad_pop = 0.0d0

      select case(rank)
      case(1)
         do istate = 1, nstates
            diab_pop(istate) = braket_1d(wfx(istate, :), wfx(istate, :))
            ad_pop(istate) = braket_1d(wfx_ad(istate, :), wfx_ad(istate, :))
         end do
      case(2)
         do istate = 1, nstates
            diab_pop(istate) = braket_2d(wf2x(istate, :, :), wf2x(istate, :, :))
         end do
      case(3)
         do istate = 1, nstates
            diab_pop(istate) = braket_3d(wf3x(istate, :, :, :), wf3x(istate, :, :, :))
         end do
      end select

   end subroutine update_pop

   !=== Diagonalize H and build expH ===!
   subroutine build_expH1()
      ! diagonalization of H and building of expH1
      integer :: istate, jstate
      real(DP) :: H_el(nstates, nstates, xngrid) ! H1 + V_int, this Hamiltonian is diagonalized
      real(DP) :: V_int(xngrid)
      ! auxiliary variables to calculate diagonal expH for 2 states
      real(DP) :: ampl, D
      complex(DP) :: prefactor
      ! auxiliary variables for numerical diagonalization
      real(DP), dimension(:, :, :), allocatable :: U1, invU1
      real(DP), dimension(:, :), allocatable :: eigvals
      integer :: ioerr, lwork, dim_work, liwork, dim_iwork, ldwork, dim_dwork
      integer, allocatable :: iwork(:), ipivot(:)
      real(kind = 8), allocatable :: work(:)

      ! setting H_el to H1
      H_el = H1

      ! if the field is on, we add the field couplings to the H_el matrix
      ! creating H_el=H1-mu*E(t)
      if (field_coupling) then
         do istate = 1, nstates
            do jstate = 1, nstates
               V_int(:) = - dipole_coupling(istate, jstate, :) * elmag_field(time)
               H_el(istate, jstate, :) = H_el(istate, jstate, :) + V_int(:)
            end do
         end do
      end if

      ! building expH1 for different dimensions
      if (nstates.eq.1) then ! no diagonalization is necessary

         do i = 1, xngrid
            expH1(1, 1, i) = dcmplx(dcos(-H_el(1, 1, i) * dt / 2.0d0), dsin(-H_el(1, 1, i) * dt / 2.0d0)) !exp(-i H(x) tau/(2 h_bar))
         end do

      else if (nstates.eq.2) then ! diagonalized H and directly constructed exp, see Tannor chapter 11.7 (eq 11.204)

         do i = 1, xngrid
            ampl = - (H_el(1, 1, i) + H_el(2, 2, i)) * dt / 4.0d0
            prefactor = dcmplx(dcos(ampl), dsin(ampl))
            !D = sqrt(4 * |V21|^2 + (V1-V2)^2)
            D = dsqrt(4.0d0 * H_el(2, 1, i)**2.0d0 + (H_el(1, 1, i) - H_el(2, 2, i))**2.0d0)

            expH1(1, 1, i) = prefactor * dcmplx(dcos(D * dt / 4.0d0), dsin(D * dt / 4.0d0) / D * (H_el(2, 2, i) - H_el(1, 1, i)))
            expH1(2, 2, i) = prefactor * dcmplx(dcos(D * dt / 4.0d0), dsin(D * dt / 4.0d0) / D * (H_el(1, 1, i) - H_el(2, 2, i)))
            expH1(1, 2, i) = prefactor * dcmplx(0.0d0, dsin(D * dt / 4.0d0) / D * (-2.0d0 * H_el(1, 2, i)))
            expH1(2, 1, i) = prefactor * dcmplx(0.0d0, dsin(D * dt / 4.0d0) / D * (-2.0d0 * H_el(2, 1, i)))

         end do

      elseif (nstates.gt.2) then ! use numerical diagonalization of H_el for more than 2 states

         ! allocate U1, invU1 and eigvals
         allocate(U1(nstates, nstates, xngrid), invU1(nstates, nstates, xngrid), eigvals(nstates, xngrid))

         !--- U1 CALCULATION ---!
         ! U1 needs to contain the matrix that will be diagonalized, it will contain output eigenvectors
         U1 = H_el

         ! allocating dimensions for dsyevd
         allocate(work(1), iwork(1))
         lwork = -1
         liwork = -1
         call dsyevd('V', 'U', nstates, U1(:, :, 1), nstates, eigvals(:, 1), &
               work, lwork, iwork, liwork, ioerr)
         if(ioerr/=0) write(*, *) 'ERROR: preparing matrix diagonalization in expH1 construction failed'
         dim_work = work(1)
         dim_iwork = iwork(1)
         deallocate(work, iwork)
         allocate(work(dim_work), iwork(dim_iwork))

         ! performing real diagonalization
         U1 = H_el
         do i = 1, xngrid
            lwork = 1 + 6 * nstates + 2 * nstates**2
            liwork = 3 + 5 * nstates
            call dsyevd('V', 'U', nstates, U1(:, :, i), nstates, eigvals(:, i), &
                  work, lwork, iwork, liwork, ioerr)
            if(ioerr/=0) write(*, *) 'ERROR: matrix diagonalization in expH1 construction failed'
            if(i.gt.1) call check_overlap(U1(:, :, i), U1(:, :, i - 1))
         end do
         deallocate(work, iwork)

         !--- INVERSION OF U1 ---!
         ! first getting the right dimensions for arrays
         dim_dwork = 1
         allocate(work(dim_dwork), ipivot(nstates))
         ldwork = -1
         invU1(:, :, 1) = U1(:, :, 1)
         call dgetri(nstates, invU1, nstates, ipivot, work, ldwork, ioerr)
         if(ioerr/=0) print*, 'ERROR: preparation for inverting U in expH1 construction failed'
         dim_dwork = work(1)
         ldwork = nstates
         deallocate(work)
         allocate(work(dim_dwork))

         ! calculating the inversion
         do i = 1, xngrid
            invU1(:, :, i) = U1(:, :, i)
            ! first factorizing the U matrix so that it can be inverted later
            call dgetrf(nstates, nstates, invU1(:, :, i), nstates, ipivot, ioerr)
            if(ioerr/=0) print*, 'ERROR: factorizing error when inverting U in expH1 construction failed'
            ! now inverting the matrix using the factorized form
            call dgetri(nstates, invU1(:, :, i), nstates, ipivot, work, ldwork, ioerr)
            if(ioerr/=0) print*, 'ERROR: inverting transformation matrix U in expH1 construction failed'
         end do
         deallocate(work, ipivot)

         !--- expH1 CALCULATION ---!
         ! expH1 = U * exp(-iE dt/2) * U^(-1)
         expH1 = 0.0d0 ! initializing expH1 such that I can use it for eigenvalues
         ! put eigenvalues at the diagonal, offdiagonal elements are zero
         do i = 1, xngrid
            ! First prepare exp(-iE dt/2) in the expH1 matrix
            do istate = 1, nstates
               expH1(istate, istate, i) = dcmplx(dcos(-eigvals(istate, i) * dt / 2.0d0), dsin(-eigvals(istate, i) * dt / 2.0d0))
            end do
            expH1(:, :, i) = matmul(U1(:, :, i), matmul(expH1(:, :, i), invU1(:, :, i)))
         end do

      end if
   end subroutine build_expH1

   !=== ADIABATIC TRANSFORMATION ===!
   subroutine adiab_trans_matrix()
      integer :: ioerr, lwork, dim_work, liwork, dim_iwork, ldwork, dim_dwork
      integer, allocatable :: iwork(:), ipivot(:)
      real(kind = 8), allocatable :: work(:)

      write(*, *) "Adiabatic energies and transformation matrix U calculated."
      ! U1 needs to contain the matrix that will be diagonalized, it will contain output eigenvectors
      U1 = H1

      ! allocating dimensions for dsyevd
      allocate(work(1), iwork(1))
      lwork = -1
      liwork = -1
      call dsyevd('V', 'U', nstates, U1(:, :, 1), nstates, H1_ad(:, 1), &
            work, lwork, iwork, liwork, ioerr)
      if(ioerr/=0) write(*, *) 'ERROR: preparing matrix diagonalization'
      dim_work = work(1)
      dim_iwork = iwork(1)
      deallocate(work, iwork)
      allocate(work(dim_work), iwork(dim_iwork))

      ! performing real diagonalization
      U1 = H1
      do i = 1, xngrid
         lwork = 1 + 6 * nstates + 2 * nstates**2
         liwork = 3 + 5 * nstates
         call dsyevd('V', 'U', nstates, U1(:, :, i), nstates, H1_ad(:, i), &
               work, lwork, iwork, liwork, ioerr)
         if(ioerr/=0) write(*, *) 'ERROR: matrix diagonalization'
         if(i.gt.1) call check_overlap(U1(:, :, i), U1(:, :, i - 1))
      end do
      deallocate(work, iwork)

      ! calculating inverse of U
      ! U = from adiabatic to diabatic
      ! invU = from diabatic to adiabatic
      ! first getting the right dimensions for arrays
      dim_dwork = 1
      allocate(work(dim_dwork), ipivot(nstates))
      ldwork = -1
      invU1(:, :, 1) = U1(:, :, 1)
      call dgetri(nstates, invU1, nstates, ipivot, work, ldwork, ioerr)
      if(ioerr/=0) print*, 'ERROR: preparation for inverting U'
      dim_dwork = work(1)
      ldwork = nstates
      deallocate(work)
      allocate(work(dim_dwork))

      ! calculating the inversion
      do i = 1, xngrid
         invU1(:, :, i) = U1(:, :, i)
         ! first factorizing the U matrix so that it can be inverted later
         call dgetrf(nstates, nstates, invU1(:, :, i), nstates, ipivot, ioerr)
         if(ioerr/=0) print*, 'ERROR: factorizing error when inverting U'
         ! now inverting the matrix using the factorized form
         call dgetri(nstates, invU1(:, :, i), nstates, ipivot, work, ldwork, ioerr)
         if(ioerr/=0) print*, 'ERROR: inverting transformation matrix U'
      end do
      deallocate(work, ipivot)

   end subroutine adiab_trans_matrix

   subroutine check_overlap(eigenv1, eigenv2)

      real(kind = 8), intent(inout) :: eigenv1(nstates, nstates)
      real(kind = 8), intent(in) :: eigenv2(nstates, nstates)
      integer :: istate

      do istate = 1, nstates
         if (eigenv1(1, istate) * eigenv2(1, istate) + eigenv1(2, istate) * eigenv2(2, istate)<0.0d0) then
            eigenv1(:, istate) = -eigenv1(:, istate)
         end if
      end do

   end subroutine check_overlap

   subroutine wf_adiab_trans()
      integer :: istate, jstate

      select case(rank)
      case(1)
         wfx_ad = 0.0d0
         do i = 1, xngrid
            wfx_ad(:, i) = matmul(invU1(:, :, i), wfx(:, i))
         end do
      end select

   end subroutine wf_adiab_trans

   !=== PRINTING ===!
   subroutine printwf_1d(state)

      integer, intent(in) :: state
      integer :: i
      real(DP) :: v

      file_unit = 200 + state

      write(file_unit, '(A,F10.3,A)') "#time ", time, " a.u."

      do i = 1, size(x)
         select case(run)
         case(0)
            v = H1(state, state, i)
         case(1)
            v = v1(i)  ! V is same for all states
         end select

         write(file_unit, *) x(i), real(wfx(state, i)), aimag(wfx(state, i)), real(conjg(wfx(state, i)) * wfx(state, i)), v
      end do

      write(file_unit, *)

   end subroutine

   subroutine printwf_ad_1d(state)

      integer, intent(in) :: state
      integer :: i
      real(DP) :: v

      file_unit = 400 + state

      write(file_unit, '(A,F10.3,A)') "#time ", time, " a.u."

      do i = 1, size(x)

         v = H1_ad(state, i)

         write(file_unit, *) x(i), real(wfx_ad(state, i)), aimag(wfx_ad(state, i)), &
               real(conjg(wfx_ad(state, i)) * wfx_ad(state, i)), v
      end do

      write(file_unit, *)

   end subroutine

   subroutine printwf_2d(state)

      integer, intent(in) :: state
      integer :: i, j

      file_unit = 200 + state

      write(file_unit, '(A,F10.3,A)') "#time ", time, " a.u."

      do i = 1, size(x)
         do j = 1, size(y)
            write(file_unit, *) x(i), y(j), real(wf2x(state, i, j)), aimag(wf2x(state, i, j)), &
                  real(conjg(wf2x(state, i, j)) * wf2x(state, i, j)), v2(i, j)
         end do
      end do

      write(file_unit, *)

   end subroutine

   subroutine printwf_3d(state)

      integer, intent(in) :: state
      integer :: i, j, k

      file_unit = 200 + state

      write(file_unit, '(A,F10.3,A)') "#time ", time, " a.u."

      do i = 1, size(x)
         do j = 1, size(y)
            do k = 1, size(z)
               write(file_unit, *) x(i), y(j), z(k), real(wf3x(state, i, j, k)), aimag(wf3x(state, i, j, k)), &
                     real(conjg(wf3x(state, i, j, k)) * wf3x(state, i, j, k)), v3(i, j, k)
            end do
         end do
      end do

      write(file_unit, *)

   end subroutine

   subroutine print_chk(state)
      integer, intent(in) :: state

      open(666, file = 'wf.chk', action = 'WRITE', iostat = iost)
      write(666, '(A,I3)') "#QDYN wave function file after IT propagation for state ", state
      write(666, *) "#Rank:", rank, "Ngrid:", xngrid, yngrid, zngrid
      if(rank .eq. 1) write(666, *) wfx(state, :)
      if(rank .eq. 2) write(666, *) wf2x(state, :, :)
      if(rank .eq. 3) write(666, *) wf3x(state, :, :, :)
      close(666)

   end subroutine

   subroutine printen()

      write(101, '(F8.1,5(F16.9))') time, energy, energy_diff, norm

   end subroutine

   subroutine printen_state(state)

      integer, intent(in) :: state
      file_unit = 300 + state
      write(file_unit, '(F8.1,5(F16.9))') time, energy, energy_diff

   end subroutine

   subroutine print_field()

      write(102, '(F10.3,F14.9)') time, elmag_field(time)

   end subroutine

   subroutine print_pop()

      call update_pop()

      write(103, '(F10.3,10F14.9)') time, diab_pop, norm
      write(104, '(F10.3,10F14.9)') time, ad_pop, norm

   end subroutine

   !=== ENERGIES ===!
   !TODO: components of energy here
   function kinetic_energy_1d(wfx)
      implicit none
      complex(DP), intent(in) :: wfx(:)
      complex(DP), allocatable :: wft(:)
      real(DP) :: kinetic_energy_1d
      integer :: i

      allocate(wft(xngrid))

      ! calculating <T>
      ! iFFT -> k
      call dfftw_plan_dft_1d(plan_backward, xngrid, wfx, wfp, FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(plan_backward, wfx, wfp)
      call dfftw_destroy_plan(plan_backward)

      wfp = wfp / dsqrt(real(xngrid, kind = DP))

      ! p(t)
      do i = 1, xngrid
         wft(i) = wfp(i) * px(i)**2 / (2 * mass_x)
      end do

      ! calculating <T> in the momentum representation, back FFT to coordinate represetation was skipped
      kinetic_energy_1d = braket_1d(wfp, wft)

   end function

   subroutine update_energy_1d(wfx)
      implicit none
      complex(DP), intent(in) :: wfx(:)
      real(DP) :: old_energy

      ! calculating <T>
      energy(3) = kinetic_energy_1d(wfx)

      ! calculating <V>
      energy(2) = braket_1d(wfx, v1 * wfx) / braket_1d(wfx, wfx)

      ! saving old energy
      old_energy = energy(1)

      ! calculating <E> = <V> + <T>
      energy(1) = energy(2) + energy(3)

      ! energy difference from the last step
      energy_diff = energy(1) - old_energy

   end subroutine

   subroutine update_energy_2d(wf2x)
      implicit none
      complex(DP), intent(in) :: wf2x(:, :)
      real(DP) :: old_energy
      complex(DP), allocatable :: wf2t(:, :)
      integer :: i, j

      allocate(wf2t(xngrid, yngrid))

      ! calculating T(psi)
      ! iFFT -> k
      call dfftw_plan_dft_2d(plan_backward, xngrid, yngrid, wf2x, wf2p, FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(plan_backward, wf2x, wf2p)
      call dfftw_destroy_plan(plan_backward)

      wf2p = wf2p / dsqrt(real(xngrid * yngrid, kind = DP))

      ! p(t)
      do i = 1, xngrid
         do j = 1, yngrid
            wf2t(i, j) = wf2p(i, j) * (px(i)**2 / (2 * mass_x) + py(j)**2 / (2 * mass_y))
         end do
      end do

      ! calculating <T> in the momentum representation, back FFT to coordinate represetation was skipped
      energy(3) = braket_2d(wf2p, wf2t) / braket_2d(wf2p, wf2p)
      ! Note that the braket_1d uses dx for the integration instead of dp, which is not correct
      ! but in the fraction, the dx/dx will cancel out so it works. For more complex stuff, a
      ! problem can appear here!

      ! calculating <V>
      energy(2) = braket_2d(wf2x, v2 * wf2x) / braket_2d(wf2x, wf2x)

      ! saving old energy
      old_energy = energy(1)

      ! calculating <E> = <V> + <T>
      energy(1) = energy(2) + energy(3)

      ! energy difference from the last step
      energy_diff = energy(1) - old_energy

   end subroutine

   subroutine update_energy_3d(wf3x)
      implicit none
      complex(DP), intent(in) :: wf3x(:, :, :)
      real(DP) :: old_energy
      complex(DP), allocatable :: wf3t(:, :, :)
      integer :: i, j, k

      allocate(wf3t(xngrid, yngrid, zngrid))


      ! calculating T(psi)
      ! iFFT -> k
      call dfftw_plan_dft_3d(plan_backward, xngrid, yngrid, zngrid, wf3x, wf3p, FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(plan_backward, wf3x, wf3p)
      call dfftw_destroy_plan(plan_backward)

      wf3p = wf3p / dsqrt(real(xngrid * yngrid * zngrid, kind = DP))

      ! p(t)
      do i = 1, xngrid
         do j = 1, yngrid
            do k = 1, zngrid
               wf3t(i, j, k) = wf3p(i, j, k) * (px(i)**2 / (2 * mass_x) + py(j)**2 / (2 * mass_y) + pz(k)**2 / (2 * mass_z))
            end do
         end do
      end do

      ! calculating <T> in the momentum representation, back FFT to coordinate represetation was skipped
      energy(3) = braket_3d(wf3p, wf3t) / braket_3d(wf3p, wf3p)
      ! Note that the braket_1d uses dx for the integration instead of dp, which is not correct
      ! but in the fraction, the dx/dx will cancel out so it works. For more complex stuff, a
      ! problem can appear here!

      ! calculating <V>
      energy(2) = braket_3d(wf3x, v3 * wf3x) / braket_3d(wf3x, wf3x)

      ! saving old energy
      old_energy = energy(1)

      ! calculating <E> = <V> + <T>
      energy(1) = energy(2) + energy(3)

      ! energy difference from the last step
      energy_diff = energy(1) - old_energy

   end subroutine

   subroutine update_total_energy_1d()
      implicit none
      real(DP) :: old_energy, local_norm, thresh = 1.0d-6
      complex(DP) :: temp_en
      integer :: i, istate, jstate

      ! saving old energy
      old_energy = energy(1)

      energy = 0.0d0

      ! calculating <T>
      do istate = 1, nstates
         local_norm = braket_1d(wfx(istate, :), wfx(istate, :))
         if (local_norm.gt.thresh) then
            energy(3) = energy(3) + kinetic_energy_1d(wfx(istate, :))
         end if
      end do

      ! calculating <V>
      do istate = 1, nstates
         do jstate = istate, nstates
            ! if i=j then <i|H_el|j> is real number
            if (istate.eq.jstate) energy(2) = energy(2) + &
                  braket_1d(wfx(istate, :), H1(istate, jstate, :) * wfx(jstate, :))
            ! if i/=j then <i|H_el|j> is a complex number and <i|H_el|j> + <j|H_el|i> = 2*Re{<i|H_el|j>}
            if (istate.ne.jstate) energy(2) = energy(2) + &
                  2 * real(braket_1d(wfx(istate, :), H1(istate, jstate, :) * wfx(jstate, :)))
         end do
      end do

      ! calculating <E> = <V> + <T>
      energy(1) = energy(2) + energy(3)

      energy(:) = energy(:) / norm

      ! energy difference from the last step
      energy_diff = energy(1) - old_energy

   end subroutine

   function elmag_field(t)
      real(DP) :: t, elmag_field
      real(DP) :: point(1) !This is just a simple trick I use to create a 1D array out of t which is necessary for evalf

      point(1) = t
      if ((t.lt.field_on).or.(t.gt.field_off)) then
         elmag_field = 0.0d0
      else
         elmag_field = evalf(2, point)                   !Evaluating field at time t
         if (EvalErrType > 0) then
            WRITE(*, *)'*** Error evaluating potential: ', EvalErrMsg ()
         end if
      end if

   end function

end module
