module mod_energies

use FFTW3
use mod_vars
use mod_utils,  only: braket_1d, braket_2d, braket_3d, normalize_1d


implicit none

CONTAINS

!=== ENERGIES ===!
subroutine update_energy_1d(wfx, energy)
implicit none
  complex(DP), intent(in)       :: wfx(:)
  real(DP), intent(inout)       :: energy
  complex(DP), allocatable      :: h_wfx(:), wfp(:), wft(:)
  integer                       :: i

  allocate(h_wfx(ngrid))
  allocate(wfp(ngrid))
  allocate(wft(ngrid))

  ! calculating T(psi)
  ! FFT -> K
  call dfftw_plan_dft_1d(plan_forward, ngrid, wfx, wfp, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wfx, wfp)
  call dfftw_destroy_plan(plan_forward)

  wfp = wfp / dsqrt(real(ngrid, kind=DP))

  ! p(t)
  do i=1, ngrid
    wfp(i) = wfp(i)*px(i)**2/(2*mass)
  end do

  ! FFT -> x
  call dfftw_plan_dft_1d(plan_backward, ngrid, wfp, wft, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wfp, wft)
  call dfftw_destroy_plan(plan_backward)
  ! heck
  wft = wft / dsqrt(real(ngrid, kind=DP))

  ! Hamiltonian applied to the wavefunction
  ! H(psi) = T(psi) + V*psi
  do i=1, ngrid
    h_wfx(i) = wft(i) + wfx(i)*v1(i)
  end do

  energy = braket_1d(wfx, h_wfx)

end subroutine

subroutine update_energy_2d(wf2x, energy)
implicit none
  complex(DP), intent(in)       :: wf2x(:, :)
  real(DP), intent(inout)       :: energy
  complex(DP), allocatable      :: h_wf2x(:, :), wf2p(:, :), wf2t(:, :)
  integer                       :: i, j

  allocate(h_wf2x(ngrid, ngrid))
  allocate(wf2p(ngrid, ngrid))
  allocate(wf2t(ngrid, ngrid))


  ! calculating T(psi)
  ! FFT -> K
  call dfftw_plan_dft_2d(plan_forward, ngrid, ngrid, wf2x, wf2p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf2x, wf2p)
  call dfftw_destroy_plan(plan_forward)

  wf2p = wf2p / ngrid

  ! p(t)
  do i=1, ngrid
    do j=1, ngrid
      wf2p(i,j) = wf2p(i,j)*(px(i)**2 + py(j)**2)/(2*mass)
    end do
  end do

  ! FFT -> x
  call dfftw_plan_dft_2d(plan_backward, ngrid, ngrid, wf2p, wf2t, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wf2p, wf2t)
  call dfftw_destroy_plan(plan_backward)

  wf2t = wf2t / ngrid

  ! Hamiltonian applied to the wavefunction
  ! H(psi) = T(psi) + V*psi
  do i=1, ngrid
   do j=1, ngrid
    h_wf2x(i,j) = wf2t(i,j) + wf2x(i,j)*v2(i,j)
   end do
  end do

  energy = braket_2d(wf2x, h_wf2x)

end subroutine

subroutine update_energy_3d(wf3x, energy)
implicit none
  complex(DP), intent(in)       :: wf3x(:, :, :)
  real(DP), intent(inout)       :: energy
  complex(DP), allocatable      :: h_wf3x(:, :, :), wf3p(:, :, :), wf3t(:, :, :)
  integer                       :: i, j, k

  allocate(h_wf3x(ngrid, ngrid, ngrid))
  allocate(wf3p(ngrid, ngrid, ngrid))
  allocate(wf3t(ngrid, ngrid, ngrid))


  ! calculating T(psi)
  ! FFT -> K
  call dfftw_plan_dft_3d(plan_forward, ngrid, ngrid, ngrid, wf3x, wf3p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf3x, wf3p)
  call dfftw_destroy_plan(plan_forward)

  wf3p = wf3p / dsqrt(real(ngrid, kind=DP)**3)

  ! p(t)
  do i=1, ngrid
    do j=1, ngrid
      do k=1, ngrid
        wf3p(i,j,k) = wf3p(i,j,k)*(px(i)**2 + py(j)**2 + pz(k)**2)/(2*mass)
      end do
    end do
  end do

  ! FFT -> x
  call dfftw_plan_dft_3d(plan_backward, ngrid, ngrid, ngrid, wf3p, wf3t, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wf3p, wf3t)
  call dfftw_destroy_plan(plan_backward)

  wf3t = wf3t / dsqrt(real(ngrid, kind=DP)**3)

  ! Hamiltonian applied to the wavefunction
  ! H(psi) = T(psi) + V*psi
  do i=1, ngrid
   do j=1, ngrid
     do k=1, ngrid
       h_wf3x(i,j,k) = wf3t(i,j,k) + wf3x(i,j,k)*v3(i,j,k)
     end do
    end do
  end do

  energy = braket_3d(wf3x, h_wf3x)

end subroutine

end module
