module mod_utils

use mod_vars
use FFTW3

  implicit none

CONTAINS

!=== SCALAR PRODUCT ===!
! <bra|ket> (bra wavefunction is converted to conjugate in this function)
function braket_1d(bra, ket)

  complex(DP), intent(in)       :: bra(:), ket(:)
  real(DP)                      :: braket_1d
  integer                       :: i

  braket_1d=0.0d0

  do i=1,ngrid
    braket_1d = braket_1d + dx * conjg(bra(i))*ket(i)
  end do

end function braket_1d

function braket_2d(bra, ket )

  complex(DP), intent(in)       :: bra(:,:), ket(:,:)
  real(DP)                      :: braket_2d
  integer                       :: i, j

  braket_2d=0.0d0

  do i=1,ngrid
    do j=1, ngrid
      braket_2d = braket_2d + dx**2 * conjg(bra(i,j))*ket(i,j)
    end do
  end do

end function braket_2d

function braket_3d(bra, ket )

  complex(DP), intent(in)       :: bra(:,:,:), ket(:,:,:)
  real(DP)                      :: braket_3d
  integer                       :: i, j, k

  braket_3d=0.0d0
  
  do i=1,ngrid
    do j=1, ngrid
      do k=1, ngrid
        braket_3d = braket_3d + dx**3 * conjg(bra(i,j,k))*ket(i,j,k)
      end do
    end do
  end do

end function braket_3d

!=== PROJECT OUT STATE ===!
! |psi> = sum c_i |phi_i>
! |psi_new> = |psi> - <phi_i|psi>|phi_i> 
! Projecting out the i-th state wavefunction phi_i from the total wavefunction
subroutine project_out_1d(phi_i,wfx)

  complex(DP), intent(inout)    :: wfx(:)
  complex(DP), intent(in)       :: phi_i(:)
  real(DP)                      :: c_i

  c_i = braket_1d(phi_i, wfx)
  wfx = wfx - c_i * phi_i

end subroutine

!=== NORMALIZATION ===!
subroutine normalize_1d(wf) 

  complex(DP), intent(inout)    :: wf(:)
  real(DP)                      :: norm

  norm = braket_1d(wf, wf)

  !TODO: here should be norm check with some clever threshold
  ! currenlty very simple patch
  if (dabs(dsqrt(norm)-1.0d0) >= 1.0d-6 .and. run == 0) then
    write(*,*) "WARNING! Norm exceeded threshold", norm
  end if

  wf = wf / dsqrt(norm)

end subroutine normalize_1d

subroutine normalize_2d(wf)

  complex(DP), intent(inout)    :: wf(:,:)
  real(DP)                      :: norm

  norm = braket_2d(wf, wf)

  !TODO: here should be norm check with some clever threshold
  ! currenlty very simple patch
  if (dabs(dsqrt(norm)-1.0d0) >= 1.0d-6 .and. run == 0) then
    write(*,*) "WARNING! Norm exceeded threshold", norm
  end if

  wf = wf / dsqrt(norm)

end subroutine normalize_2d

subroutine normalize_3d(wf)

  complex(DP), intent(inout)    :: wf(:,:,:)
  real(DP)                      :: norm

  norm = braket_3d(wf, wf)

  !TODO: here should be norm check with some clever threshold
  ! currenlty very simple patch
  if (dabs(dsqrt(norm)-1.0d0) >= 1.0d-6 .and. run == 0) then
    write(*,*) "WARNING! Norm exceeded threshold", norm
  end if

  wf = wf / dsqrt(norm)

end subroutine normalize_3d

!=== PRINTING ===!
subroutine printwf_1d(wf,x,v1)

  complex(DP), intent(in)    :: wf(:)
  real(DP), intent(in)       :: x(:),v1(:)
  integer                    :: i

  do i=1, size(x)
    write(201,*) x(i), real(wf(i)), aimag(wf(i)), real(wf(i))**2+aimag(wf(i))**2, v1(i)
  end do
  write(201,*)
  write(201,*)

end subroutine


subroutine printwf_2d(wf,x,y,v2)

  complex(DP), intent(in)    :: wf(:,:)
  real(DP), intent(in)       :: x(:),y(:),v2(:,:)
  integer                    :: i,j

  do i=1, size(x)
    do j=1, size(y)
      write(202,*) x(i), y(j), real(wf(i,j)), aimag(wf(i,j)), real(wf(i,j))**2+aimag(wf(i,j))**2, v2(i,j)
    end do
  end do
  write(202,*)
  write(202,*)

end subroutine

subroutine printwf_3d(wf,x,y,z,v3)

  complex(DP), intent(in)    :: wf(:,:,:)
  real(DP), intent(in)       :: x(:),y(:),z(:),v3(:,:,:)
  integer                    :: i,j,k

  do i=1, size(x)
    do j=1, size(y)
      do k=1, size(z)
        write(203,*) x(i), y(j), z(k), real(wf(i,j,k)), aimag(wf(i,j,k)), real(wf(i,j,k))**2+aimag(wf(i,j,k))**2, v3(i,j,k)
      end do
     end do
  end do
  write(203,*)
  write(203,*)

end subroutine

subroutine print_chk()

  open(666,file='wf.chk', action='WRITE', iostat=iost)
  write(666,*) "#QDYN checkpoint file for reloading WF to program"
  write(666,*) "#Rank:",rank,"Pot:",pot,"Ngrid:",ngrid
  if(rank .eq. 1) write(666,*) wfx
  if(rank .eq. 2) write(666,*) wf2x
  if(rank .eq. 3) write(666,*) wf3x
  close(666)

end subroutine

subroutine printen()

    write(101,'(F8.1,3(F16.9))') time, energy
    
end subroutine

!=== ENERGIES ===!
!TODO: I should combine all of these. With rank, I can tell which one will be used.
subroutine update_energy_1d(wfx)
implicit none
  complex(DP), intent(in)       :: wfx(:)
  ! I probabaly do not need wfp(:)
  complex(DP), allocatable      :: wfp(:), wft(:)
  integer                       :: i

  allocate(wfp(ngrid))
  allocate(wft(ngrid))

  ! calculating <T>
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

  wft = wft / dsqrt(real(ngrid, kind=DP))

  !TODO: I should not need normalization here but I do it because I'm calculating during propagation
  ! to be deleted later
  energy(3) = braket_1d(wfx, wft)/braket_1d(wfx, wfx)

  ! calculating <V>
  energy(2) = braket_1d(wfx, v1*wfx)/braket_1d(wfx, wfx)

  ! calculating <E> = <V> + <T>
  energy(1) = energy(2)+ energy(3)

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
