module mod_utils

use mod_vars
use FFTW3
use fparser,    ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg

  implicit none
  integer,private     :: file_unit

CONTAINS

!=== SCALAR PRODUCT ===!
! <bra|ket> (bra wavefunction is converted to conjugate in this function)
function braket_1d(bra, ket)

  complex(DP), intent(in)       :: bra(:), ket(:)
  real(DP)                      :: braket_1d
  integer                       :: i

  braket_1d=0.0d0

  do i=1, xngrid
    braket_1d = braket_1d + dx * conjg(bra(i))*ket(i)
  end do

end function braket_1d

function braket_2d(bra, ket )

  complex(DP), intent(in)       :: bra(:,:), ket(:,:)
  real(DP)                      :: braket_2d
  integer                       :: i, j

  braket_2d=0.0d0

  do i=1, xngrid
    do j=1, yngrid
      braket_2d = braket_2d + dx * dy * conjg(bra(i,j))*ket(i,j)
    end do
  end do

end function braket_2d

function braket_3d(bra, ket )

  complex(DP), intent(in)       :: bra(:,:,:), ket(:,:,:)
  real(DP)                      :: braket_3d
  integer                       :: i, j, k

  braket_3d=0.0d0
  
  do i=1, xngrid
    do j=1, yngrid
      do k=1, zngrid
        braket_3d = braket_3d + dx * dy * dz * conjg(bra(i,j,k))*ket(i,j,k)
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
  complex(DP)                   :: rot(xngrid)
  real(DP)                      :: c_i

  c_i = braket_1d(phi_i, wfx)
  wfx = wfx - c_i * phi_i

  ! WARNING: numerical instability causes optimization to lower states rotated by 90 degrees in the imaginary plane. 
  !It is purely numerical and can be removed by projecting out 90 degrees rotated wf.
  if (project_rot) then
    rot = cmplx(0.0d0, 1.0d0)*phi_i
    c_i = braket_1d(rot, wfx)
    wfx = wfx - c_i * rot
  end if

end subroutine

subroutine project_out_2d(phi_i,wf2x)

  complex(DP), intent(inout)    :: wf2x(:,:)
  complex(DP), intent(in)       :: phi_i(:,:)
  complex(DP)                   :: rot(xngrid,yngrid)
  real(DP)                      :: c_i

  c_i = braket_2d(phi_i, wf2x)
  wf2x = wf2x - c_i * phi_i

  ! WARNING: numerical instability causes optimization to lower states rotated by 90 degrees in the imaginary plane. 
  !It is purely numerical and can be removed by projecting out 90 degrees rotated wf.
  if (project_rot) then
    rot = cmplx(0.0d0, 1.0d0)*phi_i
    c_i = braket_2d(rot, wf2x)
    wf2x = wf2x - c_i * rot
  end if

end subroutine

subroutine project_out_3d(phi_i,wf3x)

  complex(DP), intent(inout)    :: wf3x(:,:,:)
  complex(DP), intent(in)       :: phi_i(:,:,:)
  complex(DP)                   :: rot(xngrid,yngrid,zngrid)
  real(DP)                      :: c_i

  c_i = braket_3d(phi_i, wf3x)
  wf3x = wf3x - c_i * phi_i

  ! WARNING: numerical instability causes optimization to lower states rotated by 90 degrees in the imaginary plane. 
  !It is purely numerical and can be removed by projecting out 90 degrees rotated wf.
  if (project_rot) then
    rot = cmplx(0.0d0, 1.0d0)*phi_i
    c_i = braket_3d(rot, wf3x)
    wf3x = wf3x - c_i * rot
  end if

end subroutine

!=== NORMALIZATION ===!
subroutine update_norm() 

  select case(rank)
    case(1)
      norm = braket_1d(wfx(1,:), wfx(1,:))
    case(2)
      norm = braket_2d(wf2x(1,:,:), wf2x(1,:,:))
    case(3)
      norm = braket_3d(wf3x(1,:,:,:), wf3x(1,:,:,:))
  end select

  !TODO: here should be norm check with some clever threshold
  ! currenlty very simple patch
  if (dabs(dsqrt(norm)-1.0d0) >= 1.0d-6 .and. run == 0) then
    write(*,'(a,F10.8,a)') "WARNING! Norm (",norm,") exceeded threshold"
    select case(rank)
      case(1)
        write(*,*) "Renormalization!"
        call normalize_1d(wfx(1,:))
        norm = braket_1d(wfx(1,:), wfx(1,:))
      case(2)
        write(*,*) "Renormalization!"
        call normalize_2d(wf2x(1,:,:))
        norm = braket_2d(wf2x(1,:,:), wf2x(1,:,:))
      case(3)
        write(*,*) "Renormalization!"
        call normalize_3d(wf3x(1,:,:,:))
        norm = braket_3d(wf3x(1,:,:,:), wf3x(1,:,:,:))
    end select

  end if

end subroutine update_norm

subroutine normalize_1d(wf) 
! currently normalizing only first state
  complex(DP), intent(inout)    :: wf(:)
  real(DP)                      :: norm

  norm = braket_1d(wf, wf)

  wf = wf / dsqrt(norm)

end subroutine normalize_1d

subroutine normalize_2d(wf)

  complex(DP), intent(inout)    :: wf(:,:)
  real(DP)                      :: norm

  norm = braket_2d(wf, wf)

  wf = wf / dsqrt(norm)

end subroutine normalize_2d

subroutine normalize_3d(wf)

  complex(DP), intent(inout)    :: wf(:,:,:)
  real(DP)                      :: norm

  norm = braket_3d(wf, wf)

  wf = wf / dsqrt(norm)

end subroutine normalize_3d

!=== PRINTING ===!
subroutine printwf_1d(state)

  integer, intent(in)        :: state
  integer                    :: i
  real(DP)                   :: v

  file_unit=200+state

  write(file_unit,'(A,F10.3,A)') "#time ", time, " a.u."

  do i=1, size(x)
    select case(run)
    case(0)
      v = H1(state,state,i)
    case(1)
      v = v1(i)  ! V is same for all states
    end select

    write(file_unit,*) x(i), real(wfx(state,i)), aimag(wfx(state,i)), real(conjg(wfx(state,i))*wfx(state,i)), v
  end do

  write(file_unit,*)

end subroutine

subroutine printwf_2d(state)

  integer, intent(in)        :: state
  integer                    :: i,j

  file_unit=200+state

  write(file_unit,'(A,F10.3,A)') "#time ", time, " a.u."

  do i=1, size(x)
    do j=1, size(y)
      write(file_unit,*) x(i), y(j), real(wf2x(state,i,j)), aimag(wf2x(state,i,j)), real(conjg(wf2x(state,i,j))*wf2x(state,i,j)),&
      v2(i,j)
    end do
  end do

  write(file_unit,*)

end subroutine

subroutine printwf_3d(state)

  integer, intent(in)        :: state
  integer                    :: i,j,k

  file_unit=200+state

  write(file_unit,'(A,F10.3,A)') "#time ", time, " a.u."

  do i=1, size(x)
    do j=1, size(y)
      do k=1, size(z)
        write(file_unit,*) x(i), y(j), z(k), real(wf3x(state,i,j,k)), aimag(wf3x(state,i,j,k)), &
          real(conjg(wf3x(state,i,j,k))*wf3x(state,i,j,k)), v3(i,j,k)
      end do
     end do
  end do

  write(file_unit,*)

end subroutine

subroutine print_chk()

  open(666,file='wf.chk', action='WRITE', iostat=iost)
  write(666,*) "#QDYN checkpoint file for reloading WF to program"
  write(666,*) "#Rank:",rank,"Pot:",pot,"Ngrid:",xngrid
  if(rank .eq. 1) write(666,*) wfx
  if(rank .eq. 2) write(666,*) wf2x
  if(rank .eq. 3) write(666,*) wf3x
  close(666)

end subroutine

subroutine printen()

    write(101,'(F8.1,5(F16.9))') time, energy, energy_diff, norm
    
end subroutine

subroutine printen_state(state)

  integer, intent(in)           :: state
  file_unit = 300+state
    write(file_unit,'(F8.1,5(F16.9))') time, energy, energy_diff

end subroutine

subroutine print_field()

    write(102,'(F10.3,F14.9)') time, elmag_field(time)
    
end subroutine

!=== ENERGIES ===!
!TODO: I should combine all of these. With rank, I can tell which one will be used.
subroutine update_energy_1d(wfx)
implicit none
  complex(DP), intent(in)       :: wfx(:)
  complex(DP), allocatable      :: wft(:)
  real(DP)                      :: old_energy
  integer                       :: i

  allocate(wft(xngrid))

  ! calculating <T>
  ! FFT -> K
  call dfftw_plan_dft_1d(plan_forward, xngrid, wfx, wfp, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wfx, wfp)
  call dfftw_destroy_plan(plan_forward)

  wfp = wfp / dsqrt(real(xngrid, kind=DP))

  ! p(t)
  do i=1, xngrid
    wft(i) = wfp(i)*px(i)**2/(2*mass_x)
  end do

  ! calculating <T> in the momentum representation, back FFT to coordinate represetation was skipped
  energy(3) = braket_1d(wfp, wft)/braket_1d(wfp, wfp)
  ! Note that the braket_1d uses dx for the integration instead of dp, which is not correct
  ! but in the fraction, the dx/dx will cancel out so it works. For more complex stuff, a 
  ! problem can appear here!

  ! calculating <V>
  energy(2) = braket_1d(wfx, v1*wfx)/braket_1d(wfx, wfx)

  ! saving old energy
  old_energy = energy(1)

  ! calculating <E> = <V> + <T>
  energy(1) = energy(2)+ energy(3)
  
  ! energy difference from the last step
  energy_diff = energy(1) - old_energy

end subroutine

subroutine update_energy_2d(wf2x)
implicit none
  complex(DP), intent(in)       :: wf2x(:, :)
  real(DP)                      :: old_energy
  complex(DP), allocatable      :: wf2t(:, :)
  integer                       :: i, j

  allocate(wf2t(xngrid, yngrid))

  ! calculating T(psi)
  ! FFT -> K
  call dfftw_plan_dft_2d(plan_forward, xngrid, yngrid, wf2x, wf2p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf2x, wf2p)
  call dfftw_destroy_plan(plan_forward)

  wf2p = wf2p / dsqrt(real(xngrid*yngrid, kind=DP))

  ! p(t)
  do i=1, xngrid
    do j=1, yngrid
      wf2t(i,j) = wf2p(i,j)*(px(i)**2/(2*mass_x) + py(j)**2/(2*mass_y))
    end do
  end do

  ! calculating <T> in the momentum representation, back FFT to coordinate represetation was skipped
  energy(3) = braket_2d(wf2p, wf2t)/braket_2d(wf2p, wf2p)
  ! Note that the braket_1d uses dx for the integration instead of dp, which is not correct
  ! but in the fraction, the dx/dx will cancel out so it works. For more complex stuff, a 
  ! problem can appear here!

  ! calculating <V>
  energy(2) = braket_2d(wf2x, v2*wf2x)/braket_2d(wf2x, wf2x)

  ! saving old energy
  old_energy = energy(1)

  ! calculating <E> = <V> + <T>
  energy(1) = energy(2)+ energy(3)
  
  ! energy difference from the last step
  energy_diff = energy(1) - old_energy

end subroutine

subroutine update_energy_3d(wf3x)
implicit none
  complex(DP), intent(in)       :: wf3x(:, :, :)
  real(DP)                      :: old_energy
  complex(DP), allocatable      :: wf3t(:, :, :)
  integer                       :: i, j, k

  allocate(wf3t(xngrid, yngrid, zngrid))


  ! calculating T(psi)
  ! FFT -> K
  call dfftw_plan_dft_3d(plan_forward, xngrid, yngrid, zngrid, wf3x, wf3p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf3x, wf3p)
  call dfftw_destroy_plan(plan_forward)

  wf3p = wf3p / dsqrt(real(xngrid*yngrid*zngrid, kind=DP))

  ! p(t)
  do i=1, xngrid
    do j=1, yngrid
      do k=1, zngrid
        wf3t(i,j,k) = wf3p(i,j,k)*(px(i)**2/(2*mass_x) + py(j)**2/(2*mass_y) + pz(k)**2/(2*mass_z))
      end do
    end do
  end do

  ! calculating <T> in the momentum representation, back FFT to coordinate represetation was skipped
  energy(3) = braket_3d(wf3p, wf3t)/braket_3d(wf3p, wf3p)
  ! Note that the braket_1d uses dx for the integration instead of dp, which is not correct
  ! but in the fraction, the dx/dx will cancel out so it works. For more complex stuff, a 
  ! problem can appear here!

  ! calculating <V>
  energy(2) = braket_3d(wf3x, v3*wf3x)/braket_3d(wf3x, wf3x)

  ! saving old energy
  old_energy = energy(1)

  ! calculating <E> = <V> + <T>
  energy(1) = energy(2)+ energy(3)
  
  ! energy difference from the last step
  energy_diff = energy(1) - old_energy

end subroutine

subroutine update_energy_1d_rt()
implicit none
  complex(DP), allocatable      :: wft(:)
  real(DP)                      :: old_energy
  integer                       :: i

  allocate(wft(xngrid))

  ! calculating <T>
  ! FFT -> K
  call dfftw_plan_dft_1d(plan_forward, xngrid, wfx(1,:), wfp, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wfx(1,:), wfp)
  call dfftw_destroy_plan(plan_forward)

  wfp = wfp / dsqrt(real(xngrid, kind=DP))

  ! p(t)
  do i=1, xngrid
    wfp(i) = wfp(i)*px(i)**2/(2*mass_x)
  end do

  !TODO: this back FFT is not necessary since I still have the wf. It should be removed.
  ! FFT -> x
  call dfftw_plan_dft_1d(plan_backward, xngrid, wfp, wft, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wfp, wft)
  call dfftw_destroy_plan(plan_backward)

  wft = wft / dsqrt(real(xngrid, kind=DP))

  ! calculating <T>
  energy(3) = braket_1d(wfx(1,:), wft)/braket_1d(wfx(1,:), wfx(1,:))

  ! calculating <V>
  energy(2) = braket_1d(wfx(1,:), H1(1,1,:)*wfx(1,:))/braket_1d(wfx(1,:), wfx(1,:))

  ! saving old energy
  old_energy = energy(1)

  ! calculating <E> = <V> + <T>
  energy(1) = energy(2)+ energy(3)
  
  ! energy difference from the last step
  energy_diff = energy(1) - old_energy

end subroutine

function elmag_field(t)
  real(DP)    :: t, elmag_field
  real(DP)    :: point(1) !This is just a simple trick I use to create a 1D array out of t which is necessary for evalf

  point(1) = t
  elmag_field = evalf(2, point)                   !Evaluating field at time t
  if (EvalErrType > 0) then
    WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
  end if

end function

end module
