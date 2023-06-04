module mod_utils

use mod_vars

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

subroutine printen(time,energy)

    real(DP), intent(in)     :: time, energy

    write(101,'(F8.1,F15.8)') time, energy
    
end subroutine

end module
