module mod_utils

  implicit none
  public
  INTEGER, PARAMETER    :: DP = KIND(1.0d0)
  real(DP), parameter   :: pi = 3.14159265

CONTAINS

!=== NORMALIZATION ===!
subroutine normalize1d(wf,ngrid,dx) 

  complex(DP), intent(inout)    :: wf(:)
  integer, intent(in)           :: ngrid
  real(DP), intent(in)          :: dx
  real(DP)                      :: norm
  integer                       :: i

  norm=0.0d0

  do i=1,ngrid
    norm = norm + dx * (real(wf(i))**2+aimag(wf(i))**2)
  end do

  wf = wf / sqrt(norm)

end subroutine normalize1d

subroutine normalize2d(wf,ngrid,dx)

  complex(DP), intent(inout)    :: wf(:,:)
  integer, intent(in)           :: ngrid
  real(DP), intent(in)          :: dx
  real(DP)                      :: norm
  integer                       :: i, j

  norm=0.0d0

  do i=1,ngrid
    do j=1, ngrid
      norm = norm + dx**2 * (real(wf(i,j))**2+aimag(wf(i,j))**2)
    end do
  end do

  wf = wf / sqrt(norm)

end subroutine normalize2d

subroutine normalize3d(wf,ngrid,dx)

  complex(DP), intent(inout)    :: wf(:,:,:)
  integer, intent(in)           :: ngrid
  real(DP), intent(in)          :: dx
  real(DP)                      :: norm
  integer                       :: i, j, k

  norm=0.0d0

  do i=1,ngrid
    do j=1, ngrid
      do k=1, ngrid
        norm = norm + dx**3 * (real(wf(i,j,k))**2+aimag(wf(i,j,k))**2)
      end do
    end do
  end do

  wf = wf / sqrt(norm)

end subroutine normalize3d

!=== PRINTING ===!
subroutine printwf1d(wf,x,v1)

  complex(DP), intent(in)    :: wf(:)
  real(DP), intent(in)       :: x(:),v1(:)
  integer                    :: i

  do i=1, size(x)
    write(201,*) x(i), real(wf(i)), aimag(wf(i)), real(wf(i))**2+aimag(wf(i))**2, v1(i)
  end do
  write(201,*)
  write(201,*)

end subroutine


subroutine printwf2d(wf,x,y,v2)

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

subroutine printwf3d(wf,x,y,z,v3)

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

end module
