module mod_propag
use FFTW3
use mod_utils 
use mod_vars

implicit none

CONTAINS

!--REAL TIME PROPAGATION
subroutine propag_rt_1d()

  ! Propagating with H matrix: half step
  call propag_H_rt_1d()

  ! Propagating istate with T matrix: full step
  call propag_T_1d()

  ! Propagating with H matrix: half step
  call propag_H_rt_1d()

end subroutine propag_rt_1d

subroutine propag_H_rt_1d()
! This function propagates wave function with the H_el hamiltonian
implicit none
  integer                                  :: i, istate, jstate
  complex(DP), dimension(nstates,xngrid)   :: wfx_tmp ! temporary wf for propagation

  wfx_tmp=0.0d0

  if (field_coupling) call build_expH1()
  
  do istate=1,nstates ! states to be propagated
    do jstate=1,nstates ! this loop goes through all the states contributing to propagation of istate
      do i=1, xngrid
        wfx_tmp(istate,i) = wfx_tmp(istate,i)+expH1(istate,jstate,i)*wfx(jstate,i)
      end do
    end do
  end do

  ! Saving propagated wave function
  wfx = wfx_tmp

end subroutine propag_H_rt_1d

subroutine propag_rt_2d(wf2x)

implicit none
  complex(DP), intent(inout)    :: wf2x(:,:)
  integer                       :: i,j

  ! V(t/2)
  do i=1, xngrid
   do j=1, yngrid
    wf2x(i,j) = wf2x(i,j)*expV2(i,j)
   end do
  end do

  ! FFT -> K
  call dfftw_plan_dft_2d(plan_forward, xngrid, yngrid, wf2x, wf2p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf2x, wf2p)
  call dfftw_destroy_plan(plan_forward)

  wf2p = wf2p / dsqrt(real(xngrid*yngrid, kind=DP))

  ! p(t)
  do i=1, xngrid
    do j=1, yngrid
      wf2p(i,j) = wf2p(i,j)*expT2(i,j)
    end do
  end do

  ! FFT -> x
  call dfftw_plan_dft_2d(plan_backward, xngrid, yngrid, wf2p, wf2x, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wf2p, wf2x)
  call dfftw_destroy_plan(plan_backward)

  wf2x = wf2x / dsqrt(real(xngrid*yngrid, kind=DP))

  ! V(t/2)
  do i=1, xngrid
    do j=1, yngrid
      wf2x(i,j) = wf2x(i,j)*expV2(i,j)
    end do
  end do

end subroutine propag_rt_2d

subroutine propag_rt_3d(wf3x)

implicit none
  complex(DP), intent(inout)    :: wf3x(:,:,:)
  integer                       :: i,j,k

  ! V(t/2)
  do i=1, xngrid
   do j=1, yngrid
     do k=1, zngrid
       wf3x(i,j,k) = wf3x(i,j,k)*expV3(i,j,k)
     end do
    end do
  end do

  ! FFT -> K
  call dfftw_plan_dft_3d(plan_forward, xngrid, yngrid, zngrid, wf3x, wf3p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf3x, wf3p)
  call dfftw_destroy_plan(plan_forward)

  wf3p = wf3p / dsqrt(real(xngrid*yngrid*zngrid, kind=DP))

  ! p(t)
  do i=1, xngrid
    do j=1, yngrid
      do k=1, zngrid
        wf3p(i,j,k) = wf3p(i,j,k)*expT3(i,j,k)
      end do
    end do
  end do

  ! FFT -> x
  call dfftw_plan_dft_3d(plan_backward, xngrid, yngrid, zngrid, wf3p, wf3x, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wf3p, wf3x)
  call dfftw_destroy_plan(plan_backward)

  wf3x = wf3x / dsqrt(real(xngrid*yngrid*zngrid, kind=DP))

  ! V(t/2)
  do i=1, xngrid
    do j=1, yngrid
      do k=1, zngrid
        wf3x(i,j,k) = wf3x(i,j,k)*expV3(i,j,k)
      end do
    end do
  end do

end subroutine propag_rt_3d


!--IMAGINARY TIME PROPAGATION
subroutine propag_it_1d(wfx)

implicit none
  complex(DP), intent(inout)    :: wfx(:)
  integer                       :: i

  !TODO: half steps are demanding.. full steps would be better
  ! V(t/2)
  do i=1, xngrid
    wfx(i) = wfx(i)*expV1(i)
  end do

  ! FFT -> K
  call dfftw_plan_dft_1d(plan_forward, xngrid, wfx, wfp, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wfx, wfp)
  call dfftw_destroy_plan(plan_forward)

  wfp = wfp / dsqrt(real(xngrid, kind=DP))

  ! p(t)
  do i=1, xngrid
    wfp(i) = wfp(i)*expT1(i)
  end do

  ! FFT -> x
  call dfftw_plan_dft_1d(plan_backward, xngrid, wfp, wfx, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wfp, wfx)
  call dfftw_destroy_plan(plan_backward)

  wfx = wfx / dsqrt(real(xngrid, kind=DP))

  ! V(t/2)
  do i=1, xngrid
    wfx(i) = wfx(i)*expV1(i)
  end do

end subroutine propag_it_1d

subroutine propag_it_2d(wf2x)

implicit none
  complex(DP), intent(inout)    :: wf2x(:,:)
  integer                       :: i,j

  ! V(t/2)
  do i=1, xngrid
   do j=1, yngrid
    wf2x(i,j) = wf2x(i,j)*expV2(i,j)
   end do
  end do

  ! FFT -> K
  call dfftw_plan_dft_2d(plan_forward, xngrid, yngrid, wf2x, wf2p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf2x, wf2p)
  call dfftw_destroy_plan(plan_forward)

  wf2p = wf2p / dsqrt(real(xngrid*yngrid, kind=DP))

  ! p(t)
  do i=1, xngrid
    do j=1, yngrid
      wf2p(i,j) = wf2p(i,j)*expT2(i,j)
    end do
  end do

  ! FFT -> x
  call dfftw_plan_dft_2d(plan_backward, xngrid, yngrid, wf2p, wf2x, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wf2p, wf2x)
  call dfftw_destroy_plan(plan_backward)

  wf2x = wf2x / dsqrt(real(xngrid*yngrid, kind=DP))

  ! V(t/2)
  do i=1, xngrid
    do j=1, yngrid
      wf2x(i,j) = wf2x(i,j)*expV2(i,j)
    end do
  end do

end subroutine propag_it_2d

subroutine propag_it_3d(wf3x)

implicit none
  complex(DP), intent(inout)    :: wf3x(:,:,:)
  integer                       :: i,j,k

  ! V(t/2)
  do i=1, xngrid
   do j=1, yngrid
     do k=1, zngrid
       wf3x(i,j,k) = wf3x(i,j,k)*expV3(i,j,k)
     end do
    end do
  end do

  ! FFT -> K
  call dfftw_plan_dft_3d(plan_forward, xngrid, yngrid, zngrid, wf3x, wf3p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf3x, wf3p)
  call dfftw_destroy_plan(plan_forward)

  wf3p = wf3p / dsqrt(real(xngrid*yngrid*zngrid, kind=DP))

  ! p(t)
  do i=1, xngrid
    do j=1, yngrid
      do k=1, zngrid
        wf3p(i,j,k) = wf3p(i,j,k)*expT3(i,j,k)
      end do
    end do
  end do

  ! FFT -> x
  call dfftw_plan_dft_3d(plan_backward, xngrid, yngrid, zngrid, wf3p, wf3x, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wf3p, wf3x)
  call dfftw_destroy_plan(plan_backward)

  wf3x = wf3x / dsqrt(real(xngrid*yngrid*zngrid, kind=DP))

  ! V(t/2)
  do i=1, xngrid
    do j=1, yngrid
      do k=1, zngrid
        wf3x(i,j,k) = wf3x(i,j,k)*expV3(i,j,k)
      end do
    end do
  end do

end subroutine propag_it_3d

!--GENERAL TIME PROPAGATION FUNCTIONS
subroutine propag_T_1d()
! This function propagates the wave function with kinetic operator T
implicit none
  integer    :: i, istate

  do istate=1,nstates ! states to be propagated
    ! FFT -> K
    call dfftw_plan_dft_1d(plan_forward, xngrid, wfx(istate,:), wfp, FFTW_FORWARD, FFTW_ESTIMATE )
    call dfftw_execute_dft(plan_forward, wfx(istate,:), wfp)
    call dfftw_destroy_plan(plan_forward)

    wfp = wfp / dsqrt(real(xngrid, kind=DP))

    ! p(t)
    do i=1, xngrid
      wfp(i) = wfp(i)*expT1(i)
    end do

    ! FFT -> x
    call dfftw_plan_dft_1d(plan_backward, xngrid, wfp, wfx(istate,:), FFTW_BACKWARD, FFTW_ESTIMATE )
    call dfftw_execute_dft(plan_backward, wfp, wfx(istate,:))
    call dfftw_destroy_plan(plan_backward)

    wfx(istate,:) = wfx(istate,:) / dsqrt(real(xngrid, kind=DP))
  end do

end subroutine propag_T_1d

end module 
