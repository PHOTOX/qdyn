module mod_propag
use FFTW3
use mod_utils 
use mod_vars

implicit none

CONTAINS

!--REAL TIME PROPAGATION
subroutine propag_rt_1d()

implicit none
  integer                       :: i, istate, jstate
  !TODO: Can I modify the previous wavefunciton? I probably cannot do that for the offdiagonal elements

  do istate=1,nstates
    do jstate=1,nstates
      ! Diagonal elements in the V matrix: half step
      if (istate.eq.jstate) then
        !TODO: half steps are demanding.. full steps would be better
        ! V(t/2)
        do i=1, ngrid
        !TODO: expV1 must have index of state too
          wfx(istate,i) = wfx(istate,i)*expV1(i)
        end do

        !>jj
        ! Field propagation here
        if (field_coupling) then
          ! TODO: I need to add the dipole moment
          wfx(istate,:) = wfx(istate,:)*cmplx(cos(dipole_coupling(istate,jstate,:)*elmag_field(time)*dt/2.0d0),&
            sin(dipole_coupling(istate,jstate,:)*elmag_field(time)*dt/2.0d0))
        end if
        !<jj
      end if

      !if (istate.not.jstate) then
        !OFFdiagonal elements here
        !wfx(istate,:) = wfx(istate,:) + wfx(jstate,:)*coupling
      !end if

      ! Diagonal elements in the T matrix: full step
      if (istate.eq.jstate) then
        ! FFT -> K
        call dfftw_plan_dft_1d(plan_forward, ngrid, wfx(istate,:), wfp, FFTW_FORWARD, FFTW_ESTIMATE )
        call dfftw_execute_dft(plan_forward, wfx(istate,:), wfp)
        call dfftw_destroy_plan(plan_forward)

        ! AFTER FFT, divide by sqrt(ngrid). This was found empirically to work, but checking the package would be desirable.
        ! Using this, energy and norm are conserved. Applied throughout the whole code.
        wfp = wfp / dsqrt(real(ngrid, kind=DP))

        ! p(t)
        do i=1, ngrid
          wfp(i) = wfp(i)*expT1(i)
        end do

        ! FFT -> x
        call dfftw_plan_dft_1d(plan_backward, ngrid, wfp, wfx(istate,:), FFTW_BACKWARD, FFTW_ESTIMATE )
        call dfftw_execute_dft(plan_backward, wfp, wfx(istate,:))
        call dfftw_destroy_plan(plan_backward)

        wfx(istate,:) = wfx(istate,:) / dsqrt(real(ngrid, kind=DP))
      end if

      ! Diagonal elements in the V matrix: half step
      if (istate.eq.jstate) then
        ! V(t/2)
        do i=1, ngrid
          wfx(istate,i) = wfx(istate,i)*expV1(i)
        end do

        !>jj
        ! Field propagation here
        if (field_coupling) then
          wfx(istate,:) = wfx(istate,:)*cmplx(cos(dipole_coupling(istate,jstate,:)*elmag_field(time)*dt/2.0d0),&
            sin(dipole_coupling(istate,jstate,:)*elmag_field(time)*dt/2.0d0))
        end if
        !<jj
      end if
    end do
  end do

end subroutine propag_rt_1d

subroutine propag_rt_2d(wf2x)

implicit none
  complex(DP), intent(inout)    :: wf2x(:,:)
  integer                       :: i,j

  ! V(t/2)
  do i=1, ngrid
   do j=1, ngrid
    wf2x(i,j) = wf2x(i,j)*expV2(i,j)
   end do
  end do

  ! FFT -> K
  call dfftw_plan_dft_2d(plan_forward, ngrid, ngrid, wf2x, wf2p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf2x, wf2p)
  call dfftw_destroy_plan(plan_forward)

  wf2p = wf2p / real(ngrid, kind=DP)

  ! p(t)
  do i=1, ngrid
    do j=1, ngrid
      wf2p(i,j) = wf2p(i,j)*expT2(i,j)
    end do
  end do

  ! FFT -> x
  call dfftw_plan_dft_2d(plan_backward, ngrid, ngrid, wf2p, wf2x, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wf2p, wf2x)
  call dfftw_destroy_plan(plan_backward)

  wf2x = wf2x / real(ngrid, kind=DP)

  ! V(t/2)
  do i=1, ngrid
    do j=1, ngrid
      wf2x(i,j) = wf2x(i,j)*expV2(i,j)
    end do
  end do

end subroutine propag_rt_2d

subroutine propag_rt_3d(wf3x)

implicit none
  complex(DP), intent(inout)    :: wf3x(:,:,:)
  integer                       :: i,j,k

  ! V(t/2)
  do i=1, ngrid
   do j=1, ngrid
     do k=1, ngrid
       wf3x(i,j,k) = wf3x(i,j,k)*expV3(i,j,k)
     end do
    end do
  end do

  ! FFT -> K
  call dfftw_plan_dft_3d(plan_forward, ngrid, ngrid, ngrid, wf3x, wf3p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf3x, wf3p)
  call dfftw_destroy_plan(plan_forward)

  wf3p = wf3p / dsqrt(real(ngrid, kind=DP)**3)

  ! p(t)
  do i=1, ngrid
    do j=1, ngrid
      do k=1, ngrid
        wf3p(i,j,k) = wf3p(i,j,k)*expT3(i,j,k)
      end do
    end do
  end do

  ! FFT -> x
  call dfftw_plan_dft_3d(plan_backward, ngrid, ngrid, ngrid, wf3p, wf3x, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wf3p, wf3x)
  call dfftw_destroy_plan(plan_backward)

  wf3x = wf3x / dsqrt(real(ngrid, kind=DP)**3)

  ! V(t/2)
  do i=1, ngrid
    do j=1, ngrid
      do k=1, ngrid
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
  do i=1, ngrid
    wfx(i) = wfx(i)*expV1(i)
  end do

  ! FFT -> K
  call dfftw_plan_dft_1d(plan_forward, ngrid, wfx, wfp, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wfx, wfp)
  call dfftw_destroy_plan(plan_forward)

  wfp = wfp / dsqrt(real(ngrid, kind=DP))

  ! p(t)
  do i=1, ngrid
    wfp(i) = wfp(i)*expT1(i)
  end do

  ! FFT -> x
  call dfftw_plan_dft_1d(plan_backward, ngrid, wfp, wfx, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wfp, wfx)
  call dfftw_destroy_plan(plan_backward)

  wfx = wfx / dsqrt(real(ngrid, kind=DP))

  ! V(t/2)
  do i=1, ngrid
    wfx(i) = wfx(i)*expV1(i)
  end do

end subroutine propag_it_1d

subroutine propag_it_2d(wf2x)

implicit none
  complex(DP), intent(inout)    :: wf2x(:,:)
  integer                       :: i,j

  ! V(t/2)
  do i=1, ngrid
   do j=1, ngrid
    wf2x(i,j) = wf2x(i,j)*expV2(i,j)
   end do
  end do

  ! FFT -> K
  call dfftw_plan_dft_2d(plan_forward, ngrid, ngrid, wf2x, wf2p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf2x, wf2p)
  call dfftw_destroy_plan(plan_forward)

  wf2p = wf2p / real(ngrid, kind=DP)

  ! p(t)
  do i=1, ngrid
    do j=1, ngrid
      wf2p(i,j) = wf2p(i,j)*expT2(i,j)
    end do
  end do

  ! FFT -> x
  call dfftw_plan_dft_2d(plan_backward, ngrid, ngrid, wf2p, wf2x, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wf2p, wf2x)
  call dfftw_destroy_plan(plan_backward)

  wf2x = wf2x / real(ngrid, kind=DP)

  ! V(t/2)
  do i=1, ngrid
    do j=1, ngrid
      wf2x(i,j) = wf2x(i,j)*expV2(i,j)
    end do
  end do

end subroutine propag_it_2d

subroutine propag_it_3d(wf3x)

implicit none
  complex(DP), intent(inout)    :: wf3x(:,:,:)
  integer                       :: i,j,k

  ! V(t/2)
  do i=1, ngrid
   do j=1, ngrid
     do k=1, ngrid
       wf3x(i,j,k) = wf3x(i,j,k)*expV3(i,j,k)
     end do
    end do
  end do

  ! FFT -> K
  call dfftw_plan_dft_3d(plan_forward, ngrid, ngrid, ngrid, wf3x, wf3p, FFTW_FORWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_forward, wf3x, wf3p)
  call dfftw_destroy_plan(plan_forward)

  wf3p = wf3p / dsqrt(real(ngrid, kind=DP)**3)

  ! p(t)
  do i=1, ngrid
    do j=1, ngrid
      do k=1, ngrid
        wf3p(i,j,k) = wf3p(i,j,k)*expT3(i,j,k)
      end do
    end do
  end do

  ! FFT -> x
  call dfftw_plan_dft_3d(plan_backward, ngrid, ngrid, ngrid, wf3p, wf3x, FFTW_BACKWARD, FFTW_ESTIMATE )
  call dfftw_execute_dft(plan_backward, wf3p, wf3x)
  call dfftw_destroy_plan(plan_backward)

  wf3x = wf3x / dsqrt(real(ngrid, kind=DP)**3)

  ! V(t/2)
  do i=1, ngrid
    do j=1, ngrid
      do k=1, ngrid
        wf3x(i,j,k) = wf3x(i,j,k)*expV3(i,j,k)
      end do
    end do
  end do

end subroutine propag_it_3d


end module 
