program qdyn
  use FFTW3
  use mod_init
  use mod_propag

! -------------------------------------------------------------------!
!                               / Qdyn \                             !
! Program for N-dimensional numerical quantum propagation on a grid. !
! Authors: Jiri Suchan, Petr Slavicek (2017)                         !
! -------------------------------------------------------------------!

   implicit none
   integer     :: n
!  real(DP)    :: 

!--Initialization

  call init()

!--Propagation mode--

select case(run)
  ! IMAGINARY TIME PROPAGATION
  case(0)
    write(*,*) 
    write(*,*) "RUN: 0 - IMAGINARY TIME PROPAGATION"
    write(*,*)
    write(*,*) "------ time -----"

if (rank .eq. 1) then
    do n=1, nstep 

      ! V(t/2)      
      do i=1, ngrid
        wfx(i) = wfx(i)*theta_v1(i)
      end do
      ! FFT -> K
      call dfftw_plan_dft_1d(plan_forward, ngrid, wfx, wfp, FFTW_FORWARD, FFTW_ESTIMATE )
      call dfftw_execute_dft(plan_forward, wfx, wfp)
      call dfftw_destroy_plan(plan_forward)

      ! p(t)
      do i=1, ngrid
        wfp(i) = wfp(i)*kin_p1(i)
      end do

      ! FFT -> x
      call dfftw_plan_dft_1d(plan_backward, ngrid, wfp, wfx, FFTW_BACKWARD, FFTW_ESTIMATE )
      call dfftw_execute_dft(plan_backward, wfp, wfx)
      call dfftw_destroy_plan(plan_backward)

      ! V(t/2)
      do i=1, ngrid
        wfx(i) = wfx(i)*theta_v1(i)
      end do

      !Normalize, print
      call normalize1d(wfx,ngrid,dx)
     
      if (modulo(dt*n,dtwrite) .eq. 0 ) then
        write(*,'(F5.1)') n*dt
        call printwf1d(wfx,x,v1)
      end if

    end do

elseif (rank .eq. 2) then

  do n=1, nstep

      ! V(t/2)      
      do i=1, ngrid
       do j=1, ngrid
        wf2x(i,j) = wf2x(i,j)*theta_v2(i,j)
       end do
      end do
      ! FFT -> K
      call dfftw_plan_dft_2d(plan_forward, ngrid, ngrid, wf2x, wf2p, FFTW_FORWARD, FFTW_ESTIMATE )
      call dfftw_execute_dft(plan_forward, wf2x, wf2p)
      call dfftw_destroy_plan(plan_forward)

      ! p(t)
      do i=1, ngrid
        do j=1, ngrid
          wf2p(i,j) = wf2p(i,j)*kin_p2(i,j)
        end do
      end do

      ! FFT -> x
      call dfftw_plan_dft_2d(plan_backward, ngrid, ngrid, wf2p, wf2x, FFTW_BACKWARD, FFTW_ESTIMATE )
      call dfftw_execute_dft(plan_backward, wf2p, wf2x)
      call dfftw_destroy_plan(plan_backward)

      ! V(t/2)
      do i=1, ngrid
        do j=1, ngrid
          wf2x(i,j) = wf2x(i,j)*theta_v2(i,j)
        end do
      end do

      !Normalize, print
      call normalize2d(wf2x,ngrid,dx)

      if (modulo(dt*n,dtwrite) .eq. 0 ) then
        write(*,'(F5.1)') n*dt
        call printwf2d(wf2x,x,y,v2)
      end if

    end do


endif

end select


end program qdyn

