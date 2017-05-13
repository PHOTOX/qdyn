module mod_init
  use FFTW3
  use fparser,    ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg
  use mod_utils

  implicit none
  public
  real(DP)              :: dt, xmin, xmax, xmean, stddev, dx, k_0, mass, dtwrite
  real(DP), dimension(:), allocatable    :: x,y,z,point
  complex(DP), dimension(:), allocatable :: wfx, wfp, theta_v1, kin_p1
  complex(DP), dimension(:,:), allocatable :: wf2x, wf2p, theta_v2, kin_p2
  complex(DP), dimension(:,:,:), allocatable :: wf3x, wf3p, theta_v3, kin_p3
  integer               :: run, nstep, ngrid, wf, rank, iost, i, j, k
  integer ( kind = 8 )  :: plan_forward, plan_backward

  real(DP), dimension(:), allocatable     :: v1, px, py, pz
  real(DP), dimension(:,:), allocatable   :: v2
  real(DP), dimension(:,:,:), allocatable :: v3
  character(len=50)             :: pot=''
  character(len=*),dimension(1),parameter :: var1 = (/'x'/)
  character(len=*),dimension(2),parameter :: var2 = (/'x','y'/)
  character(len=*),dimension(3),parameter :: var3 = (/'x','y','z'/)

  namelist /general/run,nstep,dt,dtwrite,ngrid,rank,xmin,xmax,mass,wf,pot

CONTAINS

subroutine init()
  implicit none

!-- Reading input file
  open(100,file='input.q', status='OLD', action='READ',delim='APOSTROPHE', iostat=iost)
  read(100, general, iostat=iost)
  if (iost.ne.0) then
    write(*,*)'ERROR: input.q file must be provided'
    write(*,*) iost
    stop 1
  end if
  close(100)

!-- Input check
  call check()

!-- GRID set-up

dx=(xmax-xmin)/(ngrid-1)

! allocating arrays
if(rank .eq. 1) allocate(x(ngrid),v1(ngrid),px(ngrid),point(1))
if(rank .eq. 2) allocate(x(ngrid),y(ngrid),v2(ngrid,ngrid),px(ngrid),py(ngrid),point(2))
if(rank .eq. 3) allocate(x(ngrid),y(ngrid),z(ngrid),v3(ngrid,ngrid,ngrid),px(ngrid),py(ngrid),pz(ngrid),point(3))

if(rank .eq. 1) allocate(wfx(ngrid), wfp(ngrid), theta_v1(ngrid), kin_p1(ngrid))
if(rank .eq. 2) allocate(wf2x(ngrid,ngrid), wf2p(ngrid,ngrid), theta_v2(ngrid,ngrid), kin_p2(ngrid,ngrid))
if(rank .eq. 3) allocate(wf3x(ngrid,ngrid,ngrid), wf3p(ngrid,ngrid,ngrid), theta_v3(ngrid,ngrid,ngrid), kin_p3(ngrid,ngrid,ngrid))

! setting up grid poitns
x(1)=xmin
if(rank .gt. 1) y(1)=xmin
if(rank .gt. 2) z(1)=xmin

do i=2, ngrid
  x(i) = x(i-1) + dx
  if(rank .gt. 1) y(i) = y(i-1) + dx
  if(rank .gt. 2) z(i) = z(i-1) + dx
end do

!-- Initialization of WafeFunction

if(wf .eq. 0) then                                                           ! generating gaussian wavepacket
  write(*,*) "Generating gaussian wave packet at the center of the grid."
  xmean=(xmax-xmin)/2.0d0+xmin
  stddev=(xmax-xmean)/20.0d0                                                  ! 5sigma - 96% of gaussian is on the grid 
  k_0 = sqrt(2*mass*0.5)                                                   ! sqrt(2*m*E)/h = k0
  do i=1, ngrid
    select case (rank)
      case (1)
        wfx(i) = cmplx(exp((-1.0d0*(x(i)-xmean)**2)/(2*stddev**2)) * cos(k_0 * x(i)), &
                       exp((-1.0d0*(x(i)-xmean)**2)/(2*stddev**2)) * sin(k_0 * x(i)) )  
      case (2)
        do j=1, ngrid
          wf2x(i,j) =  cmplx(exp(- (((x(i)-xmean)**2)/(2*stddev**2)) - (((y(j)-xmean)**2)/(2*stddev**2))) * cos(k_0*x(i)), &
                             exp(- (((x(i)-xmean)**2)/(2*stddev**2)) - (((y(j)-xmean)**2)/(2*stddev**2))) * sin(k_0*x(i)) )
        end do
      case (3)
        do j=1, ngrid
          do k=1, ngrid
            wf3x(i,j,k) =  cmplx(exp(- (((x(i)-xmean)**2)/(2*stddev**2)) - (((y(j)-xmean)**2)/(2*stddev**2)) &
                           - (((z(k)-xmean)**2)/(2*stddev**2))), &
                           exp(- (((x(i)-xmean)**2)/(2*stddev**2)) - (((y(j)-xmean)**2)/(2*stddev**2)) &
                           - (((z(k)-xmean)**2)/(2*stddev**2))))
          end do
        end do
    end select
  end do
elseif(wf .eq. 1) then
!Procedure for loading WF from file
end if

if(rank .eq. 1) call normalize1d(wfx,ngrid,dx) 
if(rank .eq. 2) call normalize2d(wf2x,ngrid,dx)
if(rank .eq. 3) call normalize3d(wf3x,ngrid,dx) 

!-- Initialization of FFT procedures
if(rank .eq. 1) then

call dfftw_plan_dft_1d(plan_forward, ngrid, wfx, wfp, FFTW_FORWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_forward, wfx, wfp)
call dfftw_destroy_plan(plan_forward)

call dfftw_plan_dft_1d(plan_backward, ngrid, wfp, wfx, FFTW_BACKWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_backward, wfp, wfx)
call dfftw_destroy_plan(plan_backward)

elseif(rank .eq. 2) then

call dfftw_plan_dft_2d(plan_forward, ngrid, ngrid, wf2x, wf2p, FFTW_FORWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_forward, wf2x, wf2p)
call dfftw_destroy_plan(plan_forward)

call dfftw_plan_dft_2d(plan_backward, ngrid, ngrid, wf2p, wf2x, FFTW_BACKWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_backward, wf2p, wf2x)
call dfftw_destroy_plan(plan_backward)

elseif(rank .eq. 3) then

call dfftw_plan_dft_3d(plan_forward, ngrid, ngrid, ngrid, wf3x, wf3p, FFTW_FORWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_forward, wf3x, wf3p)
call dfftw_destroy_plan(plan_forward)

call dfftw_plan_dft_3d(plan_backward, ngrid, ngrid, ngrid, wf3p, wf3x, FFTW_BACKWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_backward, wf3p, wf3x)
call dfftw_destroy_plan(plan_backward)

end if

!-- POTENTIAL energy init
  write(*,*)"Potential: ",pot
  call initf (1)                                                     !Initialization of parser

  if(rank .eq. 1) then
    call parsef (1, pot, var1)                                       !Bytcompiling function string  
    do i=1, ngrid
      point = x(i)
      v1(i) = evalf (1, point)                                       !Evaluating potential for grid
      if (EvalErrType > 0) then
        WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
        stop 1
      end if
      theta_v1(i) = cmplx(cos(-v1(i)*dt/2.0d0),sin(-v1(i)*dt/2.0d0))   !exp(-i V(x) tau/(2 h_bar))
    end do

!2D   
  elseif(rank .eq. 2) then
    call parsef (1, pot, var2)                                      
    do i=1, ngrid
      do j=1, ngrid
        point = (/x(i),y(j)/)
        v2(i,j) = evalf (1, point) 
        if (EvalErrType > 0) then
          WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
          stop 1
        end if
        theta_v2(i,j) = cmplx(cos(v2(i,j)*dt/2.0d0),sin(v2(i,j)*dt/2.0d0))   !exp(-i V(x) tau/(2 h_bar))
      end do
    end do
  
!3D
  elseif(rank .eq. 3) then
    call parsef (1, pot, var3)
    do i=1, ngrid
      do j=1, ngrid
        do k=1, ngrid
          point = (/x(i),y(j),z(k)/)
          v3(i,j,k) = evalf (1, point)
          if (EvalErrType > 0) then
            WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
            stop 1
          end if
          theta_v3(i,j,k) = cmplx(cos(v3(i,j,k)*dt/2.0d0),sin(v3(i,j,k)*dt/2.0d0))  
        end do
      end do
    end do 
  end if

!-- KINETIC energy init        exp[-iT/h_bar tau]

select case(rank)
  case(1)
    do i=1, ngrid
     if(i .lt. ngrid/2) then
       px(i) = 2*pi*i/(ngrid*dx)
     else
       px(i) = 2*pi*(i-ngrid)/(ngrid*dx)
     end if
     kin_p1(i) = cmplx(cos(-px(i)**2*dt/(2*mass)),sin(-px(i)**2*dt/(2*mass)))
    end do
  !2D
  case(2)

   do i=1, ngrid
     if(i .lt. ngrid/2) then
       px(i) = 2*pi*i/(ngrid*dx)
     else
       px(i) = 2*pi*(i-ngrid)/(ngrid*dx)
     end if
   end do

   do j=1, ngrid
     if(j .lt. ngrid/2) then
       py(j) = 2*pi*j/(ngrid*dx)
     else
       py(j) = 2*pi*(j-ngrid)/(ngrid*dx)
     end if
   end do

    do i=1, ngrid
      do j=1, ngrid
        kin_p2(i,j) = cmplx(cos(-(px(i)**2+py(j)**2)*dt/(2*mass)),sin(-(px(i)**2+py(j)**2)*dt/(2*mass)))
      end do
    end do

  !3D
  case(3)

   do i=1, ngrid
     if(i .lt. ngrid/2) then
       px(i) = 2*pi*i/(ngrid*dx)
     else
       px(i) = 2*pi*(i-ngrid)/(ngrid*dx)
     end if
   end do

   do j=1, ngrid
     if(j .lt. ngrid/2) then
       py(j) = 2*pi*j/(ngrid*dx)
     else
       py(j) = 2*pi*(j-ngrid)/(ngrid*dx)
     end if
   end do

   do k=1, ngrid
     if(k .lt. ngrid/2) then
       pz(k) = 2*pi*k/(ngrid*dx)
     else
       pz(k) = 2*pi*(k-ngrid)/(ngrid*dx)
     end if
   end do
 
    do i=1, ngrid
      do j=1, ngrid
        do k=1, ngrid
          kin_p3(i,j,k) = cmplx(cos(-(px(i)**2+py(j)**2+pz(k)**2)*dt/(2*mass)),&
                                sin(-(px(i)**2+py(j)**2+pz(k)**2)*dt/(2*mass)))
        end do
      end do
    end do
end select

!--Printing of WF

if(rank .eq. 1) then
  call normalize1d(wfx,ngrid,dx)

  open(201,file='wf1d.out', action='WRITE', iostat=iost)
  write(201,*) "#WF - QDYN output"
  write(201,*) "#x   REAL   IMAG   PROBABILITY-DENSITY POTENTIAL"
  close(201)
  open(201,file='wf1d.out', status='old', position='append', action='WRITE', iostat=iost)

  call printwf1d(wfx,x,v1)
  write(*,*)"Outputing WF to file wf1d.out"

elseif(rank .eq. 2) then
  call normalize2d(wf2x,ngrid,dx)

  open(202,file='wf2d.out', action='WRITE', iostat=iost)
  write(202,*) "#WF - QDYN output"
  write(202,*) "#x  y   REAL   IMAG   PROBABILITY-DENSITY  POTENTIAL"
  close(202)
  open(202,file='wf2d.out', status='old', position='append', action='WRITE', iostat=iost)

  call printwf2d(wf2x,x,y,v2)
  write(*,*)"Outputing WF to file wf2d.out"

elseif(rank .eq. 3) then
  call normalize3d(wf3x,ngrid,dx)

  open(203,file='wf3d.out', action='WRITE', iostat=iost)
  write(203,*) "#WF - QDYN output"
  write(203,*) "#x  y  z   REAL   IMAG   PROBABILITY-DENSITY POTENTIAL"
  close(203)
  open(203,file='wf3d.out', status='old', position='append', action='WRITE', iostat=iost)

  call printwf3d(wf3x,x,y,z,v3)
  write(*,*)"Outputing WF to file wf3d.out"
endif

end subroutine init

subroutine check()

write(*,*) "====== Qdyn ====="

! ngrid is power of 2
if ((ngrid .ne. 0) .and. (IAND(ngrid, ngrid-1) .eq. 0))  then
  write(*,*) "Grid size: ",ngrid
else
  write(*,*) "ERR: Grid size must be power of two."
  stop 1
end if 

! params of grid
if ((rank .lt. 1) .or. (rank .gt. 3)) then
  write(*,*) "ERR: Dimensionality must be 1,2 or 3."
  stop 1
else
  write(*,*) "Number of dimensions: ",rank
end if

if (xmin .gt. xmax) then
  write(*,*) "ERR: xmin must be smaller than xmax."
  stop 1
else
  write(*,*) "xmin, xmax:", xmin, xmax
end if

select case (wf)
  case (0)
    write(*,*) "WF will be generated by program."
  case (1)
    write(*,*) "WF will be read from wf.in file."
    open(200,file='wf.in', status='OLD', action='READ',delim='APOSTROPHE', iostat=iost)
    if (iost.ne.0) then
      write(*,*)'ERROR: wf.in file must be provided'
      write(*,*) iost
      stop 1
    end if
  case default
    write(*,*) "ERR: wf must be either 0 or 1."
end select

if (pot == "") then
  write(*,*) "Potential not provided! Use analytical form. x,y,z for corresponding rank "
  stop 1
end if

if (run .eq. 0) then
!   write(*,*) "0 - IMAGINARY TIME PROPAGATION"
else
   write(*,*) "ERR: run must be set to 0"
   stop 1
end if

write(*,*) 'All checked.'

end subroutine check

end module

