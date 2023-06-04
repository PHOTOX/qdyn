module mod_vars

use fparser,    ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg

  implicit none
  public
  INTEGER, PARAMETER    :: DP = KIND(1.0d0)
  real(DP), parameter   :: pi = 3.14159265
  real(DP)              :: dt, xmin, xmax, xmean, stddev, dx, k_0, mass, dtwrite, energy
  real(DP), dimension(:), allocatable    :: x,y,z,point
  complex(DP), dimension(:), allocatable :: wfx, wfp, theta_v1, kin_p1
  !jj - add up
  complex(DP), dimension(:), allocatable :: wfxgs
  complex(DP), dimension(:,:), allocatable :: wf2x, wf2p, theta_v2, kin_p2
  complex(DP), dimension(:,:,:), allocatable :: wf3x, wf3p, theta_v3, kin_p3
  integer               :: run, nstep, ngrid, wf, rank, iost, i, j, k, nstates
  integer ( kind = 8 )  :: plan_forward, plan_backward

  real(DP), dimension(:), allocatable     :: v1, px, py, pz
  real(DP), dimension(:,:), allocatable   :: v2
  real(DP), dimension(:,:,:), allocatable :: v3
  character(len=50)             :: pot=''
  !jj
  character(len=50)             :: file_name
  character(len=*),dimension(1),parameter :: var1 = (/'x'/)
  character(len=*),dimension(2),parameter :: var2 = (/'x','y'/)
  character(len=*),dimension(3),parameter :: var3 = (/'x','y','z'/)

  namelist /general/run,nstep,dt,dtwrite,ngrid,rank,xmin,xmax,mass,wf,pot,nstates

CONTAINS

subroutine read_input()
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
!jj
if(rank .eq. 1) allocate(wfxgs(ngrid))
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
      if(run .eq. 0) theta_v1(i) = cmplx(cos(-v1(i)*dt/2.0d0),sin(-v1(i)*dt/2.0d0))   !exp(-i V(x) tau/(2 h_bar))
      if(run .eq. 1) theta_v1(i) = cmplx(exp(-v1(i)*dt/2.0d0),exp(-v1(i)*dt/2.0d0))   !exp(-i V(x) tau/(2 h_bar))
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
        if(run .eq. 0) theta_v2(i,j) = cmplx(cos(-v2(i,j)*dt/2.0d0),sin(-v2(i,j)*dt/2.0d0))   !exp(-i V(x) tau/(2 h_bar))
        if(run .eq. 1) theta_v2(i,j) = cmplx(exp(-v2(i,j)*dt/2.0d0),exp(-v2(i,j)*dt/2.0d0))  
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
          if(run .eq. 0) theta_v3(i,j,k) = cmplx(cos(-v3(i,j,k)*dt/2.0d0),sin(-v3(i,j,k)*dt/2.0d0))  
          if(run .eq. 1) theta_v3(i,j,k) = cmplx(exp(-v3(i,j,k)*dt/2.0d0),exp(-v3(i,j,k)*dt/2.0d0))
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
     if(run .eq. 0) kin_p1(i) = cmplx(cos(-px(i)**2*dt/(2*mass)),sin(-px(i)**2*dt/(2*mass)))
     if(run .eq. 1) kin_p1(i) = cmplx(exp(-px(i)**2*dt/(2*mass)),exp(-px(i)**2*dt/(2*mass)))
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
        if(run .eq. 0) kin_p2(i,j) = cmplx(cos(-(px(i)**2+py(j)**2)*dt/(2*mass)),sin(-(px(i)**2+py(j)**2)*dt/(2*mass)))
        if(run .eq. 1) kin_p2(i,j) = cmplx(exp(-(px(i)**2+py(j)**2)*dt/(2*mass)),exp(-(px(i)**2+py(j)**2)*dt/(2*mass)))
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
          if(run .eq. 0) kin_p3(i,j,k) = cmplx(cos(-(px(i)**2+py(j)**2+pz(k)**2)*dt/(2*mass)),&
                                         sin(-(px(i)**2+py(j)**2+pz(k)**2)*dt/(2*mass)))
          if(run .eq. 1) kin_p3(i,j,k) = cmplx(exp(-(px(i)**2+py(j)**2+pz(k)**2)*dt/(2*mass)),&
                                         exp(-(px(i)**2+py(j)**2+pz(k)**2)*dt/(2*mass)))
        end do
      end do
    end do
end select

!--Printing of WF

end subroutine read_input

subroutine check()

write(*,*) "====== Qdyn ====="

! ngrid is power of 2
if ((ngrid .ne. 0) .and. (IAND(ngrid, ngrid-1) .eq. 0))  then
  write(*,'(A,I5)') " Grid size: ",ngrid
else
  write(*,*) "ERR: Grid size must be power of two."
  stop 1
end if 

! dimensionality
if ((rank .lt. 1) .or. (rank .gt. 3)) then
  write(*,*) "ERR: Dimensionality must be 1,2 or 3."
  stop 1
else
  write(*,'(A,I1)') " Number of dimensions: ",rank
end if

! params of grid
if (xmin .gt. xmax) then
  write(*,*) "ERR: xmin must be smaller than xmax."
  stop 1
else
  write(*,'(A,F8.4,F8.4)') " xmin, xmax: ", xmin, xmax
end if

! initial wf selection
select case (wf)
  case (0)
    write(*,*) "WF will be generated by program."
  case (1)
    write(*,*) "WF will be read from wf.chk file."
    open(666,file='wf.chk', status='OLD', action='READ',delim='APOSTROPHE', iostat=iost)
    if (iost.ne.0) then
      write(*,*)'ERROR: wf.chk file must be provided'
      write(*,*) iost
      stop 1
    end if
  case default
    write(*,*) "ERR: wf must be either 0 or 1."
end select

! potential
if (pot == "") then
  write(*,*) "Potential not provided! Use analytical form. x,y,z for corresponding rank "
  stop 1
end if

! run case (imag/real)
select case(run)
  ! CLASSICAL REAL TIME PROPAGATION
  case(0)
    write(*,*) "Propagation: REAL TIME"
  ! CLASSICAL IMAGINARY TIME PROPAGATION
  case(1)
    write(*,*) "Propagation: IMAGINARY TIME"
  case default
    write(*,*) "ERR: Unrecongnized run option. Exiting"
    stop 1
end select

! number of states
if (nstates < 1) then
  write(*,*) "ERR: number of states must be 1 or more."
  stop 1
else
  if (run .eq. 0) then
    write(*,*) "ERR: nstates > 1 available only for imag propagation."
    stop 1
  else
    write(*,'(A,I1)') " nstates: ", nstates
  end if
end if

!>jj because code is not ready for more states now
if (nstates > 2) then
  write(*,*) "Too many states!"
  stop 1
else
  if (nstates.eq.2 .and. rank.gt.1) then
        write(*,*) "CODE NOT READY FOR IMAG PROP WITH MORE THAN TWO STATES!"
        stop 1
  end if
end if
!<jj

end subroutine check

end module

