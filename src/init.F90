module mod_init
  use mod_vars
  use FFTW3
  use fparser,    ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg
  use mod_utils

  implicit none
  public

CONTAINS

subroutine init()
  implicit none

!-- Initialization of Wave Function
! initializing just for the 1. state in RT propagation
! initializing all states for IT propagation



!!!!

!-- GRID set-up
dx=(xmax-xmin)/(ngrid-1)

! allocating arrays
if(rank .eq. 1) allocate(x(ngrid),v1(ngrid),px(ngrid),point(1))
if(rank .eq. 2) allocate(x(ngrid),y(ngrid),v2(ngrid,ngrid),px(ngrid),py(ngrid),point(2))
if(rank .eq. 3) allocate(x(ngrid),y(ngrid),z(ngrid),v3(ngrid,ngrid,ngrid),px(ngrid),py(ngrid),pz(ngrid),point(3))

!if(rank .eq. 1) allocate(wfx(ngrid), wfp(ngrid), theta_v1(ngrid), kin_p1(ngrid))
!jj
if(rank .eq. 1) allocate(wfx(nstates,ngrid), wfp(ngrid), theta_v1(ngrid), kin_p1(ngrid))
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

!TODO: move this to init.F90
!-- POTENTIAL energy init
if (analytic) then
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
      if(run .eq. 1) theta_v1(i) = cmplx(exp(-v1(i)*dt/2.0d0),0)   !exp(-i V(x) tau/(2 h_bar))
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
        if(run .eq. 1) theta_v2(i,j) = cmplx(exp(-v2(i,j)*dt/2.0d0),0)  
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
          if(run .eq. 1) theta_v3(i,j,k) = cmplx(exp(-v3(i,j,k)*dt/2.0d0),0)
        end do
      end do
    end do 
  end if
! reading potential from file pot.dat
else 
  write(*,*) "Potential read from file: pot.dat"
  !call interpolation_1d(...)
  open(667,file='pot.dat', status='OLD', action='READ',delim='APOSTROPHE', iostat=iost)
  select case(rank)
  case(1)
    read(667,*)v1
    do j=1, ngrid
      if(run .eq. 0) theta_v1(j) = cmplx(cos(-v1(j)*dt/2.0d0),sin(-v1(j)*dt/2.0d0))   !exp(-i V(x) tau/(2 h_bar))
      if(run .eq. 1) theta_v1(j) = cmplx(exp(-v1(j)*dt/2.0d0),0)   !exp(-i V(x) tau/(2 h_bar))
    end do
  case default
    write(*,*) "cannot read files with rank bigger than 1."
    stop 1
  end select
  close(667)
end if 

!-- KINETIC energy init        exp[-iT/h_bar tau]

select case(rank)
  case(1)
    do i=1, ngrid
!     if(i .lt. ngrid/2) then
     if(i .le. ngrid/2) then
!       px(i) = 2*pi*i/(ngrid*dx)
       !jj
       px(i) = 2*pi*(i-1)/(ngrid*dx)
     else
!       px(i) = 2*pi*(i-ngrid)/(ngrid*dx)
       !jj
       px(i) = 2*pi*(i-1-ngrid)/(ngrid*dx)
     end if
     if(run .eq. 0) kin_p1(i) = cmplx(dcos(-px(i)**2*dt/(2*mass)),dsin(-px(i)**2*dt/(2*mass)))
     if(run .eq. 1) kin_p1(i) = cmplx(dexp(-px(i)**2*dt/(2*mass)),0)
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
        if(run .eq. 1) kin_p2(i,j) = cmplx(exp(-(px(i)**2+py(j)**2)*dt/(2*mass)),0)
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
          if(run .eq. 1) kin_p3(i,j,k) = cmplx(exp(-(px(i)**2+py(j)**2+pz(k)**2)*dt/(2*mass)),0)
        end do
      end do
    end do
end select

!--Initialize field function
if (use_field) then
  call initf (2)                                                     !Initialization of parser
  call parsef (2, field, (/'t'/))                                       !Bytcompiling function string  
  if (EvalErrType > 0) then
    WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
  end if
end if

!--Generating wavepacket
if(wf .eq. 0) then                                                           ! generating gaussian wavepacket
  write(*,*) "Generating gaussian wave packet at the center of the grid."
  !TODO: make these variables modifiable in input
  xmean=(xmax-xmin)/2.0d0+xmin*1.1d0 ! on purpose a bit shifted
  stddev=(xmax-xmean)/20.0d0                                                  ! 5sigma - 96% of gaussian is on the grid 
  k_0 = sqrt(2*mass*0.5)                                                   ! sqrt(2*m*E)/h = k0
  !jj - No initial momentum set for the wavepacket
  k_0 = k_0*0.0d0
  do i=1, ngrid
    select case (rank)
      case (1)
        !jj - for imag propagation, I should loop jstate and copy the value
        wfx(1,i) = cmplx(exp((-1.0d0*(x(i)-xmean)**2)/(2*stddev**2)) * cos(k_0 * x(i)), &
                       exp((-1.0d0*(x(i)-xmean)**2)/(2*stddev**2)) * sin(k_0 * x(i)) )  
        ! for the imag time, I create the came wave packet for all states
        if (run.eq.1 .and. nstates.ge.2) then
          do jstate=2,nstates
            wfx(jstate,i) = wfx(1,i)
          end do
        end if
            
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
read(666,*) !two empty lines
read(666,*)
if(rank .eq. 1) read(666,*)wfx
if(rank .eq. 2) read(666,*)wf2x
if(rank .eq. 3) read(666,*)wf3x
close(666)
end if

!Normalize wf
if(rank .eq. 1) call normalize_1d(wfx(1,:)) 
!>jj - normalization of the other states
! I should do this in the normalize 
if(run.eq.1 .and. nstates.ge.2) then
  do jstate=2,nstates
    call normalize_1d(wfx(jstate,:))
  end do
end if
!<jj
if(rank .eq. 2) call normalize_2d(wf2x)
if(rank .eq. 3) call normalize_3d(wf3x) 

!-- Initialization of FFT procedures
if(rank .eq. 1) then

call dfftw_plan_dft_1d(plan_forward, ngrid, wfx(1,:), wfp, FFTW_FORWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_forward, wfx(1,:), wfp)
call dfftw_destroy_plan(plan_forward)
wfp = wfp / dsqrt(real(ngrid, kind=DP))

call dfftw_plan_dft_1d(plan_backward, ngrid, wfp, wfx(1,:), FFTW_BACKWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_backward, wfp, wfx(1,:))
call dfftw_destroy_plan(plan_backward)
wfx(1,:) = wfx(1,:) / dsqrt(real(ngrid, kind=DP))

elseif(rank .eq. 2) then

call dfftw_plan_dft_2d(plan_forward, ngrid, ngrid, wf2x, wf2p, FFTW_FORWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_forward, wf2x, wf2p)
call dfftw_destroy_plan(plan_forward)
wf2p = wf2p / dsqrt(real(ngrid, kind=DP)**2)

call dfftw_plan_dft_2d(plan_backward, ngrid, ngrid, wf2p, wf2x, FFTW_BACKWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_backward, wf2p, wf2x)
call dfftw_destroy_plan(plan_backward)
wf2x = wf2x / dsqrt(real(ngrid, kind=DP)**2)

elseif(rank .eq. 3) then

call dfftw_plan_dft_3d(plan_forward, ngrid, ngrid, ngrid, wf3x, wf3p, FFTW_FORWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_forward, wf3x, wf3p)
call dfftw_destroy_plan(plan_forward)
wf3p = wf3p / dsqrt(real(ngrid, kind=DP)**3)

call dfftw_plan_dft_3d(plan_backward, ngrid, ngrid, ngrid, wf3p, wf3x, FFTW_BACKWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_backward, wf3p, wf3x)
call dfftw_destroy_plan(plan_backward)
wf3x = wf3x / dsqrt(real(ngrid, kind=DP)**3)

end if

!--Printing of WF

if(rank .eq. 1) then
  call normalize_1d(wfx(1,:))

  ! creating file name
  do jstate=1,nstates

    file_unit = 200+jstate
    write(file_name,*) jstate
    file_name='wf1d.'//trim(adjustl(file_name))//'.out'
    open(file_unit,file=file_name, action='WRITE', iostat=iost)
    ! opening file unii
    write(file_unit,*) "#WF - QDYN output"
    write(file_unit,*) "#x   REAL   IMAG   PROBABILITY-DENSITY POTENTIAL"
    close(file_unit)
    open(file_unit,file=file_name, status='old', position='append', action='WRITE', iostat=iost)

    call printwf_1d(jstate,x,v1)
    write(*,*)"Outputing WF to file "//file_name

  end do

elseif(rank .eq. 2) then
  call normalize_2d(wf2x)

  open(202,file='wf2d.out', action='WRITE', iostat=iost)
  write(202,*) "#WF - QDYN output"
  write(202,*) "#x  y   REAL   IMAG   PROBABILITY-DENSITY  POTENTIAL"
  close(202)
  open(202,file='wf2d.out', status='old', position='append', action='WRITE', iostat=iost)

  call printwf_2d(wf2x,x,y,v2)
  write(*,*)"Outputing WF to file wf2d.out"

elseif(rank .eq. 3) then
  call normalize_3d(wf3x)

  open(203,file='wf3d.out', action='WRITE', iostat=iost)
  write(203,*) "#WF - QDYN output"
  write(203,*) "#x  y  z   REAL   IMAG   PROBABILITY-DENSITY POTENTIAL"
  close(203)
  open(203,file='wf3d.out', status='old', position='append', action='WRITE', iostat=iost)

  call printwf_3d(wf3x,x,y,z,v3)
  write(*,*)"Outputing WF to file wf3d.out"
endif

!--Open file with energies
select case(run)
case(0)
  open(101,file='energies.dat', action='WRITE', iostat=iost)
  write(101,*) "#  time     total energy    potential       kinetic         energy diff     norm"
  close(101)
  open(101,file='energies.dat', status='old', position='append', action='WRITE', iostat=iost)

  !--Writing energies
  select case(rank)
  case(1)
    call update_energy_1d(wfx(1,:))
  case(2)
    call update_energy_2d(wf2x, energy(1))
  case(3)
    call update_energy_3d(wf3x, energy(1))
  end select
  call update_norm()
  call printen()
case(1)
  do jstate=1,nstates
    file_unit=300+jstate
    write(file_name,*) jstate
    file_name='energies.'//trim(adjustl(file_name))//'.dat'
    open(file_unit,file=file_name, action='WRITE', iostat=iost)
    write(file_unit,*) "#  time     total energy    potential       kinetic         energy diff"
    close(file_unit)
    open(file_unit,file=file_name, status='old', position='append', action='WRITE', iostat=iost)

    !--Writing energies
    select case(rank)
    case(1)
      call update_energy_1d(wfx(1,:))
    case(2)
      call update_energy_2d(wf2x, energy(1))
    case(3)
      call update_energy_3d(wf3x, energy(1))
    end select
    
    call printen_state(jstate)
  end do
end select

! printing beggining of the output
write(*,*) 

select case(run)
  case(0)
    write(*,*) "RUN: 0 - REAL TIME PROPAGATION"
  case(1)
    write(*,*) "RUN: 1 - IMAGINARY TIME PROPAGATION"
  case default
    stop 1
end select

write(*,*)
write(*,*) "------ time -----"

end subroutine init

end module

