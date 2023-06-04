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
read(666,*) !two empty lines
read(666,*)
if(rank .eq. 1) read(666,*)wfx
if(rank .eq. 2) read(666,*)wf2x
if(rank .eq. 3) read(666,*)wf3x
close(666)
end if

!>jj
!Reading ground state wf
if (nstates.eq.2) then
open(667,file='wfgs.chk', status='OLD', action='READ',delim='APOSTROPHE', iostat=iost)
if (iost.ne.0) then
  write(*,*)'ERROR: wf.chk file must be provided'
  write(*,*) iost
  stop 1
end if
read(667,*) !two empty lines
read(667,*)
if(rank .eq. 1) read(667,*)wfxgs
close(667)
write(*,*)"*GS wf read"
end if
!<jj

!Normalize wf
if(rank .eq. 1) call normalize_1d(wfx) 
if(rank .eq. 2) call normalize_2d(wf2x)
if(rank .eq. 3) call normalize_3d(wf3x) 

!-- Initialization of FFT procedures
if(rank .eq. 1) then

call dfftw_plan_dft_1d(plan_forward, ngrid, wfx, wfp, FFTW_FORWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_forward, wfx, wfp)
call dfftw_destroy_plan(plan_forward)
wfp = wfp / dsqrt(real(ngrid, kind=DP))

call dfftw_plan_dft_1d(plan_backward, ngrid, wfp, wfx, FFTW_BACKWARD, FFTW_ESTIMATE )
call dfftw_execute_dft(plan_backward, wfp, wfx)
call dfftw_destroy_plan(plan_backward)
wfx = wfx / dsqrt(real(ngrid, kind=DP))

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
  call normalize_1d(wfx)

  ! creating file name
  write(file_name,*) nstates
  file_name='wf1d.'//trim(adjustl(file_name))//'.out'
  write(*,*) file_name
  open(201,file=file_name, action='WRITE', iostat=iost)
  ! opening file unit
  write(201,*) "#WF - QDYN output"
  write(201,*) "#x   REAL   IMAG   PROBABILITY-DENSITY POTENTIAL"
  close(201)
  open(201,file=file_name, status='old', position='append', action='WRITE', iostat=iost)

  call printwf_1d(wfx,x,v1)
  write(*,*)"Outputing WF to file "//file_name

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
open(101,file='energies.dat', action='WRITE', iostat=iost)
write(101,*) "# time    energy"
close(101)
open(101,file='energies.dat', status='old', position='append', action='WRITE', iostat=iost)

!--Writing energies
select case(rank)
case(1)
  call update_energy_1d(wfx, energy)
  call printen(time, energy)
case(2)
  call update_energy_2d(wf2x, energy)
  call printen(time, energy)
case(3)
  call update_energy_3d(wf3x, energy)
  call printen(time, energy)
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

