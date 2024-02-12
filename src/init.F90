module mod_init
  use mod_vars
  use FFTW3
  use fparser,    ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg
  use mod_utils

  implicit none
  integer,private     :: file_unit, jstate
  character(len=50)   :: file_name
  real(DP), dimension(:), allocatable    :: point
  
CONTAINS

subroutine init()

write(*,*)
write(*,*) "### Initialization ###"

!-- Initialization of Wave Function
! initializing just for the 1. state in RT propagation
! initializing all states for IT propagation

!-- GRID set-up
dx=(xmax-xmin)/(ngrid-1)
dy=(ymax-ymin)/(ngrid-1)
dz=(zmax-zmin)/(ngrid-1)

!-- Allocating arrays
! grid
if(rank .eq. 1) allocate(x(ngrid),px(ngrid),point(1))
if(rank .eq. 2) allocate(x(ngrid),y(ngrid),px(ngrid),py(ngrid),point(2))
if(rank .eq. 3) allocate(x(ngrid),y(ngrid),z(ngrid),px(ngrid),py(ngrid),pz(ngrid),point(3))

! wf
if(rank .eq. 1) allocate(wfx(nstates,ngrid), wfp(ngrid))
if(rank .eq. 2) allocate(wf2x(nstates,ngrid,ngrid), wf2p(ngrid,ngrid))
if(rank .eq. 3) allocate(wf3x(nstates,ngrid,ngrid,ngrid), wf3p(ngrid,ngrid,ngrid))

! Hamiltonian
if(run .eq. 0) then
  if(rank .eq. 1) allocate(expV1(ngrid),v1_matrix(1,1,ngrid))
  if(rank .eq. 2) allocate(expV2(ngrid,ngrid),v2(ngrid,ngrid))
  if(rank .eq. 3) allocate(expV3(ngrid,ngrid,ngrid),v3(ngrid,ngrid,ngrid))

  if(field_coupling) then
    if(rank .eq. 1) allocate(dipole_coupling(nstates,nstates,ngrid))
  end if
else if (run .eq. 1) then
  if(rank .eq. 1) allocate(expV1(ngrid),v1(ngrid))
  if(rank .eq. 2) allocate(expV2(ngrid,ngrid),v2(ngrid,ngrid))
  if(rank .eq. 3) allocate(expV3(ngrid,ngrid,ngrid),v3(ngrid,ngrid,ngrid))
end if

if(rank .eq. 1) allocate(expT1(ngrid))
if(rank .eq. 2) allocate(expT2(ngrid,ngrid))
if(rank .eq. 3) allocate(expT3(ngrid,ngrid,ngrid))

!-- Setting up grid poitns
x(1)=xmin
if(rank .gt. 1) y(1)=ymin
if(rank .gt. 2) z(1)=zmin

do i=2, ngrid
  x(i) = x(i-1) + dx
  if(rank .gt. 1) y(i) = y(i-1) + dy
  if(rank .gt. 2) z(i) = z(i-1) + dy
end do

!-- POTENTIAL energy init

if (run.eq.1) then
  ! creating potential
  if (analytic) then
    write(*,*)"Potential: ",pot
    call initf (1)                                                     !Initialization of parser

    select case(rank)
    case(1)
      call parsef (1, pot, (/'x'/))                                       !Bytcompiling function string  
      do i=1, ngrid
        point = x(i)
        v1(i) = evalf (1, point)                                       !Evaluating potential for grid
        if (EvalErrType > 0) then
          WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
          stop 1
        end if
      end do

    case(2)
      call parsef (1, pot, (/'x','y'/))                                      
      do i=1, ngrid
        do j=1, ngrid
          point = (/x(i),y(j)/)
          v2(i,j) = evalf (1, point) 
          if (EvalErrType > 0) then
            WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
            stop 1
          end if
        end do
      end do

    case(3)
      call parsef (1, pot, (/'x','y','z'/))
      do i=1, ngrid
        do j=1, ngrid
          do k=1, ngrid
            point = (/x(i),y(j),z(k)/)
            v3(i,j,k) = evalf (1, point)
            if (EvalErrType > 0) then
              WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
              stop 1
            end if
          end do
        end do
      end do 
    end select
  ! reading potential from file pot.dat
  else 
    write(*,*) "Potential read from file: pot.dat"
    open(667,file='pot.dat', status='OLD', action='READ',delim='APOSTROPHE', iostat=iost)
    select case(rank)
    case(1)
      read(667,*)v1
    case(2)
      read(667,*)v2
    case(3)
      read(667,*)v3
    end select
    close(667)
  end if 

! IT: creating operator
  if(rank .eq. 1) then
    do i=1, ngrid
      expV1(i) = cmplx(exp(-v1(i)*dt/2.0d0),0)   !exp(-i V(x) tau/(2 h_bar))
    end do
  elseif(rank .eq. 2) then
    do i=1, ngrid
      do j=1, ngrid
        expV2(i,j) = cmplx(exp(-v2(i,j)*dt/2.0d0),0)  
      end do
    end do
  elseif(rank .eq. 3) then
    do i=1, ngrid
      do j=1, ngrid
        do k=1, ngrid
          expV3(i,j,k) = cmplx(exp(-v3(i,j,k)*dt/2.0d0),0)
        end do
      end do
    end do 
  end if

else if (run.eq.0) then
  ! creating potential
  if (analytic) then
    write(*,*)"Potential: ",pot
    call initf (1)                                                     !Initialization of parser

    select case(rank)
    case(1)
      call parsef (1, pot, (/'x'/))                                       !Bytcompiling function string  
      do i=1, ngrid
        point = x(i)
        v1_matrix(1,1,i) = evalf (1, point)                                       !Evaluating potential for grid
        if (EvalErrType > 0) then
          WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
          stop 1
        end if
      end do

    case(2)
      call parsef (1, pot, (/'x','y'/))                                      
      do i=1, ngrid
        do j=1, ngrid
          point = (/x(i),y(j)/)
          v2(i,j) = evalf (1, point) 
          if (EvalErrType > 0) then
            WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
            stop 1
          end if
        end do
      end do

    case(3)
      call parsef (1, pot, (/'x','y','z'/))
      do i=1, ngrid
        do j=1, ngrid
          do k=1, ngrid
            point = (/x(i),y(j),z(k)/)
            v3(i,j,k) = evalf (1, point)
            if (EvalErrType > 0) then
              WRITE(*,*)'*** Error evaluating potential: ',EvalErrMsg ()
              stop 1
            end if
          end do
        end do
      end do 
    end select
  ! reading potential from file pot.dat
  else 
    write(*,*) "Potential read from file: pot.dat"
    open(667,file='pot.dat', status='OLD', action='READ',delim='APOSTROPHE', iostat=iost)
    select case(rank)
    case(1)
      read(667,*)v1_matrix(1,1,:)
    case(2)
      read(667,*)v2
    case(3)
      read(667,*)v3
    end select
    close(667)
  end if 

  ! RT: creating operator
  if(rank .eq. 1) then
    do i=1, ngrid
      expV1(i) = cmplx(cos(-v1_matrix(1,1,i)*dt/2.0d0),sin(-v1_matrix(1,1,i)*dt/2.0d0))   !exp(-i V(x) tau/(2 h_bar))
    end do
  elseif(rank .eq. 2) then
    do i=1, ngrid
      do j=1, ngrid
        expV2(i,j) = cmplx(cos(-v2(i,j)*dt/2.0d0),sin(-v2(i,j)*dt/2.0d0))   !exp(-i V(x) tau/(2 h_bar))
      end do
    end do
  elseif(rank .eq. 3) then
    do i=1, ngrid
      do j=1, ngrid
        do k=1, ngrid
          expV3(i,j,k) = cmplx(cos(-v3(i,j,k)*dt/2.0d0),sin(-v3(i,j,k)*dt/2.0d0))  
        end do
      end do
    end do 
  end if
end if

! RT pot energy


!-- KINETIC energy init        exp[-iT/h_bar tau]

select case(rank)
  case(1)
    do i=1, ngrid
     if(i .le. ngrid/2) then
       px(i) = 2*pi*(i-1)/(ngrid*dx)
     else
       px(i) = 2*pi*(i-1-ngrid)/(ngrid*dx)
     end if
     if(run .eq. 0) expT1(i) = cmplx(dcos(-px(i)**2/(2*mass_x)*dt),dsin(-px(i)**2/(2*mass_x)*dt))
     if(run .eq. 1) expT1(i) = cmplx(dexp(-px(i)**2/(2*mass_x)*dt),0)
    end do
  !2D
  case(2)

   do i=1, ngrid
     if(i .le. ngrid/2) then
       px(i) = 2*pi*(i-1)/(ngrid*dx)
     else
       px(i) = 2*pi*(i-1-ngrid)/(ngrid*dx)
     end if
   end do

   do j=1, ngrid
     if(j .le. ngrid/2) then
       py(j) = 2*pi*(j-1)/(ngrid*dy)
     else
       py(j) = 2*pi*(j-1-ngrid)/(ngrid*dy)
     end if
   end do

    do i=1, ngrid
      do j=1, ngrid
        if(run .eq. 0) expT2(i,j) = cmplx(cos(-(px(i)**2/(2*mass_x)+py(j)**2/(2*mass_y))*dt),&
                                      sin(-(px(i)**2/(2*mass_x)+py(j)**2/(2*mass_y))*dt))
        if(run .eq. 1) expT2(i,j) = cmplx(exp(-(px(i)**2/(2*mass_x)+py(j)**2/(2*mass_y))*dt),0)
      end do
    end do

  !3D
  case(3)

   do i=1, ngrid
     if(i .le. ngrid/2) then
       px(i) = 2*pi*(i-1)/(ngrid*dx)
     else
       px(i) = 2*pi*(i-1-ngrid)/(ngrid*dx)
     end if
   end do

   do j=1, ngrid
     if(j .le. ngrid/2) then
       py(j) = 2*pi*(j-1)/(ngrid*dy)
     else
       py(j) = 2*pi*(j-1-ngrid)/(ngrid*dy)
     end if
   end do

   do k=1, ngrid
     if(k .le. ngrid/2) then
       pz(k) = 2*pi*(k-1)/(ngrid*dz)
     else
       pz(k) = 2*pi*(k-1-ngrid)/(ngrid*dz)
     end if
   end do
 
    do i=1, ngrid
      do j=1, ngrid
        do k=1, ngrid
          if(run .eq. 0) expT3(i,j,k) = cmplx(cos(-(px(i)**2/(2*mass_x)+py(j)**2/(2*mass_y)+pz(k)**2/(2*mass_z))*dt),&
                                         sin(-(px(i)**2/(2*mass_x)+py(j)**2/(2*mass_y)+pz(k)**2/(2*mass_z))*dt))
          if(run .eq. 1) expT3(i,j,k) = cmplx(exp(-(px(i)**2/(2*mass_x)+py(j)**2/(2*mass_y)+pz(k)**2/(2*mass_z))*dt),0)
        end do
      end do
    end do
end select

!-- Field preparations
if (field_coupling) then
  ! Initialize field function
  write(*,*) "Parsing field"
  call initf (2)                                                     !Initialization of parser
  call parsef (2, field, (/'t'/))                                       !Bytcompiling function string  
  if (EvalErrType > 0) then
    WRITE(*,*)'*** Error parsing field: ',EvalErrMsg ()
  end if

  ! Read dipole couplings
  write(*,*) "Reading dipole couplings (file names expected to be dipole_coup.$i.$j.dat)"
  call read_dipole_coupling()
end if


!--Generating wavepacket
call init_wavepacket()

!Normalize wf
if(rank .eq. 1) call normalize_1d(wfx(1,:)) 
if(rank .eq. 2) call normalize_2d(wf2x(1,:,:))
if(rank .eq. 3) call normalize_3d(wf3x(1,:,:,:)) 

!TODO: this is normalization just for the IT, where all initial states are normalized. I need to split initialization of RT and IT
if(run.eq.1 .and. nstates.ge.2) then
  do jstate=2,nstates
    if(rank .eq. 1) call normalize_1d(wfx(jstate,:))
    if(rank .eq. 2) call normalize_2d(wf2x(jstate,:,:))
    if(rank .eq. 3) call normalize_3d(wf3x(jstate,:,:,:)) 
  end do
end if

!--Printing of WF
if (print_wf) then
if(rank .eq. 1) then
  ! creating file name
  do jstate=1,nstates

    file_unit = 200+jstate
    write(file_name,*) jstate
    file_name='wf1d.'//trim(adjustl(file_name))//'.out'
    open(file_unit,file=file_name, action='WRITE', iostat=iost)
    ! opening file unit
    write(file_unit,*) "#x   REAL   IMAG   NORM    POTENTIAL"
    close(file_unit)
    open(file_unit,file=file_name, status='old', position='append', action='WRITE', iostat=iost)

    call printwf_1d(jstate)
    write(*,*)"Outputing WF to file "//file_name

  end do

elseif(rank .eq. 2) then
  ! creating file name
  do jstate=1,nstates

    file_unit = 200+jstate
    write(file_name,*) jstate
    file_name='wf2d.'//trim(adjustl(file_name))//'.out'
    open(file_unit,file=file_name, action='WRITE', iostat=iost)
    ! opening file unit
    write(file_unit,*) "#x  y   REAL   IMAG   NORM    POTENTIAL"
    close(file_unit)
    open(file_unit,file=file_name, status='old', position='append', action='WRITE', iostat=iost)

    call printwf_2d(jstate)
    write(*,*)"Outputing WF to file "//file_name

  end do

elseif(rank .eq. 3) then
  ! creating file name
  do jstate=1,nstates

    file_unit = 200+jstate
    write(file_name,*) jstate
    file_name='wf3d.'//trim(adjustl(file_name))//'.out'
    open(file_unit,file=file_name, action='WRITE', iostat=iost)
    ! opening file unit
    write(file_unit,*) "#x  y   z   REAL   IMAG   NORM    POTENTIAL"
    close(file_unit)
    open(file_unit,file=file_name, status='old', position='append', action='WRITE', iostat=iost)

    call printwf_3d(jstate)
    write(*,*)"Outputing WF to file wf3d.out"
  end do

end if
end if

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
    call update_energy_1d_rt()
  case(2)
    call update_energy_2d(wf2x(1,:,:))
  case(3)
    call update_energy_3d(wf3x(1,:,:,:))
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
      call update_energy_2d(wf2x(1,:,:))
    case(3)
      call update_energy_3d(wf3x(1,:,:,:))
    end select
    
    call printen_state(jstate)
  end do
end select

!--Open file with field
if (field_coupling) then
  file_unit = 102
  file_name = 'field.dat'
  open(file_unit,file=file_name, action='WRITE', iostat=iost)
  write(file_unit,*) "#  time     field amplitude"
  close(file_unit)
  open(file_unit,file=file_name, status='old', position='append', action='WRITE', iostat=iost)

  call print_field()

end if

!closing input.q
close(100)

write(*,*) 

end subroutine init

subroutine init_wavepacket()

  real(DP)    :: x0, y0, z0, xsigma, ysigma, zsigma, px0, py0, pz0
  real(DP)    :: prefactor, gauss, momenta, norm

  namelist /init_wf/x0,y0,z0,xsigma,ysigma,zsigma,px0,py0,pz0

if(wf .eq. 0) then                                           
  write(*,*) "Generating initial gaussian wave packet using section &init_wf."

  ! Default values
  x0 = (xmax-xmin)/2.0d0+xmin+(xmax-xmin)/10                      ! shifted a little bit so that the wave packet is not symmetric
  y0 = (ymax-ymin)/2.0d0+ymin+(ymax-ymin)/10                      ! shifted a little bit so that the wave packet is not symmetric
  z0 = (zmax-zmin)/2.0d0+zmin+(zmax-zmin)/10                      ! shifted a little bit so that the wave packet is not symmetric
  xsigma = (xmax-x0)/20.0d0                                       ! 5sigma - 96% of gaussian is on the grid 
  ysigma = (ymax-y0)/20.0d0                                       ! 5sigma - 96% of gaussian is on the grid 
  zsigma = (zmax-z0)/20.0d0                                       ! 5sigma - 96% of gaussian is on the grid 
  px0 = 0.0d0
  py0 = 0.0d0
  pz0 = 0.0d0

  ! reading init_wf section in the input.q
  rewind(100)
  read(100, init_wf, iostat=iost)
  if (iost.ne.0) then
    write(*,*)'ERROR: &init_wf section must be provided in input.q'
    write(*,*) iost
    stop 1
  end if
  close(100)

  write(*,'(A,F10.5)') " x0 = ", x0
  write(*,'(A,F10.5)') " y0 = ", y0
  write(*,'(A,F10.5)') " z0 = ", z0
  write(*,'(A,F10.5)') " xsigma = ", xsigma
  write(*,'(A,F10.5)') " ysigma = ", ysigma
  write(*,'(A,F10.5)') " zsigma = ", zsigma
  write(*,'(A,F10.5)') " px0 = ", px0
  write(*,'(A,F10.5)') " py0 = ", py0
  write(*,'(A,F10.5)') " pz0 = ", pz0
  
  ! normalization prefactor of the gaussian, same value for the whole grid
  select case (rank)
  case (1)
    prefactor = xsigma**(-0.5)*pi**(-0.25)
  case (2)
    prefactor = (xsigma*ysigma)**(-0.5)*pi**(-0.5)
  case (3)
    prefactor = (xsigma*ysigma*zsigma)**(-0.5)*pi**(-0.75)
  end select

  do i=1, ngrid
    select case (rank)
      case (1)
        ! Gaussian part of the wave packet depending only on positions
        gauss = prefactor * exp(-(x(i)-x0)**2/(2*xsigma**2))
        ! momentum calculation p0*(x-x0)
        momenta = px0*(x(i)-x0)
        ! whole imaginary wf
        wfx(1,i) = cmplx(gauss * cos(momenta), gauss * sin(momenta))  
        ! for the imag time, I create the came wave packet for all states
        if (run.eq.1 .and. nstates.ge.2) then
          do jstate=2,nstates
            wfx(jstate,i) = wfx(1,i)
          end do
        end if
            
      case (2)
        do j=1, ngrid
          ! Gaussian part of the wave packet depending only on positions
          gauss = prefactor * exp(-(x(i)-x0)**2/(2*xsigma**2) -(y(j)-y0)**2/(2*ysigma**2))
          ! momentum calculation p0*(x-x0)
          momenta = px0*(x(i)-x0) + py0*(y(j)-y0)
          ! whole imaginary wf
          wf2x(1,i,j) = cmplx(gauss * cos(momenta), gauss * sin(momenta))

        ! for the imag time, I create the came wave packet for all states
        if (run.eq.1 .and. nstates.ge.2) then
          do jstate=2,nstates
            wf2x(jstate,i,j) = wf2x(1,i,j)
          end do
        end if

        end do
            
      case (3)
        do j=1, ngrid
          do k=1, ngrid
            ! Gaussian part of the wave packet depending only on positions
            gauss = prefactor * exp(-(x(i)-x0)**2/(2*xsigma**2)-(y(j)-y0)**2/(2*ysigma**2)-(z(k)-z0)**2/(2*zsigma**2))
            ! momentum calculation p0*(x-x0)
            momenta = px0*(x(i)-x0) + py0*(y(j)-y0) + pz0*(z(k)-z0)
            ! whole imaginary wf
            wf3x(1,i,j,k) = cmplx(gauss * cos(momenta), gauss * sin(momenta)) 

          ! for the imag time, I create the came wave packet for all states
          if (run.eq.1 .and. nstates.ge.2) then
            do jstate=2,nstates
              wf3x(jstate,i,j,k) = wf3x(1,i,j,k)
            end do
          end if

          end do
        end do
    end select
  end do

  ! printing norm of generated wf, it is probably point less since I normalize later
  ! but just to check we have a correct initial guess
  select case(rank)
  case(1)
    norm = braket_1d(wfx(1,:), wfx(1,:))
  case(2)
    norm = braket_2d(wf2x(1,:,:), wf2x(1,:,:))
  case(3)
    norm = braket_3d(wf3x(1,:,:,:), wf3x(1,:,:,:))
  end select
  write(*,'(A,F10.5)') " norm = ", norm 

elseif(wf .eq. 1) then
  !Procedure for loading WF from file
  read(666,*) !two empty lines
  read(666,*)
  if(rank .eq. 1) read(666,*)wfx
  if(rank .eq. 2) read(666,*)wf2x
  if(rank .eq. 3) read(666,*)wf3x
  close(666)
end if

end subroutine

subroutine read_dipole_coupling()

  integer     :: istate, jstate
  logical     :: file_exists

  do istate=1,nstates
    do jstate=1,nstates

      write(file_name,'(I1,A,I1)') istate,".",jstate
      file_name='dipole_coup.'//trim(adjustl(file_name))//'.dat'
      inquire(file=file_name, exist=file_exists)

      if (file_exists) then
        file_unit = 400
        open(file_unit,file=file_name, status='OLD', action='READ',delim='APOSTROPHE', iostat=iost)
        read(file_unit,*) dipole_coupling(istate,jstate,:)
        close(file_unit)

        write(*,'(a,i1,a,i1,a)') " * <psi_",istate,"|mu|psi_", jstate,"> = read from file"
        print *, dipole_coupling(istate,jstate,:)
        !TODO: save data to file
      else
        write(*,'(a,i1,a,i1,a)') " * <psi_",istate,"|mu|psi_", jstate,"> = 0"
        !TODO: equal it zero
        dipole_coupling(istate,jstate,:) = 0.0d0
        print *, dipole_coupling(istate,jstate,:)
      end if
    end do
  end do

end subroutine

end module

