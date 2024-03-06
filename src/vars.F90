module mod_vars

  implicit none
  public
  !-parameters
  INTEGER, PARAMETER    :: DP = KIND(1.0d0)
  real(DP), parameter   :: pi = 3.14159265358979323846
  !-propagation data
  integer               :: run, nstep, xngrid, yngrid, zngrid, rank, nstates=1
  real(DP)              :: xmin, xmax, dx, ymin, ymax, dy, zmin, zmax, dz, dtwrite, dt 
  real(DP)              :: mass_x = 0.0, mass_y = 0.0, mass_z = 0.0
  real(DP)              :: time = 0.0, norm, energy(3), energy_diff ! energy(total, potential, kinetic)
  real(DP)              :: norm_thresh = 1.0d-1
  real(DP), dimension(:), allocatable    :: x, y, z, px, py, pz 
  real(DP), dimension(:), allocatable    :: diab_pop, ad_pop
  logical               :: project_rot=.true., analytic=.true., print_wf=.true.
  character(len=10)     :: dynamics=''
  !-wave function
  complex(DP), dimension(:), allocatable :: wfp
  complex(DP), dimension(:,:), allocatable :: wfx, wfx_ad, wf2p
  complex(DP), dimension(:,:,:), allocatable :: wf2x, wf3p
  complex(DP), dimension(:,:,:,:), allocatable :: wf3x
  !-field
  logical               :: field_coupling=.false.
  character(len=100)    :: field=''
  !-auxiliary variables
  integer               :: iost, i, j, k
  integer ( kind = 8 )  :: plan_forward, plan_backward
  character(len=100)    :: pot=''
  !-exp(T) operator
  complex(DP), dimension(:), allocatable :: expT1
  complex(DP), dimension(:,:), allocatable :: expT2
  complex(DP), dimension(:,:,:), allocatable :: expT3
  !-imaginary time V operators (V is used for IT, H is used for RT)
  real(DP), dimension(:), allocatable     :: v1
  real(DP), dimension(:,:), allocatable   :: v2
  real(DP), dimension(:,:,:), allocatable :: v3
  complex(DP), dimension(:), allocatable :: expV1 !exp(V)
  complex(DP), dimension(:,:), allocatable :: expV2
  complex(DP), dimension(:,:,:), allocatable :: expV3
  !-real time H operators for diabatic H_el hamiltonian (V is used for IT, H is used for RT)
  ! TODO: these are prepared variables for RT propagation, not allocated yet in the init()
  real(DP), dimension(:,:), allocatable       :: H1_ad
  real(DP), dimension(:,:,:), allocatable     :: H1, U1, invU1, dipole_coupling  ! U = from adiabatic to diabatic, invU = from diabatic to adiabatic 
! real(DP), dimension(:,:,:), allocatable   :: H2
! real(DP), dimension(:,:,:,:), allocatable :: H3
  complex(DP), dimension(:,:,:), allocatable      :: expH1
! complex(DP), dimension(:,:,:,:), allocatable    :: expH2_matrix
! complex(DP), dimension(:,:,:,:,:), allocatable  :: expH3_matrix

  namelist /general/ dynamics, nstep, dt, dtwrite, xngrid, yngrid, zngrid, rank, &
    xmin, xmax, ymin, ymax, zmin, zmax, mass_x, mass_y, mass_z, nstates, print_wf
  namelist /it/ pot, analytic, project_rot
  namelist /rt/ pot, analytic, field_coupling, field, norm_thresh
  !TODO: rt analytic will not be use probably

CONTAINS

subroutine read_input()
  implicit none

  write(*,*)
  write(*,*) "### Reading input ###"

!-- Reading input file
  open(100,file='input.q', status='OLD', action='READ',delim='APOSTROPHE', iostat=iost)
  if (iost.ne.0) then
    write(*,*)'ERROR: input.q file must be provided'
    write(*,*) iost
    stop 1
  end if

  read(100, general, iostat=iost)
  if (iost.ne.0) then
    write(*,*)'ERROR: &general section missing or problematic'
    write(*,*) iost
    stop 1
  end if

  rewind(100)
  select case(dynamics)
  case('rt')
    read(100, rt, iostat=iost)
    if (iost.ne.0) then
      write(*,*)'ERROR: &rt section missing or problematic'
      write(*,*) iost
      stop 1
    end if
  case('it')
    read(100, it, iostat=iost)
    if (iost.ne.0) then
      write(*,*)'ERROR: &it section missing or problematic'
      write(*,*) iost
      stop 1
    end if
  end select

  ! file is closed at the end of init() subroutine because some more input is read. 

!-- Input check
  call check()


end subroutine read_input

subroutine check()

! run case (imag/real)
select case(dynamics)
  case('rt')
    write(*,*) "RUN: REAL TIME PROPAGATION"
    run = 0
  case('it')
    write(*,*) "RUN: IMAGINARY TIME PROPAGATION"
    run = 1
  case default
    write(*,*) "ERR: Unrecongnized 'dynamics' option. Choose either 'rt' or 'it'. Exiting"
    stop 1
end select

! writing time step and number of steps
write(*,'(A,F9.5,A)') " Time step: ", dt, " a.u."
write(*,'(A,I8)') " Number of steps: ", nstep 
write(*,'(A,F10.2,A)') " Total time: ", dt*nstep, " a.u."

! ngrid is power of 2
write(*,'(A,I5)') " Grid size in x: ",xngrid
if ((xngrid .le. 0) .or. (IAND(xngrid, xngrid-1) .ne. 0))  then
  write(*,*) "ERR: Grid size must be power of two."
  stop 1
end if 

if (rank .ge. 2) then
  write(*,'(A,I5)') " Grid size in y: ",yngrid
  if ((yngrid .le. 0) .or. (IAND(yngrid, yngrid-1) .ne. 0))  then
    write(*,*) "ERR: Grid size must be power of two."
    stop 1
  end if 
end if

if (rank .ge. 3) then
  write(*,'(A,I5)') " Grid size in z: ",zngrid
  if ((zngrid .le. 0) .or. (IAND(zngrid, zngrid-1) .ne. 0))  then
    write(*,*) "ERR: Grid size must be power of two."
    stop 1
  end if 
end if

! dimensionality
if ((rank .lt. 1) .or. (rank .gt. 3)) then
  write(*,*) "ERR: Dimensionality must be 1,2 or 3."
  stop 1
else
  write(*,'(A,I1)') " Number of dimensions: ",rank
end if

! number of steps
if (nstep .lt. 1) then
  write(*,*) "ERR: Number of steps must be bigger than 1."
  stop 1
else
  write(*,'(A,I8)') " Number of steps: ",rank
end if

! masses
if (mass_x .le. 0.0) then
  write(*,*) "ERR: mass_x was either not set or set negative."
  stop 1
end if

if ((rank .ge. 2) .and. (mass_y .le. 0.0)) then
  write(*,*) "ERR: mass_y was either not set or set negative."
  stop 1
end if

if ((rank .ge. 3) .and. (mass_z .le. 0.0)) then
  write(*,*) "ERR: mass_z was either not set or set negative."
  stop 1
end if

! params of grid
write(*,'(A,F8.4,F8.4)') " xmin, xmax: ", xmin, xmax
if (xmin .ge. xmax) then
  write(*,*) "ERR: xmin must be smaller than xmax."
  stop 1
end if

if (rank .ge. 2) then
write(*,'(A,F8.4,F8.4)') " ymin, ymax: ", ymin, ymax
  if (ymin .ge. ymax) then
    write(*,*) "ERR: ymin must be smaller than ymax."
    stop 1
  end if
end if

if (rank .ge. 3) then
write(*,'(A,F8.4,F8.4)') " zmin, zmax: ", zmin, zmax
  if (zmin .ge. zmax) then
    write(*,*) "ERR: zmin must be smaller than zmax."
    stop 1
  end if
end if

! potential
if (analytic) then
  if ((run .eq. 0).and.(nstates.gt.1)) then
    write(*,*) "Analytic potention cannot be used for RT dynamics with nstates>1"
    stop 1
  end if
  write(*,*) "Potential: analytic"
  if (pot == "") then
    write(*,*) "Potential not provided! Use analytical form. x,y,z for corresponding rank "
    stop 1
  end if
else
  write(*,*) "Potential: provided in file"
end if
               
! number of states
if (nstates < 1) then
  write(*,*) "ERR: number of states must be 1 or more."
  stop 1
!else
!  if (run .eq. 0 .and. nstates > 1) then
!    write(*,*) "ERR: nstates > 1 available only for imag propagation."
!    stop 1
!  end if
else
    write(*,'(A,I2)') " nstates: ", nstates
end if

!projecting out rotated wf
if (project_rot .and. (run == 1) .and. (nstates>1)) then
  write(*,*) "Projecting out also 90 degrees rotated wavefunctions due to numeric instabilities"
end if

!norm check threshold
if (run.eq.0) then
  write(*,'(a, F14.8)') "Norm conservation threshold (norm_thresh) set to ", norm_thresh
  if (norm_thresh.le.0.0d0) then
    write(*,*) "ERROR: norm_thresh cannot be smaller than 0"
    stop 1
  end if
end if

!field
if (field_coupling) then
  select case(run)
  case(0)
    write(*,*) "Field: ON"
    !jj
    write(*,*) "WARNING! Field was not tested and is not fully implemented"
    if (field == '') then
      write(*,*) "ERR: No field function specified. Set 'field' in input.q."
      stop 1
    else
      write(*,*) "|E(t)| = ", field
    end if
  case(1)
    write(*,*) "ERR: Field cannot be used with imaginary time propagation."
    stop 1
  end select
else 
  if (run == 0) then
    write(*,*) "Field: OFF"
  end if
end if

!printing wavefunction
if (print_wf) then
  write(*,*) "Printing WF: ON"
else
  write(*,*) "Printing WF: OFF"
end if

end subroutine check

end module
