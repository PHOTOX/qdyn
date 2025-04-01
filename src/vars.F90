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
  real(DP)              :: field_on=0.0d0, field_off ! turning on and off field
  character(len=9000)   :: field=''
  !-exact factorization
  !TODO: finish exact factorization derivatives
  logical               :: exact_factor=.false.
  character(len=2)      :: ef_gauge='S0' ! possible gauges: 'S0' (zero nuclear phase), 'A0' (zero TDVP)
  complex(DP), dimension(:,:,:), allocatable :: wfx_history ! used to calculate time derivatives of the wave function (hist. index, state, x)
  integer               :: efhistory=5 ! how many steps of wf we use for calcuations of wf derivative
  !-auxiliary variables
  integer               :: iost, i, j, k
  integer ( kind = 8 )  :: plan_forward, plan_backward
  character(len=9000)   :: pot=''
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
  namelist /rt/ pot, analytic, field_coupling, field, field_on, field_off, norm_thresh, &
     exact_factor, ef_gauge
  !TODO: rt analytic will not be use probably

CONTAINS

subroutine read_input()
  implicit none
  character(len=1000) :: line ! this line is used when namelist reading fails and stores the last input line

  write(*,*)
  write(*,*) "### Reading input ###"

!-- Reading input file
  open(100,file='input.q', status='OLD', action='READ',delim='APOSTROPHE', iostat=iost)
  if (iost/=0) then
    write(*,*)'ERROR: input.q file must be provided'
    write(*,*) iost
    stop 1
  end if

  read(100, general, iostat=iost)
  if (iost/=0) then
    write(*,*)'ERROR: &general section missing or problematic'
    backspace(100)
    read(100,fmt='(A)') line
    write(*,'(A)') ' Invalid line in namelist: '//trim(line)
    stop 1
  end if

  rewind(100)
  select case(dynamics)
  case('rt')
    ! setting field_off to max time before &rt is read
    field_off=nstep*dt
    ! reading &rt
    read(100, rt, iostat=iost)
    if (iost/=0) then
      write(*,*)'ERROR: &rt section missing or problematic'
      backspace(100)
      read(100,fmt='(A)') line
      write(*,'(A)') ' Invalid line in namelist: '//trim(line)
      write(*,*) iost
      stop 1
    end if
  case('it')
    read(100, it, iostat=iost)
    if (iost/=0) then
      write(*,*)'ERROR: &it section missing or problematic'
      backspace(100)
      read(100,fmt='(A)') line
      write(*,'(A)') ' Invalid line in namelist: '//trim(line)
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
    write(*,*) "* REAL TIME PROPAGATION *"
    run = 0
  case('it')
    write(*,*) "* IMAGINARY TIME PROPAGATION *"
    run = 1
  case default
    write(*,*) "ERR: Unrecongnized 'dynamics' option. Choose either 'rt' or 'it'. Exiting"
    stop 1
end select

! writing time step and number of steps
write(*,'(A25,F10.5,A)') " Time step:                ", dt, " a.u."
write(*,'(A25,I10)') " Number of steps:        ", nstep
write(*,'(A25,F10.2,A)') " Total time:              ", dt*nstep, " a.u."

! ngrid is power of 2
write(*,'(A25,I10)') " Grid size in x:            ",xngrid
if ((xngrid <= 0) .or. (IAND(xngrid, xngrid-1) /= 0))  then
  write(*,*) "ERR: Grid size must be power of two."
  stop 1
end if 

if (rank >= 2) then
  write(*,'(A25,I10)') " Grid size in y:           ",yngrid
  if ((yngrid <= 0) .or. (IAND(yngrid, yngrid-1) /= 0))  then
    write(*,*) "ERR: Grid size must be power of two."
    stop 1
  end if 
end if

if (rank >= 3) then
  write(*,'(A25,I10)') " Grid size in z:           ",zngrid
  if ((zngrid <= 0) .or. (IAND(zngrid, zngrid-1) /= 0))  then
    write(*,*) "ERR: Grid size must be power of two."
    stop 1
  end if 
end if

! dimensionality
if ((rank < 1) .or. (rank > 3)) then
  write(*,*) "ERR: Dimensionality must be 1,2 or 3."
  stop 1
else
  write(*,'(A25,I10)') " Number of dimensions:       ", rank
end if

! number of steps
if (nstep < 1) then
  write(*,*) "ERR: Number of steps must be bigger than 1."
  stop 1
else
write(*,'(A25,I10)') " Number of steps:          ",nstep
end if

! masses
if (mass_x <= 0.0) then
  write(*,*) "ERR: mass_x was either not set or set negative."
  stop 1
end if

if ((rank >= 2) .and. (mass_y <= 0.0)) then
  write(*,*) "ERR: mass_y was either not set or set negative."
  stop 1
end if

if ((rank >= 3) .and. (mass_z <= 0.0)) then
  write(*,*) "ERR: mass_z was either not set or set negative."
  stop 1
end if

! params of grid
write(*,'(A25,F10.4,A)') " xmin:                   ", xmin, " a.u."
write(*,'(A25,F10.4,A)') " xmax:                   ", xmax, " a.u."
if (xmin >= xmax) then
  write(*,*) "ERR: xmin must be smaller than xmax."
  stop 1
end if

if (rank >= 2) then
write(*,'(A25,F10.4,A)') " ymin:                   ", ymin, " a.u."
write(*,'(A25,F10.4,A)') " ymax:                   ", ymax, " a.u."
  if (ymin >= ymax) then
    write(*,*) "ERR: ymin must be smaller than ymax."
    stop 1
  end if
end if

if (rank >= 3) then
write(*,'(A25,F10.4,A)') " zmin:                   ", zmin, " a.u."
write(*,'(A25,F10.4,A)') " zmax:                   ", zmax, " a.u."
  if (zmin >= zmax) then
    write(*,*) "ERR: zmin must be smaller than zmax."
    stop 1
  end if
end if

! potential
if (analytic) then
  if ((run == 0).and.(nstates>1)) then
    write(*,*) "ERROR: Analytic potention cannot be used for RT dynamics with nstates>1"
    stop 1
  end if
  write(*,*) "Potential:                analytic"
  if (pot == "") then
    write(*,*) "Potential not provided! Use analytical form. x,y,z for corresponding rank "
    stop 1
  end if
else
  write(*,*) "Potential:                    file"
end if
               
! number of states
if (nstates < 1) then
  write(*,*) "ERR: number of states must be 1 or more."
  stop 1
else
    write(*,'(A25,I10)') " Number of states:          ", nstates
end if

!exact factorization
if (exact_factor) then
  if ((run /= 0) .or. (rank /= 1) .or. (nstates /= 2)) then
    write(*,*) "ERROR: Exact factorization available only for 1D two-state RT dynamics."
    stop 1
  end if
  write(*,*) "Exact factorization:            ON"
   if (ef_gauge == 'S0') then
      write(*,*) "Gauge:                          S0 (zero nuclear phase)"
   else if (ef_gauge == 'A0') then
      write(*,*) "Gauge:                          A0 (zero TDVP)"
   end if
end if

!printing wavefunction
if (print_wf) then
  write(*,*) "Printing WF:                    ON"
else
  write(*,*) "Printing WF:                   OFF"
end if

!field
if (field_coupling) then
  write(*,*) "Electric field:                 ON"
  write(*,'(A25,F10.2,A)') " Field on at              ", field_on, " a.t.u."
  write(*,'(A25,F10.2,A)') " Field off at             ", field_off, " a.t.u."
  if (field == '') then
    write(*,*) "ERR: No field function specified. Set 'field' in input.q."
    stop 1
  else
    write(*,*) "|E(t)| = ", trim(field)
  end if
else
  if (run == 0) then
    write(*,*) "Field: OFF"
  end if
end if

!projecting out rotated wf
if (project_rot .and. (run == 1) .and. (nstates>1)) then
  write(*,*) "Projecting out also 90 degrees rotated wavefunctions due to numeric instabilities"
end if

!norm check threshold
if (run==0) then
  write(*,'(a, F14.8)') " Norm conservation threshold (norm_thresh) set to ", norm_thresh
  if (norm_thresh<=0.0d0) then
    write(*,*) "ERROR: norm_thresh cannot be smaller than 0"
    stop 1
  end if
end if




end subroutine check

end module mod_vars