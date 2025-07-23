module mod_init
   use mod_vars
   use FFTW3
   use fparser, ONLY : initf, parsef, evalf, EvalErrType, EvalErrMsg
   use mod_utils
   use mod_exactfactor, ONLY : init_ef

   implicit none
   integer, private :: file_unit, jstate
   character(len = 50) :: file_name
   real(DP), dimension(:), allocatable :: point

CONTAINS

   subroutine init()
      implicit none
      integer :: istate, jstate

      write(*, *)
      write(*, *) "### Initialization ###"

      !-- Initialization of Wave Function
      ! initializing just for the 1. state in RT propagation
      ! initializing all states for IT propagation

      !-- GRID set-up
      dx = (xmax - xmin) / (xngrid - 1)
      dy = (ymax - ymin) / (yngrid - 1)
      dz = (zmax - zmin) / (zngrid - 1)

      !-- Allocating arrays
      ! grid
      if(rank == 1) allocate(x(xngrid), px(xngrid), point(1))
      if(rank == 2) allocate(x(xngrid), y(yngrid), px(xngrid), py(yngrid), point(2))
      if(rank == 3) allocate(x(xngrid), y(yngrid), z(zngrid), px(xngrid), py(yngrid), pz(zngrid), point(3))

      ! wf
      if(rank == 1) allocate(wfx(nstates, xngrid), wfp(xngrid))
      if(rank == 2) allocate(wf2x(nstates, xngrid, yngrid), wf2p(xngrid, yngrid))
      if(rank == 3) allocate(wf3x(nstates, xngrid, yngrid, zngrid), wf3p(xngrid, yngrid, zngrid))

      ! Hamiltonian
      if(run == 0) then
         if(rank == 1) then
            allocate(H1(nstates, nstates, xngrid), H1_ad(nstates, xngrid), U1(nstates, nstates, xngrid), &
                  invU1(nstates, nstates, xngrid))
            allocate(expH1(nstates, nstates, xngrid), wfx_ad(nstates, xngrid))
         end if
         if(rank == 2) allocate(expV2(xngrid, yngrid), v2(xngrid, yngrid))
         if(rank == 3) allocate(expV3(xngrid, yngrid, zngrid), v3(xngrid, yngrid, zngrid))

         if(field_coupling) then
            if(rank == 1) allocate(dipole_coupling(nstates, nstates, xngrid))
         end if

         if (autocorr) then
            allocate(autocorr_func(nstep + 1))
            if(rank == 1) allocate(wfx0(nstates, xngrid))
            if(rank == 2) allocate(wf2x0(nstates, xngrid, yngrid))
            if(rank == 3) allocate(wf3x0(nstates, xngrid, yngrid, zngrid))
            autocorr_func = 0.0d0
         end if

         ! allocating populations
         if (nstates>1) allocate(diab_pop(nstates), ad_pop(nstates))

      else if (run == 1) then
         if(rank == 1) allocate(expV1(xngrid), v1(xngrid))
         if(rank == 2) allocate(expV2(xngrid, yngrid), v2(xngrid, yngrid))
         if(rank == 3) allocate(expV3(xngrid, yngrid, zngrid), v3(xngrid, yngrid, zngrid))
      end if

      if(rank == 1) allocate(expT1(xngrid))
      if(rank == 2) allocate(expT2(xngrid, yngrid))
      if(rank == 3) allocate(expT3(xngrid, yngrid, zngrid))

      !-- Setting up grid poitns
      x(1) = xmin
      if(rank>1) y(1) = ymin
      if(rank>2) z(1) = zmin

      do i = 2, xngrid
         x(i) = x(i - 1) + dx
      end do

      do j = 2, yngrid
         y(j) = y(j - 1) + dy
      end do

      do k = 2, zngrid
         z(k) = z(k - 1) + dz
      end do

      !-- POTENTIAL energy init
      !- IT operators
      if (run==1) then
         ! creating potential
         if (analytic) then
            write(*, *)"Potential: ", pot
            call initf (1) !Initialization of parser

            select case(rank)
            case(1)
               call parsef (1, pot, (/'x'/)) !Bytcompiling function string
               do i = 1, xngrid
                  point = x(i)
                  v1(i) = evalf (1, point) !Evaluating potential for grid
                  if (EvalErrType > 0) then
                     WRITE(*, *)'*** Error evaluating potential: ', EvalErrMsg ()
                     stop 1
                  end if
               end do

            case(2)
               call parsef (1, pot, (/'x', 'y'/))
               do i = 1, xngrid
                  do j = 1, yngrid
                     point = (/x(i), y(j)/)
                     v2(i, j) = evalf (1, point)
                     if (EvalErrType > 0) then
                        WRITE(*, *)'*** Error evaluating potential: ', EvalErrMsg ()
                        stop 1
                     end if
                  end do
               end do

            case(3)
               call parsef (1, pot, (/'x', 'y', 'z'/))
               do i = 1, xngrid
                  do j = 1, yngrid
                     do k = 1, zngrid
                        point = (/x(i), y(j), z(k)/)
                        v3(i, j, k) = evalf (1, point)
                        if (EvalErrType > 0) then
                           WRITE(*, *)'*** Error evaluating potential: ', EvalErrMsg ()
                           stop 1
                        end if
                     end do
                  end do
               end do
            end select
            ! reading potential from file pot.dat
         else
            write(*, *) "Potential read from file:     pot.dat"
            open(667, file = 'pot.dat', status = 'OLD', action = 'READ', delim = 'APOSTROPHE', iostat = iost)
            select case(rank)
            case(1)
               read(667, *)v1
            case(2)
               read(667, *)v2
            case(3)
               read(667, *)v3
            end select
            close(667)
         end if

         ! IT: creating exponentail operator
         if(rank == 1) then
            do i = 1, xngrid
               expV1(i) = dcmplx(dexp(-v1(i) * dt / 2.0d0), 0)   !exp(-i V(x) tau/(2 h_bar))
            end do
         elseif(rank == 2) then
            do i = 1, xngrid
               do j = 1, yngrid
                  expV2(i, j) = dcmplx(dexp(-v2(i, j) * dt / 2.0d0), 0)
               end do
            end do
         elseif(rank == 3) then
            do i = 1, xngrid
               do j = 1, yngrid
                  do k = 1, zngrid
                     expV3(i, j, k) = dcmplx(dexp(-v3(i, j, k) * dt / 2.0d0), 0)
                  end do
               end do
            end do
         end if

         !- RT operators
      else if (run==0) then
         ! creating H matrix
         if (analytic) then
            write(*, *)"Potential: ", pot
            call initf (1)                                                     !Initialization of parser

            select case(rank)
            case(1)
               call parsef (1, pot, (/'x'/))                                       !Bytcompiling function string
               do i = 1, xngrid
                  point = x(i)
                  H1(1, 1, i) = evalf (1, point)                                       !Evaluating potential for grid
                  if (EvalErrType > 0) then
                     WRITE(*, *)'*** Error evaluating potential: ', EvalErrMsg ()
                     stop 1
                  end if
               end do

            case(2)
               call parsef (1, pot, (/'x', 'y'/))
               do i = 1, xngrid
                  do j = 1, yngrid
                     point = (/x(i), y(j)/)
                     v2(i, j) = evalf (1, point)
                     if (EvalErrType > 0) then
                        WRITE(*, *)'*** Error evaluating potential: ', EvalErrMsg ()
                        stop 1
                     end if
                  end do
               end do

            case(3)
               call parsef (1, pot, (/'x', 'y', 'z'/))
               do i = 1, xngrid
                  do j = 1, yngrid
                     do k = 1, zngrid
                        point = (/x(i), y(j), z(k)/)
                        v3(i, j, k) = evalf (1, point)
                        if (EvalErrType > 0) then
                           WRITE(*, *)'*** Error evaluating potential: ', EvalErrMsg ()
                           stop 1
                        end if
                     end do
                  end do
               end do
            end select
            ! reading H matrix from file
         else
            write(*, *) "Diabatic electronic Hamiltonian (H_el) read from files (H.i.j.dat)"
            select case(rank)
            case(1)
               call read_H1()
            end select
            !TODO: this is the old version for 2D and 3D that work only with 1 state now
            if (rank>1) then
               write(*, *) "Potential read from file:     pot.dat"
               open(667, file = 'pot.dat', status = 'OLD', action = 'READ', delim = 'APOSTROPHE', iostat = iost)
               select case(rank)
               case(2)
                  read(667, *)v2
               case(3)
                  read(667, *)v3
               end select
               close(667)
            end if
         end if

         ! adiabatic matrix
         if (nstates>=2) call adiab_trans_matrix()

         ! RT: creating operator
         ! if there is not field coupling (= no time-dependent part, we will built the propagators in advance)
         if (field_coupling) then
            write(*, *) "Time-dependent problem, exp(H_el) will be built every time step."
         else
            write(*, *) "Time-independent problem, exp(H_el) built."
            if(rank == 1) then

               call build_expH1()

            else if(rank == 2) then

               if (nstates==1) then
                  do i = 1, xngrid
                     do j = 1, yngrid
                        expV2(i, j) = dcmplx(dcos(-v2(i, j) * dt / 2.0d0), dsin(-v2(i, j) * dt / 2.0d0))   !exp(-i V(x) tau/(2 h_bar))
                     end do
                  end do
               else
                  write(*, *) "For rank=2, only 1 state propagation is available"
                  stop 1
               end if

            elseif(rank == 3) then
               if (nstates==1) then
                  do i = 1, xngrid
                     do j = 1, yngrid
                        do k = 1, zngrid
                           expV3(i, j, k) = dcmplx(dcos(-v3(i, j, k) * dt / 2.0d0), dsin(-v3(i, j, k) * dt / 2.0d0))
                        end do
                     end do
                  end do
               else
                  write(*, *) "For rank=3, only 1 state propagation is available"
                  stop 1
               end if
            end if
         end if
      end if

      !-- KINETIC energy init        exp[-iT/h_bar tau]
      select case(rank)
      case(1)
         do i = 1, xngrid
            if(i<=xngrid / 2) then
               px(i) = 2 * pi * (i - 1) / (xngrid * dx)
            else
               px(i) = 2 * pi * (i - 1 - xngrid) / (xngrid * dx)
            end if
            if(run == 0) expT1(i) = dcmplx(dcos(-px(i)**2 / (2 * mass_x) * dt), dsin(-px(i)**2 / (2 * mass_x) * dt))
            if(run == 1) expT1(i) = dcmplx(dexp(-px(i)**2 / (2 * mass_x) * dt), 0)
         end do
         !2D
      case(2)

         do i = 1, xngrid
            if(i <= xngrid / 2) then
               px(i) = 2 * pi * (i - 1) / (xngrid * dx)
            else
               px(i) = 2 * pi * (i - 1 - xngrid) / (xngrid * dx)
            end if
         end do

         do j = 1, yngrid
            if(j<=yngrid / 2) then
               py(j) = 2 * pi * (j - 1) / (yngrid * dy)
            else
               py(j) = 2 * pi * (j - 1 - yngrid) / (yngrid * dy)
            end if
         end do

         do i = 1, xngrid
            do j = 1, yngrid
               if(run == 0) expT2(i, j) = dcmplx(dcos(-(px(i)**2 / (2 * mass_x) + py(j)**2 / (2 * mass_y)) * dt), &
                     dsin(-(px(i)**2 / (2 * mass_x) + py(j)**2 / (2 * mass_y)) * dt))
               if(run == 1) expT2(i, j) = dcmplx(dexp(-(px(i)**2 / (2 * mass_x) + py(j)**2 / (2 * mass_y)) * dt), 0)
            end do
         end do

         !3D
      case(3)

         do i = 1, xngrid
            if(i<=xngrid / 2) then
               px(i) = 2 * pi * (i - 1) / (xngrid * dx)
            else
               px(i) = 2 * pi * (i - 1 - xngrid) / (xngrid * dx)
            end if
         end do

         do j = 1, yngrid
            if(j<=yngrid / 2) then
               py(j) = 2 * pi * (j - 1) / (yngrid * dy)
            else
               py(j) = 2 * pi * (j - 1 - yngrid) / (yngrid * dy)
            end if
         end do

         do k = 1, zngrid
            if(k<=zngrid / 2) then
               pz(k) = 2 * pi * (k - 1) / (zngrid * dz)
            else
               pz(k) = 2 * pi * (k - 1 - zngrid) / (zngrid * dz)
            end if
         end do

         do i = 1, xngrid
            do j = 1, yngrid
               do k = 1, zngrid
                  if(run == 0) expT3(i, j, k) = dcmplx(dcos(-(px(i)**2 / (2 * mass_x) + py(j)**2 / (2 * mass_y) + &
                        pz(k)**2 / (2 * mass_z)) * dt), &
                        dsin(-(px(i)**2 / (2 * mass_x) + py(j)**2 / (2 * mass_y) + pz(k)**2 / (2 * mass_z)) * dt))
                  if(run == 1) expT3(i, j, k) = dcmplx(dexp(-(px(i)**2 / (2 * mass_x) + py(j)**2 / (2 * mass_y) + &
                        pz(k)**2 / (2 * mass_z)) * dt), 0)
               end do
            end do
         end do
      end select

      !-- Field preparations
      if (field_coupling) then
         ! Initialize field function
         write(*, *) "Parsing field"
         call initf (2)                                                     !Initialization of parser
         call parsef (2, field, (/'t'/))                                       !Bytcompiling function string
         if (EvalErrType > 0) then
            WRITE(*, *)'*** Error parsing field: ', EvalErrMsg ()
         end if

         ! Read dipole couplings
         write(*, *) "Reading dipole couplings (file names expected to be dipole_coup.i.j.dat)"
         call read_dipole_coupling()
      end if


      !--Generating wavepacket
      call init_wavepacket()

      !--Calculating adiabatic wf
      if ((run==0).and.(nstates>1)) call wf_adiab_trans()

      !--Initializing EF
      if (exact_factor) call init_ef()

      !--Printing of WF
      if (print_wf) then
         write(*, *) "---------"

         if(rank == 1) then
            ! creating file name
            do jstate = 1, nstates

               file_unit = 200 + jstate
               write(file_name, *) jstate
               file_name = 'wf1d.' // trim(adjustl(file_name)) // '.dat'
               open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
               ! opening file unit
               write(file_unit, *) "#x   REAL   IMAG   NORM    POTENTIAL"
               close(file_unit)
               open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)

               call printwf_1d(jstate)
               write(*, *)"Diabatic wave function file:    " // file_name


               ! printing adiabatic
               if ((run==0).and.(nstates>1)) then
                  file_unit = 400 + jstate
                  write(file_name, *) jstate
                  file_name = 'wf1d_ad.' // trim(adjustl(file_name)) // '.dat'
                  open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
                  ! opening file unit
                  write(file_unit, *) "#x   REAL   IMAG   NORM    POTENTIAL"
                  close(file_unit)
                  open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)

                  call printwf_ad_1d(jstate)
                  write(*, *)"Adiabatic wave function file:   " // file_name

               end if

            end do

         elseif(rank == 2) then
            ! creating file name
            do jstate = 1, nstates

               file_unit = 200 + jstate
               write(file_name, *) jstate
               file_name = 'wf2d.' // trim(adjustl(file_name)) // '.dat'
               open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
               ! opening file unit
               write(file_unit, *) "#x  y   REAL   IMAG   NORM    POTENTIAL"
               close(file_unit)
               open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)

               call printwf_2d(jstate)
               write(*, *)"Diabatic wave function file:    " // file_name

            end do

         elseif(rank == 3) then
            ! creating file name
            do jstate = 1, nstates

               file_unit = 200 + jstate
               write(file_name, *) jstate
               file_name = 'wf3d.' // trim(adjustl(file_name)) // '.dat'
               open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
               ! opening file unit
               write(file_unit, *) "#x  y   z   REAL   IMAG   NORM    POTENTIAL"
               close(file_unit)
               open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)

               call printwf_3d(jstate)
               write(*, *)"Diabatic wave function file:    " // file_name
            end do

         end if
      end if

      !--Open file with energies
      select case(run)
      case(0)
         open(101, file = 'energies.dat', action = 'WRITE', iostat = iost)
         write(101, *) "#  time     total energy    potential       kinetic         energy diff     norm"
         close(101)
         open(101, file = 'energies.dat', status = 'old', position = 'append', action = 'WRITE', iostat = iost)

         !--Writing energies
         call update_norm()
         select case(rank)
         case(1)
            call update_total_energy_1d()
         case(2)
            call update_energy_2d(wf2x(1, :, :))
         case(3)
            call update_energy_3d(wf3x(1, :, :, :))
         end select
         call printen()
      case(1)
         do jstate = 1, nstates
            file_unit = 300 + jstate
            write(file_name, *) jstate
            file_name = 'energies.' // trim(adjustl(file_name)) // '.dat'
            open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
            write(file_unit, *) "#  time     total energy    potential       kinetic         energy diff"
            close(file_unit)
            open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)

            !--Writing energies
            select case(rank)
            case(1)
               call update_energy_1d(wfx(1, :))
            case(2)
               call update_energy_2d(wf2x(1, :, :))
            case(3)
               call update_energy_3d(wf3x(1, :, :, :))
            end select

            call printen_state(jstate)
         end do
      end select

      !--Open file with field
      if (field_coupling) then
         file_unit = 102
         file_name = 'field.dat'
         open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
         write(file_unit, *) "#  time     field amplitude"
         close(file_unit)
         open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)

         write(*, '(A,A)') " Electric field file:            ", file_name

         call print_field()

      end if

      ! open autocorr file and save initial wave function and initial value of S
      if (autocorr) then
         file_unit = 105
         file_name = 'autocorr_function.dat'
         open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
         write(file_unit, *) "#  time     Re[S(t)]      Im[S(t)]"
         close(file_unit)
         open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)

         write(*, '(A,A)') " Autocorrelation fucntion file:  ", file_name

         ! saving initial wave function
         select case(rank)
         case(1)
            wfx0 = wfx(:, :)
         case(2)
            wf2x0 = wf2x(:, :, :)
         case(3)
            wf3x0 = wf3x(:, :, :, :)
         end select

         ! calulating autocorrelation function and printing
         call calc_autocorr(0)
         call print_autocorr(0)

      end if

      !--Open file with diabatic populaitons
      if ((run==0).and.(nstates>1)) then

         file_unit = 103
         file_name = 'pop_diab.dat'
         open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
         write(file_unit, *) "#  time     diabatic populations (1, 2, 3, ...)   norm"
         close(file_unit)
         open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)
         write(*, '(A,A)') " Diabatic populations file:      ", trim(file_name)

         file_unit = 104
         file_name = 'pop_ad.dat'
         open(file_unit, file = file_name, action = 'WRITE', iostat = iost)
         write(file_unit, *) "#  time     adiabatic populations (1, 2, 3, ...)  norm"
         close(file_unit)
         open(file_unit, file = file_name, status = 'old', position = 'append', action = 'WRITE', iostat = iost)
         write(*, '(A,A)') " Adiabatic populations file:     ", trim(file_name)

         call print_pop()

      end if

      !closing input.q
      close(100)

      write(*, *)

   end subroutine init

   subroutine init_wavepacket()

      real(DP) :: x0, y0, z0, xsigma, ysigma, zsigma, px0, py0, pz0 ! read in input
      real(DP) :: prefactor, gauss, momenta, norm ! for generating gaussian wf
      integer :: init_state = 1 ! read from input, inital state for the dynamics
      real(DP) :: weights(nstates) ! weights for the initial wave function, read from input (input example 'weights = 0.0, 0.6, 0.4'
      ! where the number of numbers cannot exceed number of states
      complex(DP), allocatable :: aux_wf(:, :, :) ! auxiliary variable for reading wf
      logical :: gen_init_wf = .true. ! flag to generate initial wave function as Gaussian, read by namelist
      logical :: use_weights = .false. ! flag to use weights for initial wave function, determined from elements of weights
      character(len = 1000) :: line

      namelist /init_wf/x0, y0, z0, xsigma, ysigma, zsigma, px0, py0, pz0, init_state, gen_init_wf, weights

      ! Default values for generated gaussian wavepacket
      x0 = (xmax - xmin) / 2.0d0 + xmin + (xmax - xmin) / 10                      ! shifted a little bit so that the wave packet is not symmetric
      y0 = (ymax - ymin) / 2.0d0 + ymin + (ymax - ymin) / 10                      ! shifted a little bit so that the wave packet is not symmetric
      z0 = (zmax - zmin) / 2.0d0 + zmin + (zmax - zmin) / 10                      ! shifted a little bit so that the wave packet is not symmetric
      xsigma = (xmax - x0) / 20.0d0                                       ! 5sigma - 96% of gaussian is on the grid
      ysigma = (ymax - y0) / 20.0d0                                       ! 5sigma - 96% of gaussian is on the grid
      zsigma = (zmax - z0) / 20.0d0                                       ! 5sigma - 96% of gaussian is on the grid
      px0 = 0.0d0
      py0 = 0.0d0
      pz0 = 0.0d0

      ! weights for initial wavefunction
      ! this option takes the initial wf and puts it exactly as it is on all electronic states with provided weights without phase
      ! initially, we set the option to 0 for all elements, the code then checks it and if there are non-zero elements, it generates
      ! the wf by spreading it over all electronic states
      ! if not provided, the code generates the wf using init_state keyword putting the wf completely on the given state
      weights = 0.0d0

      ! allocating aux_wf
      if (rank==1)  then
         allocate(aux_wf(xngrid, 1, 1))
      elseif (rank==2) then
         allocate(aux_wf(xngrid, yngrid, 1))
      elseif (rank==3) then
         allocate(aux_wf(xngrid, yngrid, zngrid))
      end if


      ! intial printing
      write(*, *) "---------"
      write(*, *) "Generating initial wave packet."

      ! reading &init_wf section in the input.q
      rewind(100)
      read(100, init_wf, iostat = iost)
      if (iost.ne.0) then
         write(*, *)'ERROR: &init_wf section either not provided or incorrect.'
         backspace(100)
         read(100, fmt = '(A)') line
         write(*, '(A)') ' Invalid line in namelist: ' // trim(line)
         stop 1
      end if
      close(100)

      ! init_state keyword check
      if ((init_state<1).or.(init_state>nstates)) then
         write(*, *) "Initial state is greater than nstates or less than 1"
         stop 1
      end if

      ! weights checked
      if (run==1) weights = 0.0 ! for IT dynamics we always initialize in state 1 as we then copy to the other state

      ! check weights and decide whether to use them
      ! check if there is nonzero element in weights
      do i = 1, nstates
         ! check nonzaro values
         if (weights(i) /= 0.0d0) then
            use_weights = .true. ! if there is at least one non-zero element, we use weights
         end if
      end do

      if ((use_weights).and.(run==0)) then
         write(*, *) "Initializing wave function in all states according to their weights."
         do i = 1, nstates
            !check if it is a positive value
            if (weights(i) < 0.0d0) then
               write(*, *) "ERROR: Initial wave function weights must be non-negative."
               write(*, *) "       state: ", i, ",  weight: ", weights(i)
               stop 1
            end if
         end do

         ! renormalization
         write(*, *) "Renormalizing state weights."
         norm = sum(weights**2)
         weights = weights / sqrt(norm) ! renormalizing weights
      end if

      ! write what will be used for initialization
      if (run==0) then
         if (use_weights) then
            write(*, '(A17,100(F10.6))') " State weights:      ", weights
         else
            write(*, '(A17,I10)') " Initial state:      ", init_state
         end if
      end if

      ! now prepare the weights according to the calculation type since it's used to declare initial wave function in all cases
      if (run==1) then ! for IT dynamics, we initiate at all states the same wave function
         weights = 1.0d0
      else if ((run==0).and.(.not.use_weights)) then ! for RT dynamics with weights not set, use the initial state
         weights = 0.0d0 ! just to make sure its zeros
         weights(init_state) = 1.0d0 ! giving the only weight to the initial state
      end if ! else we use the weights provided by user

      if (gen_init_wf) then
         !         write(*, *) "---------"
         write(*, *) "* Gaussian wave packet parameters:"
         write(*, '(A17,F10.5)') " x0:              ", x0
         if (rank>=2) write(*, '(A17,F10.5)') " y0:               ", y0
         if (rank>=3) write(*, '(A17,F10.5)') " z0:               ", z0
         write(*, '(A17,F10.5)') " xsigma:                ", xsigma
         if (rank>=2) write(*, '(A17,F10.5)') " ysigma:                ", ysigma
         if (rank>=3) write(*, '(A17,F10.5)') " zsigma:                ", zsigma
         write(*, '(A17,F10.5)') " px0:                ", px0
         if (rank>=2) write(*, '(A17,F10.5)') " py0:                ", py0
         if (rank>=3) write(*, '(A17,F10.5)') " pz0:                ", pz0

         ! normalization prefactor of the gaussian, same value for the whole grid
         select case (rank)
         case (1)
            prefactor = xsigma**(-0.5) * pi**(-0.25)
         case (2)
            prefactor = (xsigma * ysigma)**(-0.5) * pi**(-0.5)
         case (3)
            prefactor = (xsigma * ysigma * zsigma)**(-0.5) * pi**(-0.75)
         end select

         do i = 1, xngrid
            select case (rank)
            case (1)
               ! Gaussian part of the wave packet depending only on positions
               gauss = prefactor * dexp(-(x(i) - x0)**2 / (2 * xsigma**2))
               ! momentum calculation p0*(x-x0)
               momenta = px0 * (x(i) - x0)
               ! whole complex wf
               aux_wf(i, 1, 1) = dcmplx(gauss * dcos(momenta), gauss * dsin(momenta))

               do jstate = 1, nstates
                  wfx(jstate, i) = weights(jstate) * aux_wf(i, 1, 1)
               end do

            case (2)
               do j = 1, yngrid
                  ! Gaussian part of the wave packet depending only on positions
                  gauss = prefactor * dexp(-(x(i) - x0)**2 / (2 * xsigma**2) - (y(j) - y0)**2 / (2 * ysigma**2))
                  ! momentum calculation p0*(x-x0)
                  momenta = px0 * (x(i) - x0) + py0 * (y(j) - y0)
                  ! whole complex wf
                  aux_wf(i, j, 1) = dcmplx(gauss * dcos(momenta), gauss * dsin(momenta))

                  do jstate = 1, nstates
                     wf2x(jstate, i, j) = weights(jstate)*aux_wf(i, j, 1)
                  end do
               end do

            case (3)
               do j = 1, yngrid
                  do k = 1, zngrid
                     ! Gaussian part of the wave packet depending only on positions
                     gauss = prefactor * dexp(-(x(i) - x0)**2 / (2 * xsigma**2) - (y(j) - y0)**2 / (2 * ysigma**2) - &
                           (z(k) - z0)**2 / (2 * zsigma**2))
                     ! momentum calculation p0*(x-x0)
                     momenta = px0 * (x(i) - x0) + py0 * (y(j) - y0) + pz0 * (z(k) - z0)
                     ! whole complex wf
                     aux_wf(i, j, k) = dcmplx(gauss * dcos(momenta), gauss * dsin(momenta))

                     do jstate = 1, nstates
                        wf3x(jstate, i, j, k) = weights(jstate) * aux_wf(i, j, k)
                     end do
                  end do
               end do
            end select
         end do

      else
         write(*, *) "WF will be read from wf.chk file."
         open(666, file = 'wf.chk', status = 'OLD', action = 'READ', delim = 'APOSTROPHE', iostat = iost)
         if (iost.ne.0) then
            write(*, *)'ERROR: wf.chk file must be provided'
            write(*, *) iost
            stop 1
         end if

         !Procedure for loading WF from file
         read(666, *) !two empty lines
         read(666, *)
         if(rank == 1) read(666, *) aux_wf(:, 1, 1)
         if(rank == 2) read(666, *) aux_wf(:, :, 1)
         if(rank == 3) read(666, *) aux_wf(:, :, :)
         close(666)

         do jstate = 1, nstates
            if(rank == 1) wfx(jstate, :) = weights(jstate) * aux_wf(:, 1, 1)
            if(rank == 2) wf2x(jstate, :, :) = weights(jstate) * aux_wf(:, :, 1)
            if(rank == 3) wf3x(jstate, :, :, :) = weights(jstate) * aux_wf(:, :, :)
         end do
      end if

      ! Normalize IT dynamics wf - it shouldn't be necessary as the Gaussians are normalized but renormalizatino might be
      ! necessary if the wave function is excaping the grid. While this is a problem for RT dynamics, it is not for IT dynamics
      ! since we collapse the wave function to the eigenstates and renormalize at every time step.
      if(run==1) then
         do jstate = 1, nstates
            if(rank == 1) call normalize_1d(wfx(jstate, :))
            if(rank == 2) call normalize_2d(wf2x(jstate, :, :))
            if(rank == 3) call normalize_3d(wf3x(jstate, :, :, :))
         end do
      else if (run==0) then ! update norm for RT dynamics, if the norm exceeds threshold, it will stop qdyn
         call update_norm() ! this is for RT dynamics, we normalize the initial wave function
      end if

   end subroutine

   subroutine read_dipole_coupling()

      integer :: istate, jstate
      logical :: file_exists

      do istate = 1, nstates
         do jstate = 1, nstates

            write(file_name, '(I1,A,I1)') istate, ".", jstate
            file_name = 'dipole_coup.' // trim(adjustl(file_name)) // '.dat'
            inquire(file = file_name, exist = file_exists)

            if (file_exists) then
               file_unit = 400
               open(file_unit, file = file_name, status = 'OLD', action = 'READ', delim = 'APOSTROPHE', iostat = iost)
               read(file_unit, *) dipole_coupling(istate, jstate, :)
               close(file_unit)
               write(*, '(a,i1,a,i1,a)') " * <psi_", istate, "|mu|psi_", jstate, "> = read from file"
            else
               write(*, '(a,i1,a,i1,a)') " * <psi_", istate, "|mu|psi_", jstate, "> = 0"
               dipole_coupling(istate, jstate, :) = 0.0d0
            end if
         end do
      end do

   end subroutine

   subroutine read_H1()

      integer :: istate, jstate
      logical :: file_exists

      do istate = 1, nstates
         do jstate = 1, nstates

            write(file_name, '(I1,A,I1)') istate, ".", jstate
            file_name = 'H.' // trim(adjustl(file_name)) // '.dat'
            inquire(file = file_name, exist = file_exists)

            if (file_exists) then
               file_unit = 400
               open(file_unit, file = file_name, status = 'OLD', action = 'READ', delim = 'APOSTROPHE', iostat = iost)
               read(file_unit, *) H1(istate, jstate, :)
               close(file_unit)

               write(*, '(a,i1,a,i1,a)') " * <psi_", istate, "|H_el|psi_", jstate, "> = read from file"
            else
               if (jstate==istate) then
                  write(*, '(a,i1,a,i1,a)') "ERROR: H.", istate, ".", jstate, ".dat not found!"
                  stop 1
               end if
               write(*, '(a,i1,a,i1,a)') " * <psi_", istate, "|H_el|psi_", jstate, "> = 0"
               H1(istate, jstate, :) = 0.0d0
            end if
         end do
      end do

   end subroutine

end module

