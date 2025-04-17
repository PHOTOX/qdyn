&general
    dynamics = 'rt'
    nstep = 4000 ! 4135 a.u. needed for rt dynamics
    dt = 1.0
    dtwrite = 100.0
    xngrid = 1024
    xmin = 0.5
    xmax = 8.0
    mass_x = 1728.256714
    nstates = 2
    print_wf = .false.
    rank = 1
/

&init_wf
    init_state = 1
    gen_init_wf = .false.
/

&rt
    analytic = .false.
    field_coupling = .true.
    field = '0.0715*sin(3.141592653589793*t/4135)**2*cos(0.17346*t)'    ! stationary to stationary state pulse
    field_off = 4135.0
    exact_factor = .true.
    ef_gauge = 'A0'
/
