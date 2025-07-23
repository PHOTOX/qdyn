&general
    dynamics = 'rt'
    nstep = 350
    dt = 10.0
    dtwrite = 50.0
    xngrid = 2048
    xmin = -18.0
    xmax = 32.0
    mass_x = 2000.0
    nstates = 3
    print_wf = .false.
    rank = 1
/

&init_wf
    x0 = -15.0
    xsigma = 0.7071067811865475
    px0 = 20.0
    init_state = 3
/

&rt
    analytic = .false.
    field_coupling = .false.
/
