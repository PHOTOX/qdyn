&general
    dynamics = 'rt'
    nstep = 1000
    dt = 10.0
    dtwrite = 50.0
    xngrid = 1024
    xmin = -10.0
    xmax = 20.0
    mass_x = 2000.0
    nstates = 2
    print_wf = .false.
    rank = 1
/

&init_wf
    x0 = -1.0
    xsigma = 0.7071067811865475  ! 1/sqrt(2)
    px0 = 0.0 ! np.sqrt(2*mass*0.03)
    init_state = 2
/

&rt
    analytic = .false.
    field_coupling = .false.
/
