&general
    dynamics = 'rt'
    nstep = 4000
    dt = 10.0
    dtwrite = 50.0
    xngrid = 512
    xmin = -10.0
    xmax = 25.0
    mass_x = 2000.0
    wf = 0
    nstates = 2
    print_wf = .false.
    rank = 1
/

&init_wf
    x0 = -5.0
    xsigma = 0.7071067811865475  ! 1/sqrt(2)
    px0 = 10.954451150103322 ! np.sqrt(2*mass*0.03)
    init_state = 1
/

&rt
    analytic = .false.
    field_coupling = .false.
/
