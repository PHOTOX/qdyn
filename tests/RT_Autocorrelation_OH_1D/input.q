&general
    dynamics = 'rt'
    nstep = 30000
    dt = 0.5
    dtwrite = 5000.0
    !dtwrite = 1000.0
    xngrid = 512
    xmin = 1.0
    xmax = 5.0
    mass_x = 1728.256714
    nstates = 2
    print_wf = .false.
    rank = 1
/

&init_wf
    init_state = 1
    gen_init_wf = .true.
    x0 = 2.5
    xsigma = 0.15
/

&rt
    analytic = .false.
    autocorr= .true.
/
