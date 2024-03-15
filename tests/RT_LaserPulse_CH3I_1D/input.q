&general
    dynamics = 'rt'
    nstep = 10000
    dt = 1.0
    dtwrite = 2.0
    xngrid = 256
    xmin = -2.0
    xmax = 15.0
    mass_x = 6344.666496
    nstates = 2
    print_wf = .true.
    rank = 1
/

&init_wf
    x0 = 0.05283012
    xsigma = 0.24949412507853216
    px0 = 0.0
    init_state = 1
/

&rt
    analytic = .false.
    field_coupling = .true.
    field='0.05*exp(-(t-3000.5)**2/(2*785.468**2))*(cos(0.164489*t)-sin(0.164489*t)*(t-3000.5)/(785.468**2*0.164489))'
/
