! input parameters
PARAMETER(&
    NPERIOD=80000, & ! number of periods in prop. (of fastest oscillation: kappa^-1 or om_M)
    N_total=4, & ! number of equations so 8=>1 optical mode +3D
    R0=71.5d-9, & ! sphere radius /metres

    ! NPERIOD=160000, &
    ! N_total=4, &
    ! R0=71.5d-9, &

    rho=2198.d0, & ! sphere density /kgm^-3 presumably?
    EPSR=2.1d0, &
    Epsi0=8.854d-12, &

    ! rho=2198.d0, &
    ! EPSR=1.45d0**2, &
    ! Epsi0=8.854d-12, &

    c=3.d8, & ! speed of light /ms^-1
    hbar=1.05d-34, & ! Planck's constant /Js^-1
    kB=1.4d-23, & ! Boltzmann's constant /JK^-1
    TEMP=300.d0, & ! temperature /K
    Gravity=9.8d0, & ! acceleration due to Earth's gravity /ms^-2

    WK=5.9d6, & ! 2*pi/lambda=k
    WX=1.02d-6, &
    WY=0.879d-6, &
    waist_radius=61.d-6, & ! /metres presumably?
    cavity_length=1.223d-2, & ! /metres
    Finesse=2.8d4, &

    air_pressure=1.e-6, & ! air pressure /millibars
    tweezer_input_power=0.3886d0, & ! /Watts
    detuning=-200.d3, & ! detuning of trap beam i.e. Delta? /KHz
    DelFSR=14.d9, & ! DelFSR = 1 FSR (Free Spectral Range) = 14 GHz
    theta0=0.5d0, & ! angle between tweezer polarization and cavity axis. Given as as a FRACTION of pi so pi/4  is theta0=0.25

    ! Equilibrium positions 
    ! X0=0.125*lambda, &
    ! Y0=waist_radius/sqrt(2), &
    X0=0.25d0*1.064d-6, &
    Y0=0.0d-6, &
    Z0=0.d-6)
