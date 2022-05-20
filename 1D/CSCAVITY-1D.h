! input parameters
PARAMETER(&
    N_period=80000, & ! number of periods in prop. (of fastest oscillation: kappa^-1 or om_M)
    N_total=4, & ! number of equations so 8=>1 optical mode +3D

    ! NPERIOD=160000, &
    ! N_total=4, &
    ! bead_radius=71.5d-9, &

    bead_radius=71.5d-9, & ! /metres
    bead_density=2198.d0, & ! /kgm^-3 presumably?
    epsilon_R=2.1d0, & ! relative permittivity /dimensionless
    epsilon_0=8.854d-12, &  ! permittivity of free space / Farad.m^-1

    ! rho=2198.d0, &
    ! epsilon_R=1.45d0**2, &
    ! epsilon_0=8.854d-12, &

    c=3.d8, & ! speed of light /ms^-1
    hbar=1.05d-34, & ! Planck's constant /Js^-1
    kB=1.4d-23, & ! Boltzmann's constant /JK^-1
    bath_temperature=300.d0, & ! /K
    Gravity=9.8d0, & ! acceleration due to Earth's gravity /ms^-2

    WK=5.9d6, & ! 2*pi/lambda=k
    WX=1.02d-6, &
    WY=0.879d-6, &
    waist_radius=61.d-6, & ! /metres
    cavity_length=1.223d-2, & ! /metres
    Finesse=2.8d4, &
    air_pressure=1.e-6, & ! air pressure /millibars

    tweezer_input_power=0.3886d0, & ! /Watts
    detuning_Hz=-200.d3, & ! detuning of trap beam i.e. laser frequency minus cavity frequency a.k.a. Delta /Hz
    DelFSR=14.d9, & ! DelFSR = 1 FSR (Free Spectral Range) = 14 GHz
    theta_0=0.5d0, & ! angle between tweezer polarization and cavity axis /pi

    ! Equilibrium positions 
    ! X0=0.125*lambda, &
    ! Y0=waist_radius/sqrt(2), &
    X0=0.25d0*1.064d-6, &
    Y0=0.0d-6, &
    Z0=0.d-6)
