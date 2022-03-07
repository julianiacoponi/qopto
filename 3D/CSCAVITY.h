! Vienna cavity bead parameters
! PARAMETER(&
!   N_period=160000, &
!   N_total=8, &

!   bead_radius=71.5d-9, &
!   bead_density=2198.d0, &
!   epsilon_R=1.45d0**2, &
!   epsilon_0=8.854d-12, &

! input parameters for UCL- Antonio cavity
 PARAMETER(&
    N_period=80000, & ! number of periods in prop. (of fastest oscillation: kappa^-1 or om_M)
    N_total=8, & ! number of equations so 8=>1 optical mode + 3D

    bead_radius=60.1d-9, & ! /metres
    bead_density=1850.d0, & ! /kgm^-3 presumably?

    vacuum_permittivity=8.854d-12, & ! permittivity of free space /Farad.m^-1
    relative_permittivity=1.98d0, & ! /dimensionless

    speed_of_light=3.d8, & ! /ms^-1
    hbar=1.05d-34, & ! Planck's constant /Js^-1
    kB=1.4d-23, & ! Boltzmann's constant /JK^-1
    bath_temperature=300.d0, & ! /K
    Gravity=9.8d0, & ! acceleration due to Earth's gravity /ms^-2

    linewidth_k=5.9d6, & ! k = 2π/lambda
    beam_waist_X=1.043d-6, & ! /metres
    beam_waist_Y=0.9035d-6, & ! /metres

    cavity_waist=61.d-6, & ! /metres
    cavity_length=1.223d-2, & ! /metres
    Finesse=3.d4, & ! ~ 1/kappa /dimensionless
    air_pressure=3.e-3, & ! /millibars

    tweezer_input_power=0.44d0, & ! /Watts
    DelFSR=14.d9, & ! 1 FSR (Free Spectral Range) = 14 GHz
    theta_0=0.25d0, & ! angle between tweezer polarization and cavity axis /pi

    detuning_Hz=-176.d3, & ! detuning of trap beam i.e. laser frequency minus cavity frequency a.k.a. Delta /Hz
    detuning_Hz_Antonio_0_91 = -178.36d3, & ! calculated detuning for Antonio's data 0.91
    detuning_Hz_Antonio_1_825 = -357.7d3, & ! calculated detuning for Antonio's data 1.825

    ! Equilibrium positions
    ! node is X0 = 0.25*lambda ,Y0 = cavity_waist/sqrt(2)
    X0=0.25d0*1.064d-6, &
    Y0=0.0d-6, &
    Z0=0.d-6)
