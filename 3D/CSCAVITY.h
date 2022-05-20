! Vienna cavity bead parameters
! PARAMETER(&
!   N_period=160000, &
!   N_total=8, &

!   bead_radius=71.5d-9, &
!   bead_density=2198.d0, &
!   epsilon_R=1.45d0**2, &

! input parameters for UCL- Antonio cavity
 PARAMETER(&
    N_total=8, & ! number of equations so 8=>1 optical mode + 3D

    ! use at least 6 digits, from https://physics.nist.gov/cuu/Constants/
    hbar=1.054571817d-34, & ! Planck's reduced constant /Js^-1
    speed_of_light=299792498, & ! /ms^-1
    kB=1.380649d-23, & ! Boltzmann's constant /JK^-1
    vacuum_permittivity=8.854187d-12, & ! permittivity of free space /Farad.m^-1

    Free_Spectral_Range=14.d9, & ! /Hz
    bath_temperature=300.d0, & ! /K
    air_pressure=3.d-3, & ! /millibars
    air_speed=500.d0, & ! /ms^-1

    bead_diameter=120.2d-9, & ! /metres
    bead_density=1850.d0, & ! /kgm^-3
    bead_permittivity=1.98d0, & ! /dimensionless

    tweezer_input_power=0.485d0, & ! /Watts
    beam_waist_X=1.0679d-6, & ! /metres
    beam_waist_Y=0.9276d-6, & ! /metres
    theta_tweezer_over_pi=0.233d0, & ! angle between tweezer polarisation and cavity axis /π

    cavity_waist=61.d-6, & ! /metres
    cavity_length=1.223d-2, & ! /metres
    Finesse=3.d4, & ! ~number of times photons reflect off cavity mirrors before absorption /dimensionless

    half_kappa_exp_kHz=197.3, & ! cavity linewidth ±0.6kHz (most accurately measured part of the experiment) /kHz
    tweezer_wavelength=1.064d-6) ! /metres
