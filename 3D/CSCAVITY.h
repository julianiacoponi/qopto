! input parameters for UCL- Antonio cavity
 PARAMETER(&
    NPERIOD=80000, & ! number of periods in prop. (of fastest oscillation: kappa^-1 or om_M)
    N_total=8, & ! number of equations so 8=>1 optical mode + 3D

    ! Vienna cavity bead parameters
    ! PARAMETER(&
    !   NPERIOD=160000, &
    !   N_total=8, &

    !   R0=71.5d-9, &
    !   rho=2198.d0, &
    !   EPSR=1.45d0**2, &
    !   Epsi0=8.854d-12, &

    R0=60.1d-9, & ! sphere radius /metres
    rho=1850.d0, & ! sphere density /kgm^-3 presumably?
    EPSR=1.98d0, & 
    Epsi0=8.854d-12, &

    c=3.d8, & ! speed of light /ms^-1
    hbar=1.05d-34, & ! Planck's constant /Js^-1
    kB=1.4d-23, & ! Boltzmann's constant /JK^-1
    TEMP=300.d0, & ! temperature /K
    Gravity=9.8d0, & ! acceleration due to Earth's gravity /ms^-2

    WK=5.9d6, & ! 2*pi/lambda=k
    WX=1.043d-6, &
    Wy=0.9035d-6, &
    waist=61.d-6, & ! /metres presumably?
    XL=1.223d-2, & ! /metres
    Finesse=3.d4, & ! ~ 1/kappa
    Press=3.e-3, & ! /millibars

    Pin1=0.44d0, & ! /Watts
    DET1=-176.d3, & ! detuning of trap beam i.e. Delta /KHz
    DelFSR=14.d9, & ! 1 FSR (Free Spectral Range) = 14 GHz
    theta0=0.25d0,& ! angle between tweezer polarization and cavity axis /pi

    ! Equilibrium positions
    ! node is X0 = 0.25*lambda ,Y0 = waist/sqrt(2)
    Y0=0.0d-6, &
    X0=0.25d0*1.064d-6, &
    Z0=0.d-6)
