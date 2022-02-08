! input parameters for UCL- Antonio cavity
! NPERIOD = number of periods in prop. (of fastest oscillation: kappa^-1 or om_M)
! N_total = number of equations so 8 => 1 optical mode + 3D
! R0 = sphere radius
! XL = cavity length
! Pin1 = input power in Watts tweezer beam

! Press = air pressure in millibars
! XL = cavity length
! rho = sphere density
! WK = 2pi/lambda = k
! DelFSR = 1 FSR= 14 GHz
! DET1 = detuning in KHz trap beam
! Theta0 = angle between tweezer polarization and cavity axis. Given as as
! FRACTION of pi so pi/4  is theta0 = 0.25
! waist is  waist radius

 PARAMETER(&
 NPERIOD=80000, & ! number of periods in prop. (of fastest oscillation: kappa^-1 or om_M)
 N_total=4, & ! number of equations so 8=>1 optical mode +3D

    NPERIOD=80000, &
    N_total=8, &

    ! Vienna cavity bead parameters
    ! PARAMETER(&
    !   NPERIOD=160000, &
    !   N_total=8, &

    !   R0=71.5d-9, &
    !   rho=2198.d0, &
    !   EPSR=1.45d0**2, &
    !   Epsi0=8.854d-12, &

    R0=60.1d-9, &
    rho=1850.d0, &
    EPSR=1.98d0, &
    Epsi0=8.854d-12, &

    c=3.d8, &
    hbar=1.05d-34, &
    BOLTZ=1.4d-23, &
    TEMP=300.d0, &
    Gravity=9.8d0, &

    WK=5.9d6, &
    waist=61.d-6, &
    WX=1.043d-6, &
    Wy=0.9035d-6, &
    XL=1.223d-2, &
    Finesse=3.d4, &
    Press=3.e-3, &

    Pin1=0.44d0, &
    DET1=-176.d3, &
    DelFSR=14.d9, &
    theta0=0.25d0,&
    ! Equilibrium positions
    ! node is X0 = 0.25*lambda ,Y0 = waist/sqrt(2)
    Y0=0.0d-6, &
    X0=0.25d0*1.064d-6, &
    Z0=0.d-6)
