! Code to calculate mechanical frequencies for different equilibrium positions in the x-axis

IMPLICIT NONE
integer::NPERIOD, NTOT
double precision::R0, rho, EPSR, EPSI0
double precision::C, hbar, BOLTZ, TEMP, Gravity
double precision::WK, waist_radius, WX, WY, XL, Finesse, air_pressure
double precision::tweezer_input_power, detuning, DelFSR, Theta0
double precision::X0, Y0, Z0

! parameter file with input values
INCLUDE 'CSCAVITY.h'

integer::ii
double precision::x0_sweep
double precision::omega_x, omega_y, omega_z

! calculated  parameters
double precision::detuning, Wkx0, thet, polarisability
double precision::bead_mass, Rayleigh_range
double precision::E_tweezer, E_cavity
double precision::kappa_in, kappa_nano, Gamma_M

double precision::pi, pi2


pi = dacos(-1.d0)
pi2 = 2.d0 * pi
detuning = detuning * pi2

DO ii = 1, 100
    ! increase x0 in increments of 0.005*lambda
    x0_sweep = ii * 0.005d0 * 1.064d-6
    Wkx0 = Wk * x0_sweep
    thet = theta0 * pi

    ! note BEAD routine uses x0 but not theta
    ! NOTE: does changing the above need to be done to do x0_sweep properly?
    CALL CALCULATE_BEAD_PARAMETERS(&
        Wkx0, &
        polarisability, bead_mass, Rayleigh_range, &
        E_tweezer, E_cavity, &
        kappa_in, kappa_nano, Gamma_M)

    CALL CALCULATE_MECH_FREQUENCIES(&
        thet, detuning, Wkx0, &
        polarisability, bead_mass, Rayleigh_range, &
        E_tweezer, E_cavity, &
        kappa_in, kappa_nano, Gamma_M, &
        omega_x, omega_y, omega_z, &
        x0_sweep)

    WRITE(*, *) x0_sweep, omega_x/2/pi, omega_y/2/pi, omega_z/2/pi
END DO

STOP
END

SUBROUTINE CALCULATE_BEAD_PARAMETERS(&
    thet, &
    polarisability, bead_mass, Rayleigh_range, &
    E_tweezer, E_cavity, &
    kappa_in, kappa_nano, Gamma_M)
    ! """
    ! Calculates paramaters about the bead
    ! Specifically:
    !   polarisability, bead_mass, Rayleigh_range,
    !   E_tweezer, E_cavity,
    !   kappa_in, kappa_nano, Gamma_M
    ! """
    IMPLICIT NONE 
    integer::NPERIOD, NTOT
    double precision::R0, rho, EPSR, EPSI0
    double precision::C, hbar, BOLTZ, TEMP, Gravity
    double precision::WK, waist_radius, WX, WY, XL, Finesse, air_pressure
    double precision::tweezer_input_power, detuning, DelFSR, Theta0
    double precision::X0, Y0, Z0

    INCLUDE 'CSCAVITY.h'
    double precision::E_tweezer, E_cavity, Rayleigh_range, A, thet
    double precision::omega_optical, coeff
    double precision::KAPPA, bead_mass, kappa_nano, kappa_in
    double precision::pi, pi2
    double precision::polarisability, VOL
    double precision::Gamma_M
 
    ! zero eq. initial values
    pi = dacos(-1.d0)

    polarisability = 4.*pi * EPSI0 * (EPSR - 1.)/(EPSR + 2.) * R0**3
    bead_mass = rho * 4.*pi/3. * R0**3
    Rayleigh_range = WX * WY * WK / 2.d0

    omega_optical = C * WK
    VOL = XL * pi * waist_radius**2 / 4.d0

    E_tweezer = 4. * tweezer_input_power / (Wx * Wy * pi * c * EPSI0)
    E_tweezer = sqrt(E_tweezer)
    E_cavity = hbar * omega_optical/(2. * VOL * Epsi0)
    E_cavity = sqrt(E_cavity)
    
    coeff = WK * polarisability / EPSI0 / omega_optical**2 / pi
    
    kappa_in = pi * c / finesse / XL
    ! kappa_nano = 4. * coeff**2 * DelFSR * cos(WK * x0) * cos(WK * x0)
    ! NOTE: changing this to ensure bead paramters are accurate as x0 is sweeped through
    kappa_nano = 4. * coeff**2 * DelFSR * cos(thet) * cos(thet)

    ! take usual expression eg Levitated review by Li Geraci etc
    ! 1 bar= 10^ 5 pascal; air_pressure is in mbar = 10^ 2 pascal
    ! gamma=16 P/(pi*v*rho*R)
    ! v=speed of air=500 /s
    Gamma_M = 1600.* air_pressure / pi
    Gamma_M = Gamma_M / 500 / rho / R0
    ! Fix of Feb.2016 our Gamma_M => Gamma_M/2!!
    ! gamma/2=8 P/(pi*v*rho*R)
    Gamma_M = Gamma_M / 2.d0

    RETURN
END


SUBROUTINE CALCULATE_MECH_FREQUENCIES(&
    thet, detuning, Wkx0, &
    polarisability, bead_mass, Rayleigh_range, &
    E_tweezer, E_cavity, &
    kappa_in, kappa_nano, Gamma_M, &
    omega_x, omega_y, omega_z, &
    x0_sweep)
    ! """
    ! Calculates the mechanical frequencies for a given initial position x0_sweep
    ! """
    IMPLICIT NONE
    integer::NPERIOD, NTOT
    double precision::R0, rho, EPSR, EPSI0
    double precision::C, hbar, BOLTZ, TEMP, Gravity
    double precision::WK, waist_radius, WX, WY, XL, Finesse, air_pressure
    double precision::tweezer_input_power, detuning, DelFSR, Theta0
    double precision::X0, Y0, Z0

    INCLUDE 'CSCAVITY.h'
    double precision::thet, omega_x, omega_y, omega_z
    double precision::Edip, Ediph, Rayleigh_range, ALPre, ALPim
    double precision::E_cavity, E_tweezer, polarisability

    double precision::detuning, Gamma_M
    double precision::kappa, half_kappa, kappa_in, kappa_nano, bead_mass
    double precision::pi, pi2
    double precision::Wkx0
    double precision::C1
    double precision::x0_sweep

    PI = dacos(-1.d0)
    kappa = kappa_nano + kappa_in
    half_kappa = 0.5d0 * kappa

    omega_x = polarisability * E_tweezer**2 / bead_mass / WX**2
    omega_y = polarisability * E_tweezer**2 / bead_mass / Wy**2
    omega_z = 0.5d0 * polarisability * E_tweezer**2 / bead_mass / Rayleigh_range**2

    ! Sept 5 we will use negative Edip
    Edip = -0.5d0 * polarisability * E_tweezer * E_cavity * sin(thet)
    Ediph = Edip / hbar

    ! photon number in cavity
    ! real part of photon field
    ALPRe = detuning * Ediph * cos(Wkx0) / (half_kappa**2 + detuning**2)

    ! ADD CS POTENTIAL CORRECTION to frequency squared
    C1 = -Edip / bead_mass * 2. * ALPRe * WK**2 * cos(Wkx0)
    omega_x = omega_x + C1 * sin(thet) * sin(thet)
    omega_y = omega_y + C1 * cos(thet) * cos(thet)
    omega_z = omega_z - 2. * Edip / bead_mass * ALPre * (WK - 1.d0/Rayleigh_range)**2 * cos(Wkx0)

    omega_x = sqrt(omega_x)
    omega_y = sqrt(omega_y)
    omega_z = sqrt(omega_z)

    RETURN
END


