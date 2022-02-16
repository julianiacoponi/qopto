! This code has 1 beam and 1 dimension
! ****************************************************************************************************
! Calculates quantum noise spectra.
! Uses real mean fields
! Works out position and momentum optical quadrature
! ****************************************************************************************************


IMPLICIT NONE

integer::ii
double precision::pi, pi2

! parameter file with input values (mostly unused here)
integer::N_period, N_total
double precision::bead_radius, bead_density
double precision::epsilon_0, epsilon_R
double precision::c, hbar, kB, bath_temperature, Gravity
double precision::WK, WX, WY
double precision::waist_radius, cavity_length, Finesse, air_pressure
double precision::tweezer_input_power, detuning_Hz, DelFSR, theta_0
double precision::X0, Y0, Z0
INCLUDE 'CSCAVITY-1D.h'

double precision::theta_homodyne
integer::num_points
PARAMETER(theta_homodyne=0.d0, num_points=10000)

! Gamma_M = mechanical damping due to gas,  also noise from fluctuation dissipation
double precision::Delta, kappa, omega_M, half_Gamma_M, g_coupling
double precision::half_kappa, theta_homodyne_pi

! thermal bath phonon occupancy
double precision::nbar, num_photons
double precision::num_phonons, phonons_from_optical_noise, oscillator_temperature

double precision::omega_value, max_omega, omega_range, omega_increment
double precision::S_XX_value, S_heterodyne_value, S_homodyne_value, S_homodyne_negative 
double precision::area_under_S_XX, area_under_S_heterodyne, area_under_S_homodyne
double precision, DIMENSION(num_points)::S_XX_array, S_heterodyne_array, S_homodyne_array

! for optomechanical cooling formula
double precision::C_plus, C_minus, x_cooling, x_phonons_from_cooling_formula

! which expressions to calculate spectra with
integer::expressions_choice
! which input parameters to use
integer::parameters_choice

WRITE(*, *)  "Choose input parameters: "
WRITE(*, *) '1 = original parameters'
WRITE(*, *) '2 = artificial parameters'
READ(*, *)  parameters_choice

WRITE(*, *) 'Choose expressions to calculate noise vectors with:'
WRITE(*, *) '1 = original expressions (in original scripts)'
WRITE(*, *) '2 = slightly reworked expressions (using back-action susceptibilities)'
WRITE(*, *) '3 = my own manually-derived expressions'
READ(*, *)  expressions_choice

! FORMAT can be understood from:
! https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
! http://www.personal.psu.edu/jhm/f90/lectures/23.html

! PRINT, OPEN, WRITE, CLOSE can be understood from:
! https://www.tutorialspoint.com/fortran/fortran_file_input_output.htm
! https://meteor.geol.iastate.edu/classes/mt227/lectures/Formatting.pdf


! e.g. the below is:
! '4' = 4 things written to the same line
! 'ES' = scientific exponential notation (1-9E+power, instead of 0.1-0.9)
! '23' = number of positions to be used
! '15' = number of digits to the right of the decimal point
! NOTE: this is *not* the same as the precision of the number
101 FORMAT(4ES23.15)
305 FORMAT(5ES11.3)


pi = dacos(-1.d0)
pi2 = 2.d0 * pi
theta_homodyne_pi = theta_homodyne * pi

! input parameters
IF (parameters_choice == 1) THEN
    ! formerly "FTanX.dat" - saves: omega, S_XX, S_heterodyne, S_homodyne
    OPEN(1, file="spectra.dat", status="unknown")
    OPEN(2, file="spectra-vs-omega-not-Hz.dat", status="unknown")
    Delta = -200000.0 * pi2 
    kappa = 438032.94007709384 * pi2
    omega_M = 128062.66787 * pi2
    half_Gamma_M = 0.324068D-02
    g_coupling = 53061.566105502839 * pi2
    WRITE(*, *) 'Delta, kappa, omega_M, Gamma_M, g_coupling, n = nbar * omega_M = (kB * bath_temperature)/hbar'
    WRITE(*, 500) Delta, kappa, omega_M, half_Gamma_M * 2, g_coupling, (kB * bath_temperature)/hbar

    ! thermal bath occupancy
    nbar = (kB * bath_temperature)/(hbar * omega_M)
    WRITE(*, *) 'nbar = (kB*bath_temperature)/(hbar*omega_M) = thermal bath phonon occupancy:'
    WRITE(*, 500) nbar

! dummy parameters from Mathematica
! TODO: remove these dummy parameters
ELSE IF (parameters_choice == 2) THEN
    OPEN(1, file="spectra-artificial.dat", status="unknown")
    OPEN(2, file="spectra-artificial-vs-omega-not-Hz.dat", status="unknown")
    Delta = -3.d0
    kappa = 4.5d0
    omega_M = 6.d0
    half_Gamma_M = 2.5d0 / 2.
    g_coupling = 0.25d0
    nbar = 10.d0/omega_M
    WRITE(*, *) 'Delta, kappa, omega_M, Gamma_M, g_coupling, n = nbar * omega_M = (kB * bath_temperature)/hbar'
    WRITE(*, 305) Delta, kappa, omega_M, half_Gamma_M * 2, g_coupling, nbar * omega_M

END IF
! give the .dat files column titles
WRITE(1, *) "omega/2π ", "S_XX ", "S_heterodyne ", "S_homodyne"
WRITE(2, *) "omega ", "S_XX ", "S_heterodyne ", "S_homodyne"

half_kappa = kappa * 0.5d0

! shot noise
num_photons = 0.d0

! estimate x-axis phonons using optomechanical cooling formula for comparison
C_plus = ((Delta + omega_M)**2 + half_kappa**2)**(-1)
C_minus = ((Delta - omega_M)**2 + half_kappa**2)**(-1)
x_cooling = -(g_coupling**2 * kappa) * (C_plus - C_minus)
x_phonons_from_cooling_formula = nbar * 2 * half_Gamma_M/(abs(x_cooling) + 2. * half_Gamma_M)

! open loop over frequency sweep for noise spectra
! e.g. S_XX(omega) = FT(autocorrelation <X(t)X^T(t+tau)>)

! expecting peaks at omega_M, so sample twice that range
max_omega = 2.d0 * omega_M * 1.001
! sampling both negative and positive omega
omega_range = 2.d0 * max_omega
omega_increment = omega_range/num_points

DO 11 ii = 1, num_points

    omega_value = -max_omega + (ii - 1) * omega_increment

    ! work out Power Spectral Density (PSD) of the homodyne, S_homodyne      
    ! work out PSDs for negative frequency for symmetrisation of homodyne
    ! NOTE: the call below is used for S_homodyne and S_heterodyne, not S_XX
    CALL CALCULATE_SPECTRA(&
        expressions_choice, &
        g_coupling, nbar, num_photons, theta_homodyne_pi, &
        Delta, half_kappa, omega_M, half_Gamma_M, &
        -omega_value, S_XX_value, S_homodyne_value, S_heterodyne_value)

    ! update homodyne
    S_homodyne_negative = S_homodyne_value

    ! work out the same PSDs for positive frequencies
    CALL CALCULATE_SPECTRA(&
        expressions_choice, &
        g_coupling, nbar, num_photons, theta_homodyne_pi, &
        Delta, half_kappa, omega_M, half_Gamma_M, &
        omega_value, S_XX_value, S_homodyne_value, S_heterodyne_value)

    S_XX_array(ii) = S_XX_value
    S_heterodyne_array(ii) = S_heterodyne_value

    ! symmetrise S_homodyne
    S_homodyne_value = 0.5 * (S_homodyne_negative + S_homodyne_value)
    ! to find optimal squeezing
    S_homodyne_array(ii) = S_homodyne_value

    ! save omega, S_XX, S_aa (S_heterodyne), S_(Chi_theta.Chi_theta) (S_homodyne), to spectra.dat
    ! Chi_theta = a * exp(-i*theta) + a^dagg * exp(+i*theta)

    ! NOTE: saving omega/2π here means the x-axis for any plots is frequency (Hz), not angular frequency
    WRITE(1, 101) omega_value/pi2, S_XX_value, S_heterodyne_value, S_homodyne_value
    WRITE(2, 101) omega_value, S_XX_value, S_heterodyne_value, S_homodyne_value

11 END DO

! back action limit: at zero pressure, this is the residual phonon number from optical noise
! this should be the theoretical limit if Gamma_M is set to zero
! NOTE: assumes Delta is negative?
phonons_from_optical_noise = ((omega_M + Delta)**2 + half_kappa**2)/(-4 * omega_M * Delta)
WRITE(*, *) 'Back-action limited phonons:'
WRITE(*, 500) phonons_from_optical_noise

! integrate the PSD for the (rescaled) position operator X^hat = (b + b^dagg)
! uses trapezoidal rule
CALL INTEGRATE_SPECTRUM(num_points, omega_increment, S_XX_array, area_under_S_XX)
WRITE(*, *) 'Area under S_XX:'
WRITE(*, 500) area_under_S_XX

CALL INTEGRATE_SPECTRUM(num_points, omega_increment, S_heterodyne_array, area_under_S_heterodyne)
WRITE(*, *) 'Area under S_heterodyne:'
WRITE(*, 500) area_under_S_heterodyne

CALL INTEGRATE_SPECTRUM(num_points, omega_increment, S_homodyne_array, area_under_S_homodyne)
WRITE(*, *) 'Area under S_homodyne:'
WRITE(*, 500) area_under_S_homodyne


! area_under_S_XX corresponds to 2n+1 so convert to get n = phonon occupancy of the oscillator
num_phonons = 0.5d0 * (area_under_S_XX - 1.d0)

! NOTE: original code goes straight from area "XRE" (2n+1) to get TEMP
! do more precisely with num_phonons (n) here
oscillator_temperature = num_phonons * hbar * omega_M/kB

WRITE(*, *) 'Phonons from optomechancial cooling formula:'
WRITE(*, 500) x_phonons_from_cooling_formula
WRITE(*, *) 'Phonons from S_XX sideband asymmetry:'
WRITE(*, 500) num_phonons
WRITE(*, *) 'Mechanical oscillator`s temperature /K:'
WRITE(*, 500) oscillator_temperature

500 FORMAT(ES14.5)
STOP
END

! Formerly 'ANALYT'
SUBROUTINE CALCULATE_SPECTRA(&
    expressions_choice, &
    g_coupling, nbar, num_photons, theta_homodyne_pi, &
    Delta, half_kappa, omega_M, half_Gamma_M, &
    omega_value, S_XX_value, S_homodyne_value, S_heterodyne_value)
    ! """
    ! Generic routine for calculating the noise spectra of the trap and probe beams
    ! """
    IMPLICIT NONE
    integer::N_total, expressions_choice
    PARAMETER(N_total=4)
    double precision::g_coupling, nbar, num_photons, theta_homodyne_pi
    double precision::omega_value, S_XX_value, S_homodyne_value, S_heterodyne_value
    double precision::Delta, half_kappa, omega_M, half_Gamma_M
    double complex::chi_C, chi_C_star, chi_M, chi_M_star
    double complex::BAX(1, N_total), A1(1, N_total), A1dagg(1, N_total)

    ! My alternatives with just the 1D arrays needed
    double complex, DIMENSION(N_total)::X_hat_vector, a_vector, a_dagg_vector

    ! make sure to reset this to zero
    S_XX_value = 0.d0

    ! first get the susceptibilities
    CALL CALCULATE_SUSCEPTIBILITIES(&
        omega_value, omega_M, &
        Delta, half_kappa, half_Gamma_M, &
        chi_C, chi_C_star, chi_M, chi_M_star)

    ! work out noise vectors for X^hat, a and a^dagg - to get S_XX values
    IF (expressions_choice == 1) THEN
        CALL CALCULATE_NOISE_VECTORS(&
            g_coupling, half_kappa, half_Gamma_M, &
            chi_C, chi_C_star, chi_M, chi_M_star, &
            BAX, A1, A1dagg)

        ! NOTE: this ignores num_photons (which is set to zero for now anyway)
        S_XX_value = (abs(BAX(1, 1)))**2
        S_XX_value = S_XX_value + (abs(BAX(1, 3)))**2 * (nbar + 1)
        S_XX_value = S_XX_value + (abs(BAX(1, 4)))**2 * nbar

    ELSE IF (expressions_choice == 2) THEN
        ! my alternative using the back-action expressions
        ! NOTE: putting these expressions in Mathematica match exactly with my manually derived ones
        CALL CALCULATE_NOISE_VECTORS_VIA_BACK_ACTIONS(&
            g_coupling, half_kappa, half_Gamma_M, &
            chi_C, chi_C_star, chi_M, chi_M_star, &
            X_hat_vector, a_vector, a_dagg_vector)

        S_XX_value = S_XX_value + (abs(X_hat_vector(1)))**2 * (num_photons + 1)
        S_XX_value = S_XX_value + (abs(X_hat_vector(2)))**2 * num_photons
        S_XX_value = S_XX_value + (abs(X_hat_vector(3)))**2 * (nbar + 1)
        S_XX_value = S_XX_value + (abs(X_hat_vector(4)))**2 * nbar

    ELSE IF (expressions_choice == 3) THEN
        ! my alternative using manually derived expressions
        CALL CALCULATE_NOISE_VECTORS_VIA_MANUAL_DERIVATIONS(&
            g_coupling, half_kappa, half_Gamma_M, &
            chi_C, chi_C_star, chi_M, chi_M_star, &
            X_hat_vector, a_vector, a_dagg_vector)

        S_XX_value = S_XX_value + (abs(X_hat_vector(1)))**2 * (num_photons + 1)
        S_XX_value = S_XX_value + (abs(X_hat_vector(2)))**2 * num_photons
        S_XX_value = S_XX_value + (abs(X_hat_vector(3)))**2 * (nbar + 1)
        S_XX_value = S_XX_value + (abs(X_hat_vector(4)))**2 * nbar

    END IF

    ! NOTE: Before calculating the heterodyne spectrum, the original script allows for shifting
    ! omega by a certain "heterodyne-beat-shifted" amount to calculate new susceptiblity values,
    ! which are then used to calculate new noise spectra. These new noise spectra are then used to
    ! calculate the heterodyne spectrum. This is not needed for now.

    ! work out the homodyne and heterodyne and homodyne spectra using same vectors
    IF (expressions_choice == 1) THEN
        CALL HETERODYNE(N_total, nbar, num_photons, A1dagg, S_heterodyne_value)
        CALL HOMODYNE(N_total, nbar, num_photons, theta_homodyne_pi, A1, A1dagg, S_homodyne_value)
    ELSE
        CALL HETERODYNE_2(N_total, nbar, num_photons, a_vector, a_dagg_vector, S_heterodyne_value, 1)
        CALL HOMODYNE_2(N_total, nbar, num_photons, theta_homodyne_pi, a_vector, a_dagg_vector, S_homodyne_value)
    END IF

    ! TODO: understand derivation of this expression
    ! These expression rely on theta_homodyne = 0
    ! - minus 1 is also due to the input-output relation
    ! - why only chi_C_star and not chi_C - chi_C_star like for the homodyne?
    S_heterodyne_value = (S_heterodyne_value - 1.d0)/abs(chi_C_star)**2/g_coupling**2/half_kappa/2.

    ! NOTE: minus 1 here is due to the input-output relation
    S_homodyne_value = (S_homodyne_value - 1.d0)/abs(chi_C - chi_C_star)**2/g_coupling**2/half_kappa/2.

    RETURN
END

SUBROUTINE HOMODYNE(N_total, nbar, num_photons, theta_homodyne_pi, A1, A1dagg, S_homodyne)
    ! """
    ! Work out homodyne spectra
    ! """
    IMPLICIT NONE
    integer::ii, N_total
    double precision::theta_homodyne_pi, theta, S_homodyne, nbar, num_photons
    double complex::imag_n
    double complex::homodyne_noise_vector(1, N_total)
    double complex::PY_hat_vector(1, N_total), Y_hat_vector(1, N_total)
    double complex::A1(1, N_total)
    double complex::A1dagg(1, N_total)

    theta = theta_homodyne_pi
    ! i
    imag_n = cmplx(0.d0, 1.d0)
    DO ii = 1, N_total
        ! optical field "amplitude" quadrature (photon position operator)
        ! Y^hat = (a^dagg + a)
        Y_hat_vector(1, ii) = A1(1, ii) + A1dagg(1, ii)

        ! optical field "phase" quadratrue (photon momentum operator)
        ! PY^hat = i*(a^dagg - a)
        ! NOTE: this way around (A1 - A1dagg) just defines the momentum in the opposite direction
        PY_hat_vector(1, ii) = imag_n * (A1(1, ii) - A1dagg(1, ii))

        ! since theta_homodyne_pi is set to 0, homodyne_noise_vector = Y_hat_vector
        ! Therefore S_@@ = S_YY (@=theta)
        homodyne_noise_vector(1, ii) = PY_hat_vector(1, ii) * sin(theta) + Y_hat_vector(1, ii) * cos(theta)
    END DO

    S_homodyne = 0.d0
    S_homodyne = S_homodyne + abs(homodyne_noise_vector(1, 1))**2 * num_photons
    S_homodyne = S_homodyne + abs(homodyne_noise_vector(1, 2))**2 * (num_photons + 1.d0)
    S_homodyne = S_homodyne + abs(homodyne_noise_vector(1, 3))**2 * nbar
    S_homodyne = S_homodyne + abs(homodyne_noise_vector(1, 4))**2 * (nbar + 1.d0)

    RETURN
END


! My own version with 1D vectors
SUBROUTINE HOMODYNE_2(N_total, nbar, num_photons, theta_homodyne_pi, a_vector, a_dagg_vector, S_homodyne)
    ! """
    ! Work out homodyne spectra
    ! """
    IMPLICIT NONE
    integer::ii, N_total
    double precision:: theta_homodyne_pi, S_homodyne, nbar, num_photons
    double complex::imag_n
    double complex::homodyne_noise_vector(N_total), PY_hat_vector(N_total), Y_hat_vector(N_total)
    double complex::a_vector(N_total), a_dagg_vector(N_total)

    ! i
    imag_n = cmplx(0.d0, 1.d0)
    DO ii = 1, N_total
        ! optical field "amplitude" quadrature (photon position operator)
        ! Y^hat = (a^dagg + a)
        Y_hat_vector(ii) = a_dagg_vector(ii) + a_vector(ii)

        ! optical field "phase" quadratrue (photon momentum operator)
        ! PY^hat = i*(a^dagg - a)
        ! NOTE: the original has this other way around, which just defines the momentum in
        ! the opposite dreiction
        PY_hat_vector(ii) = imag_n * (a_dagg_vector(ii) - a_vector(ii))

        ! since theta_homodyne_pi is set to 0, homodyne_noise_vector = Y_hat_vector
        ! Therefore S_@@ = S_YY (@=theta)
        homodyne_noise_vector(ii) = PY_hat_vector(ii) * sin(theta_homodyne_pi) + Y_hat_vector(ii) * cos(theta_homodyne_pi)
    END DO

    S_homodyne = 0.d0
    S_homodyne = S_homodyne + abs(homodyne_noise_vector(1))**2 * num_photons
    S_homodyne = S_homodyne + abs(homodyne_noise_vector(2))**2 * (num_photons + 1.d0)
    S_homodyne = S_homodyne + abs(homodyne_noise_vector(3))**2 * nbar
    S_homodyne = S_homodyne + abs(homodyne_noise_vector(4))**2 * (nbar + 1.d0)

    RETURN
END


SUBROUTINE HETERODYNE(N_total, nbar, num_photons, A1dagg, S_heterodyne)
    IMPLICIT NONE
    integer::ii, N_total
    double precision::S_heterodyne, nbar, num_photons
    double complex::imag_n
    double complex::A1dagg(1, N_total), heterodyne_noise_vector(1, N_total)
    
    ! (0 + i)
    imag_n = cmplx(0.d0, 1.d0)
    DO ii = 1, N_total
        heterodyne_noise_vector(1, ii) = A1dagg(1, ii)
    END DO

    S_heterodyne = 0.d0
    S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(1, 1))**2 * num_photons
    S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(1, 2))**2 * (num_photons + 1.d0)
    S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(1, 3))**2 * nbar
    S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(1, 4))**2 * (nbar + 1.d0)
    
    RETURN
END

! My own version with 1D vectors
SUBROUTINE HETERODYNE_2(N_total, nbar, num_photons, a_vector, a_dagg_vector, S_heterodyne, het_is_a_dagg)
    IMPLICIT NONE
    integer::ii, N_total, het_is_a_dagg
    double precision::S_heterodyne, nbar, num_photons
    double complex::imag_n
    double complex, DIMENSION(N_total)::a_vector, a_dagg_vector, heterodyne_noise_vector
    
    ! (0 + i)
    imag_n = cmplx(0.d0, 1.d0)
    DO ii = 1, N_total
        IF (het_is_a_dagg == 0) THEN
            heterodyne_noise_vector(ii) = a_vector(ii)
        ELSE
            heterodyne_noise_vector(ii) = a_dagg_vector(ii)
        END IF
    END DO

    S_heterodyne = 0.d0
    S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(1))**2 * num_photons
    S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(2))**2 * (num_photons + 1.d0)
    S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(3))**2 * nbar
    S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(4))**2 * (nbar + 1.d0)
    RETURN
END


! Formerly 'SUSCEPT'
SUBROUTINE CALCULATE_SUSCEPTIBILITIES(&
    omega, omega_M, &
    Delta, half_kappa, half_Gamma_M, &
    chi_C, chi_C_star, chi_M, chi_M_star)
    ! """
    ! Calculates the optical (chi_C, chi_C_star) and mechanical (chi_M, chi_M_star) susceptibilities
    ! """
    IMPLICIT NONE 
    double precision::omega, omega_M
    double precision::Delta, half_Gamma_M, half_kappa
    double complex::chi_C, chi_C_star, chi_M, chi_M_star
    double precision::chi_real, chi_imag, chi_denominator

    ! optical susceptibilities
    ! chi_C
    chi_denominator = half_kappa**2 + (omega + Delta)**2
    chi_real = half_kappa/chi_denominator
    chi_imag = (omega + Delta)/chi_denominator
    chi_C = cmplx(chi_real, chi_imag)

    ! chi_C^*(-omega) - i.e. turn omega into -omega, then chi_imag into -chi_imag.
    chi_denominator = half_kappa**2 + (-omega + Delta)**2
    chi_real = half_kappa/chi_denominator
    chi_imag = (-omega + Delta)/chi_denominator
    chi_C_star = cmplx(chi_real, -chi_imag)

    ! mechanical susceptibilities
    ! chi_M
    chi_denominator = (half_Gamma_M)**2 + (omega - omega_M)**2
    chi_real =  half_Gamma_M/chi_denominator
    chi_imag = (omega - omega_M)/chi_denominator
    chi_M = cmplx(chi_real, chi_imag)

    ! chi_M*(-omega) - i.e. turn omega into -omega, then chi_imag into -chi_imag.
    chi_denominator = (half_Gamma_M)**2 + (-omega - omega_M)**2
    chi_real =  half_Gamma_M/chi_denominator
    chi_imag = (-omega - omega_M)/chi_denominator
    chi_M_star = cmplx(chi_real, -chi_imag)

    RETURN
END

! Formerly 'Avect'
SUBROUTINE CALCULATE_NOISE_VECTORS(&
    g_coupling, half_kappa, half_Gamma_M, &
    chi_C, chi_C_star, chi_M, chi_M_star, &
    BAX, A1, A1dagg)
    ! """
    ! Work out vectors (in terms of noise operators) of:
    ! - the (rescaled) position operator X^hat and,
    ! - the optical field operators (a, a^dagg)
    ! """
    IMPLICIT NONE
    integer::ii, N_total
    double precision::Sqrtkapp
    double precision::g_coupling, half_kappa, half_Gamma_M, sqrt_Gamma_M
    PARAMETER(N_total=4)
    double complex::chi_C, chi_C_star
    double complex::chi_M, chi_M_star
    double complex::one, imag_n

    double complex::CMX, eta_C, eta_M, CA1, CA1dagg, BETX

    double complex, INTENT(OUT)::BAX(1, N_total)
    double complex, INTENT(OUT)::A1(1, N_total)
    double complex, INTENT(OUT)::A1dagg(1, N_total)
    double complex::N0X(N_total)    

    ! i
    imag_n = cmplx(0.d0, 1.d0)
    one = cmplx(1.d0, 0.0d0)

    ! Susceptibilities
    eta_C = chi_C - chi_C_star
    eta_M = chi_M - chi_M_star

    ! Normalisations
    CMX = 1.d0 + g_coupling**2 * eta_C * eta_M
    Sqrtkapp = sqrt(2.d0 * half_kappa)
    sqrt_Gamma_M = sqrt(2.d0 * half_Gamma_M)
    BETX = imag_n * g_coupling * eta_M * Sqrtkapp
    
    ! X^hat noise vector; weights of a_in, a_in^dagg, b_in, b_in^dagg
    ! f_3 i.e. for a_in
    N0X(1) = BETX * chi_C/CMX
    ! f_4 i.e. for a_in^dagg
    N0X(2) = BETX * chi_C_star/CMX
    ! f_1 i.e. for b_in
    N0X(3) = sqrt_Gamma_M * chi_M/CMX
    ! f_2 i.e. for b_in^dagg
    N0X(4) = sqrt_Gamma_M * chi_M_star/CMX  
    
    ! FOR 1D
    DO ii = 1, N_total
        BAX(1, ii) = N0X(ii) 
    END DO

    ! now work out the optical trap field, a
    CA1 = imag_n * chi_C
    ! now work out the photon field, a^dagg
    CA1dagg = -imag_n * chi_C_star

    DO ii = 1, N_total
        A1(1, ii) = CA1 * g_coupling * BAX(1, ii)
        A1dagg(1, ii) = CA1dagg * g_coupling * BAX(1, ii)
    END DO

    ! add shot or incoming noise
    ! trap beam: add cavity-filtered contribution
    A1(1, 1) = A1(1, 1) + Sqrtkapp * chi_C
    A1dagg(1, 2) = A1dagg(1, 2) + Sqrtkapp * chi_C_star
    
    ! cavity output : add incoming imprecision
    ! work out a_out = a_in - Sqrtkapp * a
    DO ii = 1, N_total
        A1(1, ii) = -A1(1, ii) * Sqrtkapp
        A1dagg(1, ii) = -A1dagg(1, ii) * Sqrtkapp
    END DO
    
    A1(1, 1) = one + A1(1, 1)
    A1dagg(1, 2) = one + A1dagg(1, 2)
    
    RETURN
END

SUBROUTINE CALCULATE_NOISE_VECTORS_VIA_BACK_ACTIONS(&
    g_coupling, half_kappa, half_Gamma_M, &
    chi_C, chi_C_star, chi_M, chi_M_star, &
    X_hat_vector, a_vector, a_dagg_vector)
    ! """
    ! Work out vectors (in terms of noise operators) of:
    ! - the (rescaled) position operator X^hat and,
    ! - the optical field operators (a, a^dagg)
    ! (alternative but equivalent expressions using the back-actions)
    ! """
    IMPLICIT NONE
    integer::ii, N_total
    double precision::g_coupling, half_kappa, half_Gamma_M, sqrt_kappa, sqrt_Gamma_M
    PARAMETER(N_total=4)
    double complex::one, imag_n, ig, nu
    double complex::chi_C, chi_C_star, chi_CBA
    double complex::chi_M, chi_M_star, chi_MBA

    double complex::a_vector(N_total)
    double complex::a_dagg_vector(N_total)
    double complex::X_hat_vector(N_total)

    ! (0 + i)
    imag_n = cmplx(0.d0, 1.d0)
    ! (1 + 0i)
    one = cmplx(1.d0, 0.0d0)

    ! i * g_coupling - appears a lot in the derivations
    ig = imag_n * g_coupling

    ! optical back-action
    ! susceptibility of optical field to back-action mechanical pressure fluctuations
    chi_CBA = ig * (chi_C - chi_C_star)

    ! mechanical back-action
    ! susceptibility of mechanical oscillator to back-action radiation pressure fluctuations
    chi_MBA = ig * (chi_M - chi_M_star)

    ! copying CMX = 1.d0 + g_coupling**2 * eta_M * eta_C
    ! = 1 - (ig)**2 * eta_M * eta_C
    ! where ig * eta_M = chi_MBA, and ig * eta_C = chi_CBA
    ! so nu = 1/CMX
    nu = 1.d0/(1.d0 - chi_CBA * chi_MBA)
    sqrt_kappa = sqrt(2.d0 * half_kappa)
    sqrt_Gamma_M = sqrt(2.d0 * half_Gamma_M)

    ! X^hat noise vector; weights of a_in, a_in^dagg, b_in, b_in^dagg
    ! f_3 = k_M3 + k_M4_star i.e. for a_in
    X_hat_vector(1) = nu * sqrt_kappa * chi_MBA * chi_C
    ! f_4 = k_M4 + k_M3_star i.e. for a_in^dagg
    X_hat_vector(2) = nu * sqrt_kappa * chi_MBA * chi_C_star
    ! f_1 = k_M1 + k_M2_star i.e. for b_in
    X_hat_vector(3) = nu * sqrt_Gamma_M * chi_M
    ! f_2 = k_M2 + k_M1_star i.e. for b_in^dagg
    X_hat_vector(4) = nu * sqrt_Gamma_M * chi_M_star

    ! calculate a, a^dagg by exploiting:
    ! a = chi_C * (ig * X + sqrt_kappa * a_in)
    ! a^dagg = chi_C_star * (-ig * X + sqrt_kappa * a_in^dagg)
    DO ii = 1, N_total
        ! now work out the optical trap field, a
        a_vector(ii) = ig * chi_C * X_hat_vector(ii)
        ! now work out the photon field, a^dagg
        a_dagg_vector(ii) = -ig * chi_C_star * X_hat_vector(ii)
    END DO

    ! add shot or incoming noise
    ! trap beam: add cavity-filtered contribution
    a_vector(1) = a_vector(1) + sqrt_kappa * chi_C
    a_dagg_vector(2) = a_dagg_vector(2) + sqrt_kappa * chi_C_star

    ! cavity output: add incoming imprecision
    ! use input-output relations:
    ! a_out = a_in - sqrt_kappa * a
    ! a_out^dagg = a_in^dagg - sqrt_kappa * a^dagg
    DO ii = 1, N_total
        a_vector(ii) = -a_vector(ii) * sqrt_kappa
        a_dagg_vector(ii) = -a_dagg_vector(ii) * sqrt_kappa
    END DO

    a_vector(1) = one + a_vector(1)
    a_dagg_vector(2) = one + a_dagg_vector(2)

    RETURN
END


SUBROUTINE CALCULATE_NOISE_VECTORS_VIA_MANUAL_DERIVATIONS(&
    g_coupling, half_kappa, half_Gamma_M, &
    chi_C, chi_C_star, chi_M, chi_M_star, &
    X_hat_vector, a_vector, a_dagg_vector)
    ! """
    ! Work out vectors (in terms of noise operators) of:
    ! - the (rescaled) position operator X^hat and,
    ! - the optical field operators (a, a^dagg)
    ! (alternative but equivalent expressions I derived manually)
    ! """
    IMPLICIT NONE
    integer::ii, N_total
    double precision::g_coupling, half_kappa, half_Gamma_M, sqrt_kappa, sqrt_Gamma_M
    PARAMETER(N_total=4)
    double complex::one, imag_n, ig
    double complex::chi_C, chi_C_star, chi_CBA
    double complex::chi_M, chi_M_star, chi_MBA

    ! expressions needed for a, a^dagg
    double complex::A_C, A_C_star
    double complex::P_C, P_C_star
    double complex::Q_C
    double complex::B_C
    double complex::k_C1, k_C1_star
    double complex::k_C2, k_C2_star
    double complex::k_C3, k_C3_star
    double complex::k_C4, k_C4_star

    ! expressions needed for b, b^dagg (hence X)
    double complex::A_M, A_M_star
    double complex::P_M, P_M_star
    double complex::Q_M
    double complex::B_M
    double complex::k_M1, k_M1_star
    double complex::k_M2, k_M2_star
    double complex::k_M3, k_M3_star
    double complex::k_M4, k_M4_star

    double complex::a_vector(N_total)
    double complex::a_dagg_vector(N_total)
    double complex::X_hat_vector(N_total)

    ! (0 + i)
    imag_n = cmplx(0.d0, 1.d0)
    ! (1 + 0i)
    one = cmplx(1.d0, 0.0d0)

    ! i * g_coupling - appears a lot in the derivations
    ig = imag_n * g_coupling

    sqrt_kappa = sqrt(2.d0 * half_kappa)
    sqrt_Gamma_M = sqrt(2.d0 * half_Gamma_M)

    ! optical back-action
    ! susceptibility of optical field to back-action thermal pressure fluctuations
    chi_CBA = ig * (chi_C - chi_C_star)

    ! mechanical back-action
    ! susceptibility of mechanical oscillator to back-action radiation pressure fluctuations
    chi_MBA = ig * (chi_M - chi_M_star)

    ! NOTE: my manual derivation assumed the following:
    ! a = k_C1.a_in + k_C2.a_in^dagg + k_C3.b_in + k_C4.b_in^dagg
    ! b = k_M1.b_in + k_M2.b_in^dagg + k_M3.a_in + k_M4.a_in^dagg
    ! a^dagg = k_C2_star.a_in + k_C1_star.a_in^dagg + k_C4_star.b_in + k_C3_star.b_in^dagg
    ! b^dagg = k_M2_star.b_in + k_M1_star.b_in^dagg + k_M4_star.a_in + k_M3_star.a_in^dagg
    ! because this left the coefficients symmetrical in {a,C,kappa}<-->{b,M,Gamma_M}

    ! define weights for a, a^dagg first
    A_C = 1.d0/(1.d0 - ig * chi_C * chi_MBA)
    A_C_star = 1.d0/(1.d0 + ig * chi_C_star * chi_MBA)

    P_C = A_C * chi_C
    P_C_star = A_C_star * chi_C_star

    Q_C = ig * P_C * P_C_star * chi_MBA
    ! Q_C_star = -Q_C

    B_C = 1.d0/(1.d0 + ig * Q_C * chi_MBA)
    ! B_C_star = B_C

    k_C1 = B_C * P_C * sqrt_kappa
    k_C2 = B_C * Q_C * sqrt_kappa
    k_C3 = ig * B_C * (P_C - Q_C) * chi_M * sqrt_Gamma_M
    k_C4 = ig * B_C * (P_C - Q_C) * chi_M_star * sqrt_Gamma_M

    k_C1_star = B_C * P_C_star * sqrt_kappa
    k_C2_star = -B_C * Q_C * sqrt_kappa
    k_C3_star = -ig * B_C * (P_C_star + Q_C) * chi_M_star * sqrt_Gamma_M
    k_C4_star = -ig * B_C * (P_C_star + Q_C) * chi_M * sqrt_Gamma_M

    a_vector(1) = k_C1
    a_vector(2) = k_C2
    a_vector(3) = k_C3
    a_vector(4) = k_C4

    a_dagg_vector(1) = k_C2_star
    a_dagg_vector(2) = k_C1_star
    a_dagg_vector(3) = k_C4_star
    a_dagg_vector(4) = k_C3_star

    ! find b, b^dagg
    A_M = 1.d0/(1.d0 - ig * chi_M * chi_CBA)
    A_M_star = 1.d0/(1.d0 + ig * chi_M_star * chi_CBA)

    P_M = A_M * chi_M
    P_M_star = A_M_star * chi_M_star

    Q_M = ig * P_M * P_M_star * chi_CBA
    ! Q_M_star = -Q_M

    B_M = 1.d0/(1.d0 + ig * Q_M * chi_CBA)
    ! B_M_star = B_M

    k_M1 = B_M * P_M * sqrt_Gamma_M
    k_M2 = B_M * Q_M * sqrt_Gamma_M
    k_M3 = ig * B_M * (P_M - Q_M) * chi_C * sqrt_kappa
    k_M4 = ig * B_M * (P_M - Q_M) * chi_C_star * sqrt_kappa

    k_M1_star = B_M * P_M_star * sqrt_Gamma_M
    k_M2_star = -B_M * Q_M * sqrt_Gamma_M
    k_M3_star = -ig * B_M * (P_M_star + Q_M) * chi_C_star * sqrt_kappa
    k_M4_star = -ig * B_M * (P_M_star + Q_M) * chi_C * sqrt_kappa

    ! X^hat = b + b^dagg
    ! X^hat = f_1.b_in + f_2.b_in^dagg + f_3.a_in + f_4.a_in^dagg
    ! f_3 = k_M3 + k_M4_star i.e. for a_in
    X_hat_vector(1) = k_M3 + k_M4_star
    ! f_4 = k_M4 + k_M3_star i.e. for a_in^dagg
    X_hat_vector(2) = k_M4 + k_M3_star
    ! f_1 = k_M1 + k_M2_star i.e. for b_in
    X_hat_vector(3) = k_M1 + k_M2_star
    ! f_2 = k_M2 + k_M1_star i.e. for b_in^dagg
    X_hat_vector(4) = k_M2 + k_M1_star
    
    ! cavity output: add incoming imprecision
    ! use input-output relations to get a final a_out used for heterodyne and homodyne:
    ! a_out = a_in - sqrt_kappa * a
    ! a_out^dagg = a_in^dagg - sqrt_kappa * a^dagg
    DO ii = 1, N_total
        a_vector(ii) = -a_vector(ii) * sqrt_kappa
        a_dagg_vector(ii) = -a_dagg_vector(ii) * sqrt_kappa
    END DO
    
    a_vector(1) = one + a_vector(1)
    a_dagg_vector(2) = one + a_dagg_vector(2)

    RETURN
END


! Formerly 'NORM1'
SUBROUTINE INTEGRATE_SPECTRUM(num_points, omega_increment, spectrum_array, area_under_spectrum)
    IMPLICIT NONE 
    integer::ii, num_points
    double precision::omega_increment, area_under_spectrum, trapeze_base, pi2
    double precision::spectrum_array(num_points)
    
    pi2 = 2.d0 * dacos(-1.d0)

    ! integrate the position spectrum of bead
    ! quick hack - use trapezoidal rule - improve later?
    ! each trapeze base is the omega_increment i.e. omega_range/num_points 

    ! NOTE: divide by 2π to convert omega to non-angular frequency (Hz)
    ! TODO: explian the division by 2π as a result of Parseval's theorem / "Sum-rule" ?
    ! this means the area is the area under the spectrum curve plotted vs. Hz
    trapeze_base = omega_increment/pi2

    ! different variables are used, but reset this here just to make sure
    area_under_spectrum = 0.d0
    DO ii = 1, num_points - 1
        area_under_spectrum = area_under_spectrum + 0.5d0 * (spectrum_array(ii) + spectrum_array(ii + 1))
    END DO
    area_under_spectrum = area_under_spectrum * trapeze_base

    RETURN
END