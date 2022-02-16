! This code has 1 beam plus 3 dimensions
! Calculates quantum noise spectra.
! USES real mean fields
! works out position and momentum optical quadrature
! Gamma_M = mechanical damping due to gas, also noise from fluctuation dissipation.


IMPLICIT NONE

! parameter file with input values
integer::N_period, N_total
double precision::bead_radius, bead_density
double precision::epsilon_0, epsilon_R
double precision::c, hbar, kB, bath_temperature, Gravity
double precision::WK, WX, WY
double precision::waist_radius, cavity_length, Finesse, air_pressure
double precision::tweezer_input_power, detuning_Hz, DelFSR, theta_0
double precision::X0, Y0, Z0
INCLUDE 'CSCAVITY.h'

integer::num_points
double precision::theta_homodyne_0
PARAMETER(num_points=10000, theta_homodyne_0=0.12d0)

integer::ii, jj, nn
double precision::pi, pi2
double precision::cavity_photons
double precision::theta, WKX0
double precision::G_matrix(6)
double precision::Rayleigh_range, bead_mass, polarisability
double precision::E_tweezer, E_cavity
double precision::detuning, kappa, Gamma_M
double precision::omega_x, omega_y, omega_z
double precision::max_omega, omega_increment, omega_sweep
double precision::theta_homodyne
double precision::half_kappa, half_Gamma_M
double precision::Ed
double precision::S_XX_value, S_YY_value, S_ZZ_value
double complex::S_XY_value, S_YX_value
double precision::S_homodyne_value, S_homodyne_negative, S_heterodyne_value
double precision::num_phonons_X, num_phonons_Y, num_phonons_Z
double precision::num_photons, num_phonons
double precision::area_under_spectrum
double precision::phonons_from_cooling_formula_array(3), phonons_array(3)
double precision::cooling_phonons_X, cooling_phonons_Y, cooling_phonons_Z
double precision, DIMENSION(num_points)::omega_array, S_XX_array, S_YY_array, S_ZZ_array, S_homodyne_array, S_heterodyne_array

! parameters to track spectra's peaks
integer::num_X0_values, x0_sweep
double precision::lambda_coeff
double precision::S_prev_XX_value, S_prev_YY_value, S_prev_ZZ_value, prev_omega
double precision::omega_at_peak_XX, omega_at_peak_YY, omega_at_peak_ZZ
integer::XX_max_found, YY_max_found, ZZ_max_found

integer::equation_choice

WRITE(6, *)  "Which equations to use?"
WRITE(6, *) '1 = Full 3D'
WRITE(6, *) '2 = Simplified Eq. 18'
READ(*, *)  equation_choice

! e.g. the below is:
! '4' = 4 things written to the same line
! 'ES' = scientific exponential notation (1-9E+power, instead of 0.1-0.9)
! '23' = number of positions to be used
! '15' = number of digits to the right of the decimal point
! NOTE: this is *not* the same as the precision of the number
101 FORMAT(4ES23.15)
305 FORMAT(5ES11.3)

! Here plot optical field amplitudes
OPEN(1, file="SXX.dat", status="unknown")
OPEN(2, file="SYY.dat", status="unknown")
OPEN(3, file="SZZ.dat", status="unknown")
! OPEN(4, file='OPTSRINGXX.dat', status='unknown')
! OPEN(5, file='OPTSRINGYY.dat', status='unknown')
! OPEN(7, file='OPTSRINGZZ.dat', status='unknown')

OPEN(9, file="FTSXY.dat", status="unknown")
OPEN(10, file="FTanOPT.dat", status="unknown")
OPEN(12, file="GCOUPLE.dat", status="unknown")
OPEN(14, file="PHONS.dat", status="unknown")

pi = dacos(-1.d0)
pi2 = 2.d0 * pi
detuning = detuning_Hz * pi2

! loop over the X0 equilibrium position

num_X0_values = 10
DO ii=1, num_X0_values
    ! increase x0 in increments of num_X0_values segments of half a wavelength (0.5 * lambda)
    ! wavelength = 1064 nanometres
    lambda_coeff = ii * 0.5d0 / num_X0_values
    x0_sweep = lambda_coeff * 1.064d-6
    WKX0 = Wk * x0_sweep
    theta = theta_0 * pi
    theta_homodyne = theta_homodyne_0 * pi

    WRITE(6, *) 'theta /π', theta/pi
    ! note BEAD routine uses x0 but not theta
    CALL CALCULATE_BEAD_PARAMETERS(&
        Rayleigh_range, bead_mass, polarisability, &
        E_tweezer, E_cavity, &
        kappa, Gamma_M)

    CALL CALCULATE_EQUILIBRIUM_PARAMETERS(&
        theta, WKX0, &
        G_matrix, &
        bead_mass, Rayleigh_range, polarisability, &
        E_tweezer, E_cavity, &
        detuning, kappa, Gamma_M, &
        omega_x, omega_y, omega_z, &
        phonons_from_cooling_formula_array)

    WRITE(12, 101) theta, (abs(G_matrix(ii)), nn=1, 6)

    Ed = 0.5d0 * polarisability * E_tweezer * E_cavity * sin(theta)
    ! Ed = Ed / hbar
    cavity_photons = Ed**2 * cos(WKX0)**2 / hbar**2

    WRITE(6, *) 'kappa /Hz', kappa / 2 / pi
    WRITE(6, *) 'Gamma_M /Hz', Gamma_M / 2 / pi

    half_kappa = kappa * 0.5d0
    half_Gamma_M = Gamma_M * 0.5d0

    cavity_photons = cavity_photons / (half_kappa**2 + detuning**2)
    WRITE(6, *) 'Number of photons in cavity', cavity_photons

    ! thermal bath occupancy
    num_phonons_X = kB * bath_temperature / hbar / omega_x
    num_phonons_Y = kB * bath_temperature / hbar / omega_y
    num_phonons_Z = kB * bath_temperature / hbar / omega_z
    ! shot noise
    num_photons = 0.d0

    ! open loop over frequency for noise spectra
    ! e.g. S_XX(omega) = FT(autocorrelation<X(t).X^dag(t + tau)>)
    max_omega = 2 * omega_x * 1.001
    omega_increment = max_omega / num_points
    S_homodyne_array = 0.d0

    XX_max_found = 0
    YY_max_found = 0
    ZZ_max_found = 0

    DO jj=1, num_points
        omega_sweep = -max_omega + 2 * (jj - 1) * omega_increment
        ! store frequency for integration
        omega_array(jj) = omega_sweep
        ! work out PSD of homodyne
        ! work out  for negative frequency for symmetrisation
        CALL CALCULATE_SPECTRA(&
            equation_choice, &
            N_total, &
            theta_homodyne, &
            G_matrix, &
            num_photons, &
            num_phonons_X, num_phonons_Y, num_phonons_Z, &
            detuning, kappa, Gamma_M, &
            -omega_sweep, &
            omega_x, omega_y, omega_z, &
            S_XX_value, S_YY_value, S_ZZ_value, &
            S_XY_value, S_YX_value, & ! cross-correlation spectra (can be complex)
            S_homodyne_value, S_heterodyne_value)

        ! update homodyne and symm disp.
        S_homodyne_negative = S_homodyne_negative + S_homodyne_value

        ! work out same for positive and symmetrise homodyne
        CALL CALCULATE_SPECTRA(&
            equation_choice, &
            N_total, &
            theta_homodyne, &
            G_matrix, &
            num_photons, &
            num_phonons_X, num_phonons_Y, num_phonons_Z, &
            detuning, kappa, Gamma_M, &
            omega_sweep, &
            omega_x, omega_y, omega_z, &
            S_XX_value, S_YY_value, S_ZZ_value, &
            S_XY_value, S_YX_value, &
            S_homodyne_value, S_heterodyne_value)

        ! update
        S_homodyne_negative = 0.5d0 * (S_homodyne_negative + S_homodyne_value)
        omega_sweep = omega_sweep / 2 / pi * 1.d-3

        ! apart from for the first loop, check if current PSD value is less than the previous PSD value
        ! this should detect the first point where the PSD drops
        ! so the previous loop's omega value is where the peak is
        IF (jj >= 1) THEN
            IF ((XX_max_found /= 1).AND.(S_XX_value < S_prev_XX_value)) THEN
                XX_max_found = 1
                WRITE(6, *) 'XX_now < XX_prev?', S_XX_value, S_prev_XX_value
                WRITE(6, *) 'omega where XX the peak is?', prev_omega
                omega_at_peak_XX = prev_omega
            END IF
            IF ((YY_max_found /= 1).AND.(S_YY_value < S_prev_YY_value)) THEN
                YY_max_found = 1
                WRITE(6, *) 'YY_now < YY_prev ?', S_YY_value, S_prev_YY_value
                WRITE(6, *) 'omega where YY the peak is?', prev_omega
                omega_at_peak_YY = prev_omega
            END IF
            IF ((ZZ_max_found /= 1).AND.(S_ZZ_value < S_prev_ZZ_value)) THEN
                ZZ_max_found = 1
                WRITE(6, *) 'ZZ_now < ZZ_prev?', S_ZZ_value, S_prev_ZZ_value
                WRITE(6, *) 'omega where ZZ the peak is?', prev_omega
                omega_at_peak_ZZ = prev_omega
            END IF
        END IF

        IF (equation_choice == 1) THEN
            IF ((omega_sweep.GE.120).AND.(omega_sweep.LE.170)) THEN
                WRITE(1, 101) omega_sweep, S_XX_value
                WRITE(2, 101) omega_sweep, S_YY_value
                WRITE(3, 101) omega_sweep, S_ZZ_value
                WRITE(9, 101) omega_sweep, S_XY_value, S_YX_value
                WRITE(10, 101) omega_sweep, S_heterodyne_value, S_homodyne_value
            END IF

        ELSE IF (equation_choice == 2) THEN
            WRITE(1, 101) omega_sweep, S_XX_value
            WRITE(2, 101) omega_sweep, S_YY_value
            WRITE(3, 101) omega_sweep, S_ZZ_value
            WRITE(9, 101) omega_sweep, S_XY_value, S_YX_value
            WRITE(10, 101) omega_sweep, S_heterodyne_value, S_homodyne_value

        END IF

        S_XX_array(jj) = S_XX_value
        S_YY_array(jj) = S_YY_value
        S_ZZ_array(jj) = S_ZZ_value
        S_heterodyne_array(jj) = S_heterodyne_value
        ! to find optimal squeezing
        S_homodyne_array(jj) = S_homodyne_negative

        ! keep track of calculated value per loop to compare to next loop value, to find maximum
        S_prev_XX_value = S_XX_value
        S_prev_YY_value = S_YY_value
        S_prev_ZZ_value = S_ZZ_value
        prev_omega = omega_sweep
    END DO

    WRITE(4, 101) lambda_coeff, (omega_at_peak_XX - omega_x)/pi2
    WRITE(5, 101) lambda_coeff, (omega_at_peak_YY - omega_y)/pi2
    WRITE(7, 101) lambda_coeff, (omega_at_peak_ZZ - omega_z)/pi2

    num_phonons = (omega_x + detuning)**2 + half_kappa**2
    num_phonons = num_phonons / 4. / omega_x / (-detuning)
    WRITE(6, *) 'X: back action limited phonons', num_phonons

    num_phonons = (omega_y + detuning)**2 + half_kappa**2
    num_phonons = num_phonons / 4. / omega_y / (-detuning)
    WRITE(6, *) 'Y: back action limited phonons', num_phonons

    num_phonons = (omega_z + detuning)**2 + half_kappa**2
    num_phonons = num_phonons / 4. / omega_z / (-detuning)
    WRITE(6, *) 'Z: back action limited phonons', num_phonons

    ! integrate and normalise the quantum noise spectra
    CALL INTEGRATE_SPECTRUM(num_points, omega_increment, S_XX_array, area_under_spectrum)
    ! Area corresponds to 2n+1 so convert to get n
    cooling_phonons_X = phonons_from_cooling_formula_array(1)
    phonons_array(1) = 0.5d0 * (area_under_spectrum - 1.d0)
    WRITE(6, *) 'X phonons from: cooling formula, S_XX FT'
    WRITE(6, 101) cooling_phonons_X, phonons_array(1)

    CALL INTEGRATE_SPECTRUM(num_points, omega_increment, S_YY_array, area_under_spectrum)
    cooling_phonons_Y = phonons_from_cooling_formula_array(2)
    phonons_array(2) = 0.5d0 * (area_under_spectrum - 1.d0)
    WRITE(6, *) 'Y phonons from: cooling formula, S_YY FT'
    WRITE(6, 101) cooling_phonons_Y, phonons_array(2)

    CALL INTEGRATE_SPECTRUM(num_points, omega_increment, S_ZZ_array, area_under_spectrum)
    cooling_phonons_Z = phonons_from_cooling_formula_array(3)
    phonons_array(3) = 0.5d0 * (area_under_spectrum - 1.d0)
    WRITE(6, *) 'Z phonons from: cooling formula, S_ZZ FT'
    WRITE(6, 101) cooling_phonons_Z, phonons_array(3)

    CALL INTEGRATE_SPECTRUM(num_points, omega_increment, S_heterodyne_array, area_under_spectrum)
    ! Area num_phonons corresponds to 2(nx+ny)+2 so DIFFERENT CONVERSION  to get nx+ny
    cooling_phonons_Z = phonons_from_cooling_formula_array(3)
    phonons_array(3) = 0.5d0 * (area_under_spectrum - 2.d0)
    WRITE(6, *) 'Z phonons from: cooling formula, S_heterodyne FT'
    WRITE(6, 101) cooling_phonons_Z, phonons_array(3)
    
    WRITE(14, 101) &
        detuning, &
        phonons_array(1), phonons_array(2), phonons_array(3), &
        phonons_from_cooling_formula_array(1), phonons_from_cooling_formula_array(2), phonons_from_cooling_formula_array(3)

! loop over x0
END DO

STOP
END


SUBROUTINE CALCULATE_SPECTRA(&
    equation_choice, &
    N_total, &
    theta, &
    G_matrix, &
    num_photons, &
    num_phonons_X, num_phonons_Y, num_phonons_Z, &
    detuning, kappa, Gamma_M, &
    omega, &
    omega_x, omega_y, omega_z, &
    S_XX_value, S_YY_value, S_ZZ_value, &
    S_XY_value, S_YX_value, &
    S_homodyne_value, S_heterodyne_value)
    ! """
    !  Generic routine for noise spectra of trap and probe beams
    ! """
    IMPLICIT NONE

    integer::equation_choice
    integer::N_total
    double precision::theta
    double precision::G_matrix(6)
    double precision::num_photons
    double precision::num_phonons_X, num_phonons_Y, num_phonons_Z
    double precision::detuning, kappa, Gamma_M
    double precision::omega
    double precision::omega_x, omega_y, omega_z
    double precision::S_XX_value, S_YY_value, S_ZZ_value
    double complex::S_XY_value, S_YX_value
    double precision::S_homodyne_value, S_heterodyne_value

    integer::ii
    double precision::pi, pi2
    double precision::G_average
    double precision::omega_heterodyne, omega_addition
    double complex::GXY
    double complex::chi_C, chi_C_star
    double complex::chi_X, chi_X_star
    double complex::chi_Y, chi_Y_star
    double complex::chi_Z, chi_Z_star
    double complex::a_hat_vector(1, N_total), a_dag_vector(1, N_total)
    double complex::b_x_vector(1, N_total), b_y_vector(1, N_total), b_z_vector(1, N_total)

    pi = dacos(-1.d0)

    S_homodyne_value = 0.d0
    b_x_vector = 0.d0
    b_y_vector = 0.d0
    b_z_vector = 0.d0

    ! WORK OUT NOISE SPECTRA
    ! First do susceptibilities
    CALL CALCULATE_SUSCEPTIBILITIES(&
        omega, &
        omega_x, omega_y, omega_z, &
        detuning, kappa, Gamma_M, &
        chi_C, chi_C_star, &
        chi_X, chi_X_star, &
        chi_Y, chi_Y_star, &
        chi_Z, chi_Z_star)

    ! work out noise vector for x,  a_hat_vector and a2
    CALL CALCULATE_NOISE_VECTORS(&
        equation_choice, &
        N_total, &
        G_matrix, &
        kappa, Gamma_M, &
        chi_C, chi_C_star, &
        chi_X, chi_X_star, &
        chi_Y, chi_Y_star, &
        chi_Z, chi_Z_star, &
        a_hat_vector, a_dag_vector, &
        b_x_vector, b_y_vector, b_z_vector)

    S_XX_value = abs(b_x_vector(1, 1))**2
    S_XX_value = S_XX_value + (num_phonons_X + 1) * abs(b_x_vector(1, 3))**2 + num_phonons_X * abs(b_x_vector(1, 4))**2
    S_XX_value = S_XX_value + (num_phonons_Y + 1) * abs(b_x_vector(1, 5))**2 + num_phonons_Y * abs(b_x_vector(1, 6))**2
    S_XX_value = S_XX_value + (num_phonons_Z + 1) * abs(b_x_vector(1, 7))**2 + num_phonons_Z * abs(b_x_vector(1, 8))**2

    S_YY_value = abs(b_y_vector(1, 1))**2
    S_YY_value = S_YY_value + (num_phonons_X + 1) * abs(b_y_vector(1, 3))**2 + num_phonons_X * abs(b_y_vector(1, 4))**2
    S_YY_value = S_YY_value + (num_phonons_Y + 1) * abs(b_y_vector(1, 5))**2 + num_phonons_Y * abs(b_y_vector(1, 6))**2
    S_YY_value = S_YY_value + (num_phonons_Z + 1) * abs(b_y_vector(1, 7))**2 + num_phonons_Z * abs(b_y_vector(1, 8))**2

    S_ZZ_value = abs(b_z_vector(1, 1))**2
    S_ZZ_value = S_ZZ_value + (num_phonons_X + 1) * abs(b_z_vector(1, 3))**2 + num_phonons_X * abs(b_z_vector(1, 4))**2
    S_ZZ_value = S_ZZ_value + (num_phonons_Y + 1) * abs(b_z_vector(1, 5))**2 + num_phonons_Y * abs(b_z_vector(1, 6))**2
    S_ZZ_value = S_ZZ_value + (num_phonons_Z + 1) * abs(b_z_vector(1, 7))**2 + num_phonons_Z * abs(b_z_vector(1, 8))**2

    ! WORK OUT CROSS SPECTRA PSD S_XY and S_YX
    S_XY_value = b_x_vector(1, 1) * conjg(b_y_vector(1, 1))
    S_XY_value = S_XY_value + (num_phonons_X + 1) * b_x_vector(1, 3) * conjg(b_y_vector(1, 3))
    S_XY_value = S_XY_value + num_phonons_X * b_x_vector(1, 4) * conjg(b_y_vector(1, 4))
    S_XY_value = S_XY_value + (num_phonons_Y + 1) * b_x_vector(1, 5) * conjg(b_y_vector(1, 5))
    S_XY_value = S_XY_value + num_phonons_Y * b_x_vector(1, 6) * conjg(b_y_vector(1, 6))
    S_XY_value = S_XY_value + (num_phonons_Z + 1) * b_x_vector(1, 7) * conjg(b_y_vector(1, 7))
    S_XY_value = S_XY_value + num_phonons_Z * b_x_vector(1, 8) * conjg(b_y_vector(1, 8))

    S_YX_value = b_y_vector(1, 1) * conjg(b_x_vector(1, 1))
    S_YX_value = S_YX_value + (num_phonons_X + 1) * b_y_vector(1, 3) * conjg(b_x_vector(1, 3))
    S_YX_value = S_YX_value + num_phonons_X * b_y_vector(1, 4) * conjg(b_x_vector(1, 4))
    S_YX_value = S_YX_value + (num_phonons_Y + 1) * b_y_vector(1, 5) * conjg(b_x_vector(1, 5))
    S_YX_value = S_YX_value + num_phonons_Y * b_y_vector(1, 6) * conjg(b_x_vector(1, 6))
    S_YX_value = S_YX_value + (num_phonons_Z + 1) * b_y_vector(1, 7) * conjg(b_x_vector(1, 7))
    S_YX_value = S_YX_value + num_phonons_Z * b_y_vector(1, 8) * conjg(b_x_vector(1, 8))

    S_XY_value = 0.5 * (S_XY_value + S_YX_value)

    IF (equation_choice == 2) THEN
        GXY = a_hat_vector(1, 4)
        ! IN THIS VERSION WE WORK OUT Eq. 18,  simplified version of cross-correlation
        S_YX_value = (S_XX_value - S_YY_value) * GXY / (omega_y - omega_x)
    END IF

    ! work out homodyne spectra using same vectors
    CALL CALCULATE_HOMODYNE_FROM_OPTICAL_NOISE_VECTORS(&
        N_total, &
        theta, &
        num_phonons_X, num_phonons_Y, num_phonons_Z, &
        num_photons, &
        a_hat_vector, a_dag_vector, &
        S_homodyne_value)
    ! rescaling of heterodyne for Gx / Gy
    G_average = (G_matrix(1) + G_matrix(2)) / 2.

    S_homodyne_value = (S_homodyne_value - 1.d0) / abs(chi_C - chi_C_star)**2 / G_average**2 / kappa / 4.

    ! work out heterodyne for positive frequency branch.
    ! Shift frequency by heterodyne beat.
    omega_addition = omega_x * 4.d0
    omega_addition = 0.d0
    omega_heterodyne = omega + omega_addition

    ! work out susceptibilities shifted in frequency
    CALL CALCULATE_SUSCEPTIBILITIES(&
        omega_heterodyne, &
        omega_x, omega_y, omega_z, &
        detuning, kappa, Gamma_M, &
        chi_C, chi_C_star, &
        chi_X, chi_X_star, &
        chi_Y, chi_Y_star, &
        chi_Z, chi_Z_star)

    ! work out noise vector again
    CALL CALCULATE_NOISE_VECTORS(&
        equation_choice, &
        N_total, &
        G_matrix, &
        kappa, Gamma_M, &
        chi_C, chi_C_star, &
        chi_X, chi_X_star, &
        chi_Y, chi_Y_star, &
        chi_Z, chi_Z_star, &
        a_hat_vector, a_dag_vector, &
        b_x_vector, b_y_vector, b_z_vector)

    CALL CALCULATE_HETERODYNE_FROM_OPTICAL_NOISE_VECTORS(&
        N_total, &
        theta, &
        num_photons, &
        num_phonons_X, num_phonons_Y, num_phonons_Z, &
        a_hat_vector, a_dag_vector, &
        S_heterodyne_value)

    ! WRITE(6, *) omega_heterodyne, S_heterodyne_value, (S_heterodyne_value) / abs(chi_C)**2 / G_matrix(1)**2
    S_heterodyne_value = (S_heterodyne_value - 1.d0) / abs(chi_C_star)**2 / G_average**2 / kappa / 4.
    ! S_heterodyne_value = S_heterodyne_value - 1.d0

    RETURN
END


SUBROUTINE CALCULATE_SUSCEPTIBILITIES(&
    omega, &
    omega_x, omega_y, omega_z, &
    detuning, kappa, Gamma_M, &
    chi_C, chi_C_star, &
    chi_X, chi_X_star, &
    chi_Y, chi_Y_star, &
    chi_Z, chi_Z_star)
    ! """
    ! work out susceptibilities for noise spectra
    ! """
    IMPLICIT NONE

    double precision::omega
    double precision::omega_x, omega_y, omega_z
    double precision::detuning, kappa, Gamma_M
    double complex::chi_C, chi_C_star
    double complex::chi_X, chi_X_star
    double complex::chi_Y, chi_Y_star
    double complex::chi_Z, chi_Z_star

    double precision::chi_real, chi_imag, chi_denominator
    double precision::half_kappa, half_Gamma_M

    half_kappa = kappa * 0.5d0
    half_Gamma_M = Gamma_M * 0.5d0

    ! optical susceptibilities
    ! chi_C
    chi_denominator = half_kappa**2 + (omega + detuning)**2
    chi_real = half_kappa/chi_denominator
    chi_imag = (omega + detuning)/chi_denominator
    chi_C = cmplx(chi_real, chi_imag)
    ! chi_C^*(-omega) - i.e. turn omega into -omega, then chi_imag into -chi_imag.
    chi_denominator = half_kappa**2 + (-omega + detuning)**2
    chi_real = half_kappa/chi_denominator
    chi_imag = (-omega + detuning)/chi_denominator
    chi_C_star = cmplx(chi_real, -chi_imag)

    ! mechanical susceptibilities
    ! chi_X
    chi_denominator = half_Gamma_M**2 + (omega - omega_x)**2
    chi_real =  half_Gamma_M/chi_denominator
    chi_imag = (omega - omega_x)/chi_denominator
    chi_X = cmplx(chi_real, chi_imag)
    ! chi_X*(-omega)
    chi_denominator = half_Gamma_M**2 + (-omega - omega_x)**2
    chi_real =  half_Gamma_M/chi_denominator
    chi_imag = (-omega - omega_x)/chi_denominator
    chi_X_star = cmplx(chi_real, -chi_imag)

    ! chi_Y
    chi_denominator = half_Gamma_M**2 + (omega - omega_y)**2
    chi_real =  half_Gamma_M/chi_denominator
    chi_imag = (omega - omega_y)/chi_denominator
    chi_Y = cmplx(chi_real, chi_imag)
    ! chi_Y*(-omega)
    chi_denominator = half_Gamma_M**2 + (-omega - omega_y)**2
    chi_real =  half_Gamma_M/chi_denominator
    chi_imag = (-omega - omega_y)/chi_denominator
    chi_Y_star = cmplx(chi_real, -chi_imag)

    ! chi_Z
    chi_denominator = half_Gamma_M**2 + (omega - omega_z)**2
    chi_real =  half_Gamma_M/chi_denominator
    chi_imag = (omega - omega_z)/chi_denominator
    chi_Z = cmplx(chi_real, chi_imag)
    ! chi_Z*(-omega)
    chi_denominator = half_Gamma_M**2 + (-omega - omega_z)**2
    chi_real =  half_Gamma_M/chi_denominator
    chi_imag = (-omega - omega_z)/chi_denominator
    chi_Z_star = cmplx(chi_real, -chi_imag)

    RETURN
END


SUBROUTINE CALCULATE_NOISE_VECTORS(&
    equation_choice, &
    N_total, &
    G_matrix, &
    kappa, Gamma_M, &
    chi_C, chi_C_star, &
    chi_X, chi_X_star, &
    chi_Y, chi_Y_star, &
    chi_Z, chi_Z_star, &
    a_hat_vector, a_dag_vector, &
    b_x_vector, b_y_vector, b_z_vector)
    ! """
    ! Work our vectors of X, Y, Z, and optical fields, in terms of noise operators
    ! """
    IMPLICIT NONE

    integer::equation_choice
    integer::N_total
    double precision::G_matrix(6)
    double precision::kappa, Gamma_M
    double complex::chi_C, chi_C_star
    double complex::chi_X, chi_X_star
    double complex::chi_Y, chi_Y_star
    double complex::chi_Z, chi_Z_star
    double complex, INTENT(OUT)::a_hat_vector(1, N_total), a_dag_vector(1, N_total)
    double complex, INTENT(OUT)::b_x_vector(1, N_total), b_y_vector(1, N_total), b_z_vector(1, N_total)

    integer::ii
    double precision::GX, GY, GZ
    double precision::GXY, GZX, GYZ
    double complex::imag_num, one
    double complex::beta_x, beta_y, beta_z
    double complex::G_coef_XY, G_coef_ZX, G_coef_YZ
    double complex::G_coef_YX, G_coef_XZ, G_coef_ZY
    double complex::RXY, RZX, RYZ
    double complex::RYX, RXZ, RZY
    double complex::M_X, M_Y, M_Z
    double complex::CSUM, CNORM
    double complex::mu_X, mu_Y, mu_Z
    double complex::eta_C_0, eta_C_neg_half_pi, eta_C_pos_half_pi
    double complex::N0X(N_total), N0Y(N_total), N0Z(N_total)

    GX = G_matrix(1)
    GY = G_matrix(2)
    GZ = G_matrix(3)
    GXY = G_matrix(4)
    GYZ = G_matrix(5)
    GZX = G_matrix(6)

    a_hat_vector = cmplx(0.d0, 0.d0)
    a_dag_vector = cmplx(0.d0, 0.d0)

    b_x_vector = cmplx(0.d0, 0.d0)
    b_y_vector = cmplx(0.d0, 0.d0)
    b_z_vector = cmplx(0.d0, 0.d0)

    N0X = cmplx(0.d0, 0.d0)
    N0Y = cmplx(0.d0, 0.d0)
    N0Z = cmplx(0.d0, 0.d0)

    imag_num = cmplx(0.d0, 1.d0)
    one = cmplx(1.d0, 0.0d0)

    ! SUSCEPTIBILITIES
    ! eta_C_phi = e^(-i.phi) * chi_C + e^(i.phi) * chi_C_star
    eta_C_0 = chi_C - chi_C_star
    eta_C_neg_half_pi = imag_num * (chi_C + chi_C_star)
    eta_C_pos_half_pi = -imag_num * (chi_C + chi_C_star)
    mu_X = chi_X - chi_X_star
    mu_Y = chi_Y - chi_Y_star
    mu_Z = chi_Z - chi_Z_star

    ! coeff of X-Y coupling- Combines direct and indirect paths
    G_coef_XY = GXY + imag_num * eta_C_0 * GX * GY
    G_coef_YX = GXY + imag_num * eta_C_0 * GX * GY
    G_coef_XZ = GZX + imag_num * eta_C_neg_half_pi * GZ * GX
    G_coef_ZX = GZX + imag_num * eta_C_pos_half_pi * GZ * GX
    G_coef_YZ = GYZ + imag_num * eta_C_neg_half_pi * GY * GZ
    G_coef_ZY = GYZ + imag_num * eta_C_pos_half_pi * GY * GZ

    ! NORMALIZATIONS
    M_X = 1.d0 + GX**2 * mu_X * eta_C_0
    M_Y = 1.d0 + GY**2 * mu_Y * eta_C_0
    M_Z = 1.d0 + GZ**2 * mu_Z * eta_C_0

    ! this is the coefficient of the amplitude (phase) quadrature's noise vector for x,y (z)
    beta_x = imag_num * sqrt(kappa) * mu_X * GX
    beta_y = imag_num * sqrt(kappa) * mu_Y * GY
    beta_z = imag_num * sqrt(kappa) * mu_Z * GZ

    ! these are the 1D positional noise vectors, i.e. qj_1D
    ! zero-th order X noise vector; weights of a_hat_vector, a_hat_vector*, bx, bx*, by, by*, bz, bz*
    N0X(1) = beta_x * chi_C / M_X
    N0X(2) = beta_x * chi_C_star / M_X
    N0X(3) = sqrt(Gamma_M) * chi_X / M_X
    N0X(4) = sqrt(Gamma_M) * chi_X_star / M_X

    ! zero-th order Y noise vector; weights of a_hat_vector, a_hat_vector*, bx, bx*, by, by*, bz, bz*
    N0Y(1) = beta_y * chi_C / M_Y
    N0Y(2) = beta_y * chi_C_star / M_Y
    N0Y(5) = sqrt(Gamma_M) * chi_Y / M_Y
    N0Y(6) = sqrt(Gamma_M) * chi_Y_star / M_Y

    ! zero-th  order Z noise vector; weights of a_hat_vector, a_hat_vector*, bx, bx*, by, by*, bz, bz*
    N0Z(1) = -imag_num * beta_z * chi_C / M_Z
    N0Z(2) = imag_num * beta_z * chi_C_star / M_Z
    N0Z(7) = sqrt(Gamma_M) * chi_Z / M_Z
    N0Z(8) = sqrt(Gamma_M) * chi_Z_star / M_Z

    ! r_3D = (x_3D, y_3D, z_3D) - column vector of position operators for each axis considered in 3D
    ! r_1D = (x_1D, y_1D, z_1D) - column vector of position operators for each axis considered in 1D

    ! qj_1D = M_j^-1 * (sqrt(Gamma_M).qj_thermal + i.mu_j.sqrt(kappa)*(gjY.Y_in + gjP.P_in)
    ! NOTE: x,y coupled to Y while z couples to P
    ! i.e. GX = gxY, GY = gyY, GZ = gzP

    ! thermal noise for each axis
    ! qj_thermal = chi_j.b_in_j + chi_j_star.b_in_j_dagger

    ! optical noises
    ! Y_in = chi_C.a_in + chi_C_star.a_in_dagger
    ! P_in = i*(chi_C_star.a_in_dagger - chi_C.a_in)

    ! Higher order, entries in matrix R from matrix equation r_3D = R.r_3D + r_1D
    RXY = imag_num * mu_X * G_coef_XY / M_X
    RYX = imag_num * mu_Y * G_coef_YX / M_Y
    RXZ = imag_num * mu_X * G_coef_XZ / M_X
    RZX = imag_num * mu_Z * G_coef_ZX / M_Z
    RYZ = imag_num * mu_Y * G_coef_YZ / M_Y
    RZY = imag_num * mu_Z * G_coef_ZY / M_Z

    ! determinant of final matrix A = 1 - R from matrix equation r_3D = A.r_1D
    CNORM = 1.d0 - ((RZX * RXZ) + (RZY * RYZ) + (RYX * RXY)) - ((RZX * RXY * RYZ) + (RYX * RXZ * RZY))

    ! ADD 3D BACK-ACTION TERMS
    DO ii=1, N_total
        ! 1st (X) column of adjoint of A
        CSUM = (1.d0 - RZY * RYZ) * N0X(ii) + (RXY + RXZ * RZY) * N0Y(ii) + (RXZ + RXY * RYZ) * N0Z(ii)
        b_x_vector(1, ii) = b_x_vector(1, ii) + CSUM / CNORM

        ! 2nd (Y) column of adjoint of A
        CSUM = (1.d0 - RZX * RXZ) * N0Y(ii) + (RYX + RYZ * RZX) * N0X(ii) + (RYZ + RYX * RXZ) * N0Z(ii)
        b_y_vector(1, ii) = b_y_vector(1, ii) + CSUM / CNORM

        ! 3rd (Z) column of adjoint of A
        CSUM = (1.d0 - RYX * RXY) * N0Z(ii) + (RZX + RZY * RYX) * N0X(ii) + (RZY + RZX * RXY) * N0Y(ii)
        b_z_vector(1, ii) = b_z_vector(1, ii) + CSUM / CNORM
    END DO

    DO ii=1, N_total
        a_hat_vector(1, ii) = imag_num * chi_C * (GX * b_x_vector(1, ii) + GY * b_y_vector(1, ii))
        a_hat_vector(1, ii) = a_hat_vector(1, ii) - GZ * chi_C * b_z_vector(1, ii)
        a_dag_vector(1, ii) = -imag_num * chi_C_star * (GX * b_x_vector(1, ii) + GY * b_y_vector(1, ii))
        a_dag_vector(1, ii) = a_dag_vector(1, ii) - GZ * chi_C_star * b_z_vector(1, ii)
    END DO

    ! add shot or incoming noise
    ! trap beam: add cavity-filtered contribution
    a_hat_vector(1, 1) = a_hat_vector(1, 1) + sqrt(kappa) * chi_C
    a_dag_vector(1, 2) = a_dag_vector(1, 2) + sqrt(kappa) * chi_C_star

    ! cavity output : add incoming imprecision
    ! work out a_out = a_in - sqrt(kappa) * a
    DO ii=1, N_total
        a_hat_vector(1, ii) = -a_hat_vector(1, ii) * sqrt(kappa)
        a_dag_vector(1, ii) = -a_dag_vector(1, ii) * sqrt(kappa)
    END DO

    a_hat_vector(1, 1) = one + a_hat_vector(1, 1)
    a_dag_vector(1, 2) = one + a_dag_vector(1, 2)

    IF (equation_choice == 2) THEN
        a_hat_vector(1, 1) = RXY
        a_hat_vector(1, 2) = RYX
        a_hat_vector(1, 3) = G_coef_XY
        a_hat_vector(1, 4) = G_coef_YX
        a_hat_vector(1, 5) = mu_X / M_X
        a_hat_vector(1, 6) = mu_Y / M_Y
    END IF

    RETURN
END


SUBROUTINE INTEGRATE_SPECTRUM(num_points, omega_increment, spectrum_array, area_under_spectrum)
    ! """
    ! Integrates a given array of spectrum values using the trapezium rule
    ! """
    IMPLICIT NONE

    integer::num_points
    double precision::omega_increment
    double precision::spectrum_array(num_points)
    double precision::area_under_spectrum

    integer::ii
    double precision::pi2, trapeze_base

    pi2 = 2.d0 * dacos(-1.d0)

    ! integrate the position spectrum of bead
    ! quick hack - use trapezoidal rule - improve later?
    ! each trapeze base is the omega_increment i.e. omega_range/num_points 

    ! NOTE: divide by 2π to convert omega to non-angular frequency (Hz)
    ! TODO: explain the division by 2π as a result of Parseval's theorem / "Sum-rule" ?
    ! this means the area is the area under the spectrum curve plotted vs. Hz
    trapeze_base = omega_increment/pi2

    area_under_spectrum = 0.d0
    DO ii=1, num_points - 1
        area_under_spectrum = area_under_spectrum + 0.5d0 * (spectrum_array(ii) + spectrum_array(ii + 1))
    END DO
    area_under_spectrum = area_under_spectrum * trapeze_base

    RETURN
END


SUBROUTINE CALCULATE_HOMODYNE_FROM_OPTICAL_NOISE_VECTORS(&
    N_total, &
    theta, &
    num_photons, &
    num_phonons_X, num_phonons_Y, num_phonons_Z, &
    a_hat_vector, a_dag_vector, &
    S_homodyne_value)
    ! """
    ! Calculates the homodyne spectrum from a_hat and a_dagger
    ! """
    IMPLICIT NONE

    integer::N_total
    double precision::theta
    double precision::num_photons
    double precision::num_phonons_X, num_phonons_Y, num_phonons_Z
    double complex::a_hat_vector(1, N_total), a_dag_vector(1, N_total)
    double precision::S_homodyne_value

    integer::ii
    double complex::imag_num
    double complex::PY_hat_vector(1, N_total), Y_hat_vector(1, N_total)
    double complex::homodyne_noise_vector(1, N_total)

    PY_hat_vector = (0.d0, 0.d0)
    Y_hat_vector = (0.d0, 0.d0)
    homodyne_noise_vector = (0.d0, 0.d0)

    imag_num = cmplx(0.d0, 1.d0)
    DO ii=1, N_total
        Y_hat_vector(1, ii) = a_hat_vector(1, ii) + a_dag_vector(1, ii)
        PY_hat_vector(1, ii) = imag_num * (a_hat_vector(1, ii) - a_dag_vector(1, ii))
        homodyne_noise_vector(1, ii) = PY_hat_vector(1, ii) * sin(theta) + Y_hat_vector(1, ii) * cos(theta)
    END DO

    S_homodyne_value = abs(homodyne_noise_vector(1, 1))**2 * num_photons
    S_homodyne_value = S_homodyne_value + abs(homodyne_noise_vector(1, 2))**2 * (num_photons + 1.d0)
    S_homodyne_value = S_homodyne_value + abs(homodyne_noise_vector(1, 3))**2 * num_phonons_X
    S_homodyne_value = S_homodyne_value + abs(homodyne_noise_vector(1, 4))**2 * (num_phonons_X + 1.d0)
    S_homodyne_value = S_homodyne_value + abs(homodyne_noise_vector(1, 5))**2 * num_phonons_Y
    S_homodyne_value = S_homodyne_value + abs(homodyne_noise_vector(1, 6))**2 * (num_phonons_Y + 1.d0)
    S_homodyne_value = S_homodyne_value + abs(homodyne_noise_vector(1, 7))**2 * num_phonons_Z
    S_homodyne_value = S_homodyne_value + abs(homodyne_noise_vector(1, 8))**2 * (num_phonons_Z + 1.d0)

    RETURN
END


SUBROUTINE CALCULATE_HETERODYNE_FROM_OPTICAL_NOISE_VECTORS(&
    N_total, &
    theta, &
    num_photons, &
    num_phonons_X, num_phonons_Y, num_phonons_Z, &
    a_hat_vector, a_dag_vector, &
    S_heterodyne_value)
    ! """
    ! Calculates the heterodyne spectrum from a_hat and a_dagger
    ! """
    IMPLICIT NONE

    integer::N_total
    double precision::theta
    double precision::num_photons
    double precision::num_phonons_X, num_phonons_Y, num_phonons_Z
    double complex::a_hat_vector(1, N_total), a_dag_vector(1, N_total)
    double precision::S_heterodyne_value

    integer::ii
    double complex::imag_num
    double complex::PY_hat_vector(1, N_total), Y_hat_vector(1, N_total)
    double complex::heterodyne_noise_vector(1, N_total)

    ! PY_hat_vector=(0.d0, 0.d0)
    ! Y_hat_vector=(0.d0, 0.d0)
    heterodyne_noise_vector = (0.d0, 0.d0)

    imag_num = cmplx(0.d0, 1.d0)
    DO ii=1, N_total
        ! Y_hat_vector(1, ii) = a_hat_vector(1, ii) + a_dag_vector(1, ii)
        ! PY_hat_vector(1, ii) = imag_num * (a_hat_vector(1, ii) - a_dag_vector(1, ii))
        ! heterodyne_noise_vector(1, ii) = PY_hat_vector(1, ii) * sin(theta) + Y_hat_vector(1, ii) * cos(theta)
        heterodyne_noise_vector(1, ii) = a_dag_vector(1, ii)
    END DO

    S_heterodyne_value = abs(heterodyne_noise_vector(1, 1))**2 * num_photons
    S_heterodyne_value = S_heterodyne_value + abs(heterodyne_noise_vector(1, 2))**2 * (num_photons + 1.d0)
    S_heterodyne_value = S_heterodyne_value + abs(heterodyne_noise_vector(1, 3))**2 * num_phonons_X
    S_heterodyne_value = S_heterodyne_value + abs(heterodyne_noise_vector(1, 4))**2 * (num_phonons_X + 1.d0)
    S_heterodyne_value = S_heterodyne_value + abs(heterodyne_noise_vector(1, 5))**2 * num_phonons_Y
    S_heterodyne_value = S_heterodyne_value + abs(heterodyne_noise_vector(1, 6))**2 * (num_phonons_Y + 1.d0)
    S_heterodyne_value = S_heterodyne_value + abs(heterodyne_noise_vector(1, 7))**2 * num_phonons_Z
    S_heterodyne_value = S_heterodyne_value + abs(heterodyne_noise_vector(1, 8))**2 * (num_phonons_Z + 1.d0)

    RETURN
END


SUBROUTINE CALCULATE_BEAD_PARAMETERS(&
    Rayleigh_range, bead_mass, polarisability, &
    E_tweezer, E_cavity, &
    kappa, Gamma_M)
    ! """
    ! subroutine below is provided by user
    ! and calculates  relevant parameters
    ! """
    IMPLICIT NONE

    integer::N_period, N_total
    double precision::bead_radius, bead_density
    double precision::epsilon_0, epsilon_R
    double precision::c, hbar, kB, bath_temperature, Gravity
    double precision::WK, WX, WY
    double precision::waist_radius, cavity_length, Finesse, air_pressure
    double precision::tweezer_input_power, detuning_Hz, DelFSR, theta_0
    double precision::X0, Y0, Z0
    INCLUDE 'CSCAVITY.h'

    double precision::Rayleigh_range, bead_mass, polarisability
    double precision::E_tweezer, E_cavity
    double precision::kappa, Gamma_M

    double precision::pi
    double precision::omega_optical, coeff, cavity_volume, amplitude
    double precision::kappa_in, kappa_nano

    ! zero eq. initial values
    pi = dacos(-1.d0)

    Rayleigh_range = WX * WY * WK / 2.d0
    WRITE(6, *) 'Rayleigh_range /m = ', Rayleigh_range

    bead_mass = bead_density * 4.*pi/3. * bead_radius**3
    WRITE(6, *) 'bead_mass /kg = ', bead_mass

    polarisability = 4.* pi * epsilon_R * (epsilon_0 - 1.) / (epsilon_0 + 2.) * bead_radius**3
    WRITE(6, *) 'bead polarisability / Farad.m^2 = ', polarisability

    omega_optical = c * WK
    ! add a factor of 4 here. Not in GALAX1-5 routines!!!
    cavity_volume = cavity_length * pi * waist_radius**2 / 4.d0

    ! Depth of cavity field. Weak and unimportant for CS case
    amplitude = omega_optical * polarisability / 2. / cavity_volume / epsilon_R
    WRITE(6, *) 'Cavity trap, amplitude/2π (zeroth shift) /Hz = '
    WRITE(6, *) amplitude / pi / 2.,  amplitude * cos(WK * x0)**2

    E_tweezer = sqrt(4. * tweezer_input_power / (Wx * Wy * pi * c * epsilon_R))
    WRITE(6, *) 'E_tweezer /Volts.m^-1 = ', E_tweezer

    E_cavity = sqrt(hbar * omega_optical / (2. * cavity_volume * epsilon_0))
    WRITE(6, *) 'E_cavity /Volts.m^-1 = ', E_cavity

    kappa_in = pi * c / Finesse / cavity_length
    WRITE(6, *) 'kappa_in /Hz = ', kappa_in / 2 / pi

    coeff = WK * polarisability / epsilon_R / omega_optical**2 / pi
    kappa_nano = 4.* coeff**2 * DelFSR * cos(WK * x0)**2
    WRITE(6, *) 'kappa_nano /Hz = ', kappa_nano

    kappa = kappa_nano + kappa_in
    WRITE(6, *) 'Optical damping: kappa /Hz = ', kappa

    ! take usual expression e.g. Levitated review by Li Geraci etc
    ! 1 bar = 10^5 pascal; air_pressure /mbar = 10^2 Pascal
    ! Gamma_M = 1600 * air_pressure / (pi * air_speed * bead_density * bead_radius)
    ! air_speed = 500 ms^-1
    Gamma_M = 1600. * air_pressure / pi / 500. / bead_density / bead_radius
    WRITE(6, *) 'Mechanical damping: Gamma_M /Hz = ', Gamma_M

    RETURN
END


SUBROUTINE CALCULATE_EQUILIBRIUM_PARAMETERS(&
    theta, WKX0, &
    G_matrix, &
    Rayleigh_range, bead_mass, polarisability, &
    E_tweezer, E_cavity, &
    detuning, kappa, Gamma_M, &
    omega_x, omega_y, omega_z, &
    phonons_from_cooling_formula_array)
    ! """
    ! subroutine below obtains the optomechanical parameters
    ! """
    IMPLICIT NONE

    integer::N_period, N_total
    double precision::bead_radius, bead_density
    double precision::epsilon_0, epsilon_R
    double precision::c, hbar, kB, bath_temperature, Gravity
    double precision::WK, WX, WY
    double precision::waist_radius, cavity_length, Finesse, air_pressure
    double precision::tweezer_input_power, detuning_Hz, DelFSR, theta_0
    double precision::X0, Y0, Z0
    INCLUDE 'CSCAVITY.h'

    double precision::theta, WKX0
    double precision::G_matrix(6)
    double precision::Rayleigh_range, bead_mass, polarisability
    double precision::E_tweezer, E_cavity
    double precision::detuning, kappa, Gamma_M
    double precision::omega_x, omega_y, omega_z
    ! analytic equilibrium phonon numbers
    double precision::phonons_from_cooling_formula_array(3)

    double precision::pi
    double precision::N_photon
    double precision::alpha_real, alpha_imag
    double precision::C1, C2
    double precision::GX, GY, GZ
    double precision::GXY, GZX, GYZ
    double precision::X_zpf, Y_zpf, Z_zpf
    double precision::Edip, Ediph
    double precision::half_kappa, kappa_nano, kappa_in
    double precision::cooling_rate_x, cooling_rate_y, cooling_rate_z

    pi = dacos(-1.d0)
    half_kappa = 0.5d0 * kappa
    WRITE(6, *) 'kappa  /kHz =', kappa / 2 / pi

    ! note that the detunings include zeroth order correction for linearised versions
    ! as a first pass work out frequencies with Wk0 = 0
    WKX0 = WK * X0

    omega_x = polarisability * E_tweezer**2 / bead_mass / WX**2
    omega_y = polarisability * E_tweezer**2 / bead_mass / Wy**2
    omega_z = 0.5d0 * polarisability * E_tweezer**2 / bead_mass / Rayleigh_range**2

    WRITE(6, *) 'WKX0 /π = ', WKX0 / pi

    ! optomechanical drive frequency
    WRITE(6, *) 'E_tweezer, E_cavity'
    WRITE(6, 100)  E_tweezer, E_cavity

    ! Sept 5 we will use negative Edip
    Edip = -0.5d0 * polarisability * E_tweezer * E_cavity * sin(theta)
    Ediph = Edip / hbar
    WRITE(6, *) 'Edip / 2π / hbar = '
    WRITE(6, 100) Ediph / 2 / pi

    ! photon number in cavity
    ! real part of photon field
    alpha_real = detuning * Ediph * cos(WKX0) / (half_kappa**2 + detuning**2)
    alpha_imag = -half_kappa * Ediph * cos(WKX0) / (half_kappa**2 + detuning**2)

    N_photon = Ediph * Ediph * cos(WKX0) * cos(WKX0)
    N_photon = N_photon / (half_kappa**2 + detuning**2)
    WRITE(6, *) 'detuning /Hz = ', detuning / 2 / pi
    WRITE(6, *) 'half_kappa /Hz = ', half_kappa / 2 / pi
    WRITE(6, *) 'Number of photons in cavity = ', N_photon

    ! ADD CS POTENTIAL CORRECTION to frequency squared
    C1 = -Edip / bead_mass * 2. * alpha_real * WK**2 * cos(WKX0)
    omega_x = omega_x + C1 * sin(theta) * sin(theta)
    omega_y = omega_y + C1 * cos(theta) * cos(theta)
    omega_z = omega_z - 2. * Edip / bead_mass * alpha_real * (WK - 1.d0 / Rayleigh_range)**2 * cos(WKX0)

    omega_x = sqrt(omega_x)
    omega_y = sqrt(omega_y)
    omega_z = sqrt(omega_z)
    WRITE(6, *) 'mech freq omega_x /Hz = ', omega_x / 2 / pi
    WRITE(6, *) 'mech freq omega_y /Hz = ', omega_y / 2 / pi
    WRITE(6, *) 'mech freq omega_z /Hz = ', omega_z / 2 / pi

    ! optomechanical couplings
    X_zpf = sqrt(hbar / (2.d0 * bead_mass * omega_x))
    Y_zpf = sqrt(hbar / (2.d0 * bead_mass * omega_y))
    Z_zpf = sqrt(hbar / (2.d0 * bead_mass * omega_z))

    GX = Ediph * WK * X_zpf * sin(theta) * sin(WKX0)
    GY = Ediph * WK * Y_zpf * cos(theta) * sin(WKX0)
    GZ = -Ediph * (WK - 1.d0 / Rayleigh_range) * Z_zpf * cos(WKX0)
    WRITE(6, *) 'GX /Hz', GX / 2 / pi
    WRITE(6, *) 'GY /Hz', GY / 2 / pi
    WRITE(6, *) 'GZ /Hz', GZ / 2 / pi

    GXY = Ediph * WK * X_zpf * WK * Y_zpf * alpha_real * sin(2. * theta) * cos(WKX0)
    GZX = 2.d0 * Ediph * (WK - 1.d0 / Rayleigh_range) * Z_zpf * WK * X_zpf * alpha_imag * sin(WKX0) * sin(theta)
    GYZ = 2.d0 * Ediph * (WK - 1.d0 / Rayleigh_range) * Z_zpf * WK * Y_zpf * alpha_imag * sin(WKX0) * cos(theta)
    WRITE (*, *) 'GXY /Hz', GXY / 2 / pi
    WRITE (*, *) 'GYZ /Hz', GYZ / 2 / pi
    WRITE (*, *) 'GZX /Hz', GZX / 2 / pi

    ! assign these couplings to an array
    G_matrix(1) = GX
    G_matrix(2) = GY
    G_matrix(3) = GZ
    G_matrix(4) = GXY
    G_matrix(5) = GYZ
    G_matrix(6) = GZX

    WRITE(6, *) 'COOLING RATES'
    CALL CALCULATE_COOLING_RATES(&
        detuning, kappa, &
        GX, GY, GZ, &
        omega_x, omega_y, omega_z, &
        cooling_rate_x, cooling_rate_y, cooling_rate_z)

    ! X
    WRITE(6, *) 'cooling_rate_x', cooling_rate_x
    WRITE(6, *) 'Gamma_M', Gamma_M
    WRITE(6, *) 'WKX0', WKX0

    WRITE(6, *) 'X phonons: at room bath_temperature; at equilibrium'
    C1 = kB * bath_temperature / hbar / omega_x
    C2 = C1 * Gamma_M / (abs(cooling_rate_x) + Gamma_M)
    WRITE(6, 100) C1, abs(C2)
    phonons_from_cooling_formula_array(1) = abs(C2)

    WRITE(6, *) 'X temperature: at room bath_temperature; at equilibrium'
    WRITE(6, *) bath_temperature, bath_temperature * Gamma_M / (Gamma_M + abs(cooling_rate_x))

    ! Y
    WRITE(6, *) 'cooling_rate_y', cooling_rate_y
    WRITE(6, *) 'Gamma_M', Gamma_M
    WRITE(6, *) 'Y0', Y0

    WRITE(6, *) 'Y phonons: at room bath_temperature; at equilibrium'
    C1 = kB * bath_temperature / hbar / omega_y
    C2 = C1 * Gamma_M / (abs(cooling_rate_y) + Gamma_M)
    WRITE(6, 100) C1, abs(C2)
    phonons_from_cooling_formula_array(2) = abs(C2)

    WRITE(6, *) 'Y temperature: at room bath_temperature; at equilibrium'
    WRITE(6, *) bath_temperature, bath_temperature * Gamma_M / (Gamma_M + abs(cooling_rate_y))

    ! Z
    WRITE(6, *) 'cooling_rate_z', cooling_rate_z
    WRITE(6, *) 'Gamma_M', Gamma_M
    WRITE(6, *) 'Z0', Z0

    WRITE(6, *) 'Z phonons: at room bath_temperature; at equilibrium'
    C1 = kB * bath_temperature / hbar / omega_z
    C2 = C1 * Gamma_M / (abs(cooling_rate_z) + Gamma_M)
    WRITE(6, 100) C1, abs(C2)
    phonons_from_cooling_formula_array(3) = abs(C2)

    WRITE(6, *) 'Z temperature: at room bath_temperature; at equilibrium'
    WRITE(6, *) bath_temperature, bath_temperature * Gamma_M / (Gamma_M + abs(cooling_rate_z))

    100 FORMAT(4D16.8)

    RETURN
END


SUBROUTINE CALCULATE_COOLING_RATES(&
    detuning, kappa, &
    GX, GY, GZ, &
    omega_x, omega_y, omega_z, &
    cooling_rate_x, cooling_rate_y, cooling_rate_z)
    ! """
    ! subroutine below obtains the optomechanical parameters
    ! Note that detuning effectively is detuning + AOPT
    ! It is corrected by zeroth order optical contribution in the linearised case
    ! """

    IMPLICIT NONE
    double precision::detuning, kappa
    double precision::GX, GY, GZ
    double precision::omega_x, omega_y, omega_z
    double precision::cooling_rate_x, cooling_rate_y, cooling_rate_z

    double precision::pi, half_kappa
    double precision::C_positive, C_negative

    pi = dacos(-1.d0)
    half_kappa = 0.5d0 * kappa
    WRITE(6, *) 'GX /Hz', GX/2/pi
    WRITE(6, *) 'GY /Hz', GY/2/pi
    WRITE(6, *) 'GZ /Hz', GZ/2/pi

    ! NOTE: here we neglect opto shift
    ! X cooling
    C_positive = 1.d0 / ((detuning + omega_x)**2 + half_kappa**2)
    C_negative = 1.d0 / ((detuning - omega_x)**2 + half_kappa**2)
    cooling_rate_x = -GX**2 * kappa * (C_positive - C_negative)

    ! Y cooling
    C_positive = 1.d0 / ((detuning + omega_y)**2 + half_kappa**2)
    C_negative = 1.d0 / ((detuning - omega_y)**2 + half_kappa**2)
    cooling_rate_y = -GY**2 * kappa * (C_positive - C_negative)

    ! Z cooling
    C_positive = 1.d0 / ((detuning + omega_z)**2 + half_kappa**2)
    C_negative = 1.d0 / ((detuning - omega_z)**2 + half_kappa**2)
    cooling_rate_z = -GZ**2 * kappa * (C_positive - C_negative)

    RETURN
END
