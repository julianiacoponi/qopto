! This code has 1 beam plus 3 dimensions
! Calculates quantum noise spectra.
! USES real mean fields
! works out position and momentum optical quadrature
! Gamma_M = mechanical damping due to gas, also noise from fluctuation dissipation.


IMPLICIT NONE

! parameter file with input values
integer::N_period, N_total
double precision::bead_radius, bead_density
double precision::vacuum_permittivity, relative_permittivity
double precision::speed_of_light, hbar, kB, bath_temperature, Gravity
double precision::linewidth_k, beam_waist_X, beam_waist_Y
double precision::cavity_waist, cavity_length, Finesse, air_pressure
double precision::tweezer_input_power, DelFSR, theta_0
double precision::detuning_Hz, detuning_Hz_Antonio_0_91, detuning_Hz_Antonio_1_825
double precision::X0, Y0, Z0
INCLUDE 'CSCAVITY.h'

integer::ii, jj, nn
double precision::pi, pi2
double precision::theta_homodyne_0
double precision::cavity_photons
double precision::theta
double precision::G_matrix(6)
double precision::Rayleigh_range, bead_mass, polarisability
double precision::E_tweezer, E_cavity
double precision::detuning, kappa, Gamma_M
double precision::omega_x, omega_y, omega_z
double precision::max_omega, min_omega, omega_increment, omega_value, omega_kHz
double precision::theta_homodyne
double precision::half_kappa, half_Gamma_M
double precision::E_drive
double precision::S_XX_value, S_YY_value, S_ZZ_value
double complex::S_XY_value, S_YX_value
double precision::S_homodyne_value, S_homodyne_negative, S_heterodyne_value
double precision::num_phonons_X, num_phonons_Y, num_phonons_Z
double precision::num_photons, num_phonons
double precision::area_under_spectrum
double precision::phonons_from_cooling_formula_array(3), phonons_array(3), x_y_phonons

! loop this many samples over omega
integer::num_X0_samples, num_omega_samples
PARAMETER(num_X0_samples=50, num_omega_samples=20000)
double precision, DIMENSION(num_omega_samples)::&
    omega_array, &
    S_XX_array, S_YY_array, S_ZZ_array, &
    S_homodyne_array, S_heterodyne_array

double precision::X0_value, lambda_coeff
double precision::omega_x_kHz, omega_y_kHz, omega_z_kHz

integer::equation_choice, detuning_choice
integer::loop_debug, extra_debug, bead_debug, equilibrium_debug

! integer labels for WRITE for the .dat output files
integer::spectra_XYZ, optical_spring_XYZ, omegas_XYZ, peaks_XYZ, cavity_parameters
integer::stdout
integer::FTSXY, FTanOPT, GCOUPLE, PHONS

! various debug flags :-)
loop_debug = 1
extra_debug = 0
bead_debug = 0
equilibrium_debug = 0

spectra_XYZ = 1
optical_spring_XYZ = 2
omegas_XYZ = 3
peaks_XYZ = 4
cavity_parameters = 5
stdout = 6
FTSXY = 7
FTanOPT = 8
GCOUPLE = 9
PHONS = 10

WRITE(stdout, *)  "Which equations to use?"
WRITE(stdout, *) '1 = Full 3D'
WRITE(stdout, *) '2 = Simplified Eq. 18'
READ(*, *)  equation_choice

WRITE(stdout, *)  "Which detuning to use?"
WRITE(stdout, *) '1 = default (-176kHz)'
WRITE(stdout, *) '2 = from Antonio`s data 0.91 (-178.36kHz)'
WRITE(stdout, *) '3 = from Antonio`s data 1.825 (-357.7kHz)'
! NOTE: the above assumes a kappa/2 = 196kHz, whereas the calcualted value here is 204kHz
READ(*, *)  detuning_choice

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
415 FORMAT(4ES23.15)
503 FORMAT(5ES11.3)

! Here plot optical field amplitudes
OPEN(spectra_XYZ, file='spectra_XYZ.dat', status='unknown')
WRITE(spectra_XYZ, *) 'kHz ', 'S_XX ', 'S_YY ', 'S_ZZ'

OPEN(optical_spring_XYZ, file='optical_spring_XYZ.dat', status='unknown')
WRITE(optical_spring_XYZ, *) 'X0/lambda ', 'X_deviation ', 'Y_deviation ', 'Z_deviation'

OPEN(omegas_XYZ, file='omegas_XYZ.dat', status='unknown')
WRITE(omegas_XYZ, *) 'X0/lambda ', 'omega_x ', 'omega_y ', 'omega_z'

OPEN(peaks_XYZ, file='peaks_XYZ.dat', status='unknown')
WRITE(peaks_XYZ, *) 'X0/lambda ', 'S_XX_peaks ', 'S_YY_peaks ', 'S_ZZ_peaks'

OPEN(cavity_parameters, file='cavity_parameters.dat', status='unknown')
WRITE(cavity_parameters, *) "bead_radius", bead_radius
WRITE(cavity_parameters, *) "bead_density", bead_density
WRITE(cavity_parameters, *) "vacuum_permittivity", vacuum_permittivity
WRITE(cavity_parameters, *) "relative_permittivity", relative_permittivity
WRITE(cavity_parameters, *) "speed_of_light", speed_of_light
WRITE(cavity_parameters, *) "hbar", hbar
WRITE(cavity_parameters, *) "kB", kB
WRITE(cavity_parameters, *) "Gravity", Gravity
WRITE(cavity_parameters, *) "linewidth_k", linewidth_k
WRITE(cavity_parameters, *) "beam_waist_X", beam_waist_X
WRITE(cavity_parameters, *) "beam_waist_Y", beam_waist_Y
WRITE(cavity_parameters, *) "cavity_waist", cavity_waist
WRITE(cavity_parameters, *) "cavity_length", cavity_length
WRITE(cavity_parameters, *) "Finesse", Finesse
WRITE(cavity_parameters, *) "air_pressure", air_pressure
WRITE(cavity_parameters, *) "tweezer_input_power", tweezer_input_power
WRITE(cavity_parameters, *) "detuning_Hz", detuning_Hz
WRITE(cavity_parameters, *) "detuning_Hz_Antonio_0_91", detuning_Hz_Antonio_0_91
WRITE(cavity_parameters, *) "detuning_Hz_Antonio_1_825", detuning_Hz_Antonio_1_825
WRITE(cavity_parameters, *) "DelFSR", DelFSR
WRITE(cavity_parameters, *) "theta_0", theta_0
WRITE(cavity_parameters, *) "X0", X0
WRITE(cavity_parameters, *) "Y0", Y0
WRITE(cavity_parameters, *) "Z0", Z0
WRITE(cavity_parameters, *) "equation_choice", equation_choice
WRITE(cavity_parameters, *) "detuning_choice", detuning_choice

IF (detuning_choice == 1) THEN
    detuning = detuning_Hz
ELSE IF (detuning_choice == 2) THEN
    detuning = detuning_Hz
ELSE IF (detuning_choice == 3) THEN
    detuning = detuning_Hz_Antonio_1_825
END IF
WRITE(cavity_parameters, *) "detuning", detuning

OPEN(FTSXY, file="FTSXY.dat", status="unknown")
OPEN(FTanOPT, file="FTanOPT.dat", status="unknown")
OPEN(GCOUPLE, file="GCOUPLE.dat", status="unknown")
OPEN(PHONS, file="PHONS.dat", status="unknown")

pi = dacos(-1.d0)
pi2 = 2.d0 * pi

theta_homodyne_0 = 0.12d0
theta_homodyne = theta_homodyne_0 * pi

detuning = detuning * pi2
theta = theta_0 * pi
! loop over the X0 equilibrium position
DO ii=1, num_X0_samples
    ! increase x0 in increments of num_X0_samples segments of a quarter wavelength (0.25 * lambda)
    lambda_coeff = ii * 0.25d0 / num_X0_samples

    ! if ignoring X0 sampling, just set it to be the node
    IF (num_X0_samples == 1) THEN
        lambda_coeff = 0.25d0
    END IF

    ! wavelength = 1064 nanometres
    X0_value = lambda_coeff * 1.064d-6

    ! note BEAD routine uses x0 but not theta
    CALL CALCULATE_BEAD_PARAMETERS(&
        Rayleigh_range, bead_mass, polarisability, &
        E_tweezer, E_cavity, &
        kappa, Gamma_M, &
        bead_debug, stdout)

    CALL CALCULATE_EQUILIBRIUM_PARAMETERS(&
        theta, X0_value, &
        G_matrix, &
        Rayleigh_range, bead_mass, polarisability, &
        E_tweezer, E_cavity, &
        detuning, kappa, Gamma_M, &
        omega_x, omega_y, omega_z, &
        num_phonons_X, num_phonons_Y, num_phonons_Z, &
        phonons_from_cooling_formula_array, &
        equilibrium_debug, stdout)

    ! only write to certain files when at the node of 0.25lambda
    IF (lambda_coeff == 0.25) THEN
        WRITE(GCOUPLE, 503) theta, (abs(G_matrix(ii)), nn=1, 6)
    END IF

    half_kappa = kappa * 0.5d0
    half_Gamma_M = Gamma_M * 0.5d0

    E_drive = 0.5d0 * polarisability * E_tweezer * E_cavity * sin(theta)
    ! Ed = Ed / hbar
    cavity_photons = (E_drive / hbar)**2 * cos(linewidth_k * X0_value)**2
    cavity_photons = cavity_photons / (half_kappa**2 + detuning**2)

    IF (loop_debug == 1) THEN
        WRITE(stdout, *)
        WRITE(stdout, *) 'X0 loop', ii, 'of', num_X0_samples
        WRITE(stdout, *) 'lambda_coeff = ', lambda_coeff
        WRITE(stdout, *) 'theta /π', theta / pi
        WRITE(stdout, *) 'kappa /Hz', kappa / pi2
        WRITE(stdout, *) 'Gamma_M /Hz', Gamma_M / pi2
        WRITE(stdout, *) 'Number of photons in cavity', cavity_photons
    END IF

    ! shot noise
    num_photons = 0.d0

    ! open loop over frequency for noise spectra
    ! e.g. S_XX(omega) = FT(autocorrelation<X(t).X^dag(t + tau)>)
    S_homodyne_array = 0.d0

    ! expecting peaks at omega_x, so sample twice that range

    max_omega = 2.d0 * omega_x * 1.001

    ! either use: postive and negative omega
    min_omega = -max_omega
    ! or: focus only on positive omega
    min_omega = 0.d0

    omega_increment = (max_omega - min_omega) / num_omega_samples

    IF (loop_debug == 1) THEN
        WRITE(stdout, *) 'omega range = ', min_omega, ' to ', max_omega
        WRITE(stdout, *) 'omega increment = ', omega_increment
    END IF

    DO jj=1, num_omega_samples
        omega_value = min_omega + (jj - 1) * omega_increment
        ! store frequency for integration
        omega_array(jj) = omega_value
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
            -omega_value, &
            omega_x, omega_y, omega_z, &
            S_XX_value, S_YY_value, S_ZZ_value, &
            S_XY_value, S_YX_value, & ! cross-correlation spectra (can be complex)
            S_homodyne_value, S_heterodyne_value)

        ! update homodyne and symm disp.
        ! TODO: clarify what the homodyne calc is doing here...
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
            omega_value, &
            omega_x, omega_y, omega_z, &
            S_XX_value, S_YY_value, S_ZZ_value, &
            S_XY_value, S_YX_value, &
            S_homodyne_value, S_heterodyne_value)

        ! update
        S_homodyne_negative = 0.5d0 * (S_homodyne_negative + S_homodyne_value)
        omega_kHz = omega_value / pi2 * 1.d-3

        ! only record the spectra at the node
        IF (lambda_coeff == 0.25) THEN
            WRITE(spectra_XYZ, 415) omega_kHz, S_XX_value, S_YY_value, S_ZZ_value
            WRITE(FTSXY, 503) omega_kHz, S_XY_value, S_YX_value
            WRITE(FTanOPT, 503) omega_kHz, S_heterodyne_value, S_homodyne_value
        END IF

        S_XX_array(jj) = S_XX_value
        S_YY_array(jj) = S_YY_value
        S_ZZ_array(jj) = S_ZZ_value
        S_heterodyne_array(jj) = S_heterodyne_value
        ! to find optimal squeezing
        S_homodyne_array(jj) = S_homodyne_negative

    ! end of loop over omega samples
    END DO

    IF (extra_debug == 1) THEN
        ! verify indexing using maxloc() gives the same value as maxval()
        WRITE(stdout, *) &
            'Peak spectral values (S_XX, S_YY, S_ZZ), using maxval()', &
            maxval(S_XX_array), &
            maxval(S_YY_array), &
            maxval(S_ZZ_array)
        WRITE(stdout, *) &
            'Peak spectral values (S_XX, S_YY, S_ZZ), using the index from maxloc()', &
            S_XX_array(maxloc(S_XX_array)), &
            S_YY_array(maxloc(S_YY_array)), &
            S_ZZ_array(maxloc(S_ZZ_array))
    END IF

    omega_x_kHz = omega_x / pi2 * 1.d-3
    omega_y_kHz = omega_y / pi2 * 1.d-3
    omega_z_kHz = omega_z / pi2 * 1.d-3

    WRITE(optical_spring_XYZ, 415) &
        lambda_coeff, &
        omega_x_kHz - omega_array(maxloc(S_XX_array)) / pi2 * 1.d-3 , &
        omega_y_kHz - omega_array(maxloc(S_YY_array)) / pi2 * 1.d-3, &
        omega_z_kHz - omega_array(maxloc(S_ZZ_array)) / pi2 * 1.d-3

    WRITE(omegas_XYZ, 415) lambda_coeff, omega_x_kHz, omega_y_kHz, omega_z_kHz
    WRITE(peaks_XYZ, 415) &
        lambda_coeff, &
        omega_array(maxloc(S_XX_array)) / pi2 * 1.d-3, &
        omega_array(maxloc(S_YY_array)) / pi2 * 1.d-3, &
        omega_array(maxloc(S_ZZ_array)) / pi2 * 1.d-3

    IF (extra_debug == 1) THEN
        num_phonons = (omega_x + detuning)**2 + half_kappa**2
        num_phonons = num_phonons / 4. / omega_x / (-detuning)
        WRITE(stdout, *) 'X: back action limited phonons', num_phonons

        num_phonons = (omega_y + detuning)**2 + half_kappa**2
        num_phonons = num_phonons / 4. / omega_y / (-detuning)
        WRITE(stdout, *) 'Y: back action limited phonons', num_phonons

        num_phonons = (omega_z + detuning)**2 + half_kappa**2
        num_phonons = num_phonons / 4. / omega_z / (-detuning)
        WRITE(stdout, *) 'Z: back action limited phonons', num_phonons

        ! TODO: figure out why this gives different areas to original for S_ZZ and S_het?
        ! integrate and normalise the quantum noise spectra
        CALL INTEGRATE_SPECTRUM(num_omega_samples, omega_increment, S_XX_array, area_under_spectrum)
        ! Area corresponds to 2n+1 so convert to get n
        phonons_array(1) = 0.5d0 * (area_under_spectrum - 1.d0)
        WRITE(stdout, *) 'X phonons from: cooling formula, S_XX FT'
        WRITE(stdout, 415) phonons_from_cooling_formula_array(1), phonons_array(1)

        CALL INTEGRATE_SPECTRUM(num_omega_samples, omega_increment, S_YY_array, area_under_spectrum)
        phonons_array(2) = 0.5d0 * (area_under_spectrum - 1.d0)
        WRITE(stdout, *) 'Y phonons from: cooling formula, S_YY FT'
        WRITE(stdout, 415) phonons_from_cooling_formula_array(2), phonons_array(2)

        CALL INTEGRATE_SPECTRUM(num_omega_samples, omega_increment, S_ZZ_array, area_under_spectrum)
        phonons_array(3) = 0.5d0 * (area_under_spectrum - 1.d0)
        WRITE(stdout, *) 'Z phonons from: cooling formula, S_ZZ FT'
        WRITE(stdout, 415) phonons_from_cooling_formula_array(3), phonons_array(3)

        CALL INTEGRATE_SPECTRUM(num_omega_samples, omega_increment, S_heterodyne_array, area_under_spectrum)
        ! Area num_phonons corresponds to 2(nx+ny)+2 so DIFFERENT CONVERSION  to get nx+ny
        x_y_phonons = 0.5d0 * (area_under_spectrum - 2.d0)
        WRITE(stdout, *) 'X+Y phonons from S_heterodyne FT'
        WRITE(stdout, 415) x_y_phonons
    END IF

    IF (lambda_coeff == 0.25) THEN
        WRITE(PHONS, 503) &
            detuning, &
            phonons_array(1), phonons_array(2), x_y_phonons, &
            phonons_from_cooling_formula_array(1), phonons_from_cooling_formula_array(2), phonons_from_cooling_formula_array(3)
    END IF

! end of loop over X0 samples
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

    double precision::pi
    double precision::G_average
    double precision::omega_heterodyne, omega_addition
    double complex::gxy
    double complex::chi_C, chi_C_star
    double complex::chi_X, chi_X_star
    double complex::chi_Y, chi_Y_star
    double complex::chi_Z, chi_Z_star
    double complex::a_hat_vector(1, N_total), a_dag_vector(1, N_total)
    double complex::X_vector(1, N_total), Y_vector(1, N_total), Z_vector(1, N_total)

    pi = dacos(-1.d0)

    S_homodyne_value = 0.d0
    X_vector = 0.d0
    Y_vector = 0.d0
    Z_vector = 0.d0

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
        X_vector, Y_vector, Z_vector)

    ! X,Y,Z_vector have 8 values, ordered (a,a*,bx,bx*,by,by*,bz,bz*)
    S_XX_value = abs(X_vector(1, 1))**2
    S_XX_value = S_XX_value + (num_phonons_X + 1) * abs(X_vector(1, 3))**2 + num_phonons_X * abs(X_vector(1, 4))**2
    S_XX_value = S_XX_value + (num_phonons_Y + 1) * abs(X_vector(1, 5))**2 + num_phonons_Y * abs(X_vector(1, 6))**2
    S_XX_value = S_XX_value + (num_phonons_Z + 1) * abs(X_vector(1, 7))**2 + num_phonons_Z * abs(X_vector(1, 8))**2

    S_YY_value = abs(Y_vector(1, 1))**2
    S_YY_value = S_YY_value + (num_phonons_X + 1) * abs(Y_vector(1, 3))**2 + num_phonons_X * abs(Y_vector(1, 4))**2
    S_YY_value = S_YY_value + (num_phonons_Y + 1) * abs(Y_vector(1, 5))**2 + num_phonons_Y * abs(Y_vector(1, 6))**2
    S_YY_value = S_YY_value + (num_phonons_Z + 1) * abs(Y_vector(1, 7))**2 + num_phonons_Z * abs(Y_vector(1, 8))**2

    S_ZZ_value = abs(Z_vector(1, 1))**2
    S_ZZ_value = S_ZZ_value + (num_phonons_X + 1) * abs(Z_vector(1, 3))**2 + num_phonons_X * abs(Z_vector(1, 4))**2
    S_ZZ_value = S_ZZ_value + (num_phonons_Y + 1) * abs(Z_vector(1, 5))**2 + num_phonons_Y * abs(Z_vector(1, 6))**2
    S_ZZ_value = S_ZZ_value + (num_phonons_Z + 1) * abs(Z_vector(1, 7))**2 + num_phonons_Z * abs(Z_vector(1, 8))**2

    ! WORK OUT CROSS SPECTRA PSD S_XY and S_YX
    S_XY_value = X_vector(1, 1) * conjg(Y_vector(1, 1))
    S_XY_value = S_XY_value + (num_phonons_X + 1) * X_vector(1, 3) * conjg(Y_vector(1, 3))
    S_XY_value = S_XY_value + num_phonons_X * X_vector(1, 4) * conjg(Y_vector(1, 4))
    S_XY_value = S_XY_value + (num_phonons_Y + 1) * X_vector(1, 5) * conjg(Y_vector(1, 5))
    S_XY_value = S_XY_value + num_phonons_Y * X_vector(1, 6) * conjg(Y_vector(1, 6))
    S_XY_value = S_XY_value + (num_phonons_Z + 1) * X_vector(1, 7) * conjg(Y_vector(1, 7))
    S_XY_value = S_XY_value + num_phonons_Z * X_vector(1, 8) * conjg(Y_vector(1, 8))

    S_YX_value = Y_vector(1, 1) * conjg(X_vector(1, 1))
    S_YX_value = S_YX_value + (num_phonons_X + 1) * Y_vector(1, 3) * conjg(X_vector(1, 3))
    S_YX_value = S_YX_value + num_phonons_X * Y_vector(1, 4) * conjg(X_vector(1, 4))
    S_YX_value = S_YX_value + (num_phonons_Y + 1) * Y_vector(1, 5) * conjg(X_vector(1, 5))
    S_YX_value = S_YX_value + num_phonons_Y * Y_vector(1, 6) * conjg(X_vector(1, 6))
    S_YX_value = S_YX_value + (num_phonons_Z + 1) * Y_vector(1, 7) * conjg(X_vector(1, 7))
    S_YX_value = S_YX_value + num_phonons_Z * Y_vector(1, 8) * conjg(X_vector(1, 8))

    S_XY_value = 0.5 * (S_XY_value + S_YX_value)

    IF (equation_choice == 2) THEN
        gxy = a_hat_vector(1, 4)
        ! IN THIS VERSION WE WORK OUT Eq. 18,  simplified version of cross-correlation
        S_YX_value = (S_XX_value - S_YY_value) * gxy / (omega_y - omega_x)
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
        X_vector, Y_vector, Z_vector)

    CALL CALCULATE_HETERODYNE_FROM_OPTICAL_NOISE_VECTORS(&
        N_total, &
        theta, &
        num_photons, &
        num_phonons_X, num_phonons_Y, num_phonons_Z, &
        a_hat_vector, a_dag_vector, &
        S_heterodyne_value)

    ! WRITE(stdout, *) omega_heterodyne, S_heterodyne_value, (S_heterodyne_value) / abs(chi_C)**2 / G_matrix(1)**2
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
    X_vector, Y_vector, Z_vector)
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
    double complex, INTENT(OUT)::X_vector(1, N_total), Y_vector(1, N_total), Z_vector(1, N_total)

    integer::ii
    double precision::g_xY, g_yY, g_zP
    double precision::gxy, gzx, gyz
    double complex::imag_num, one
    double complex::beta_x, beta_y, beta_z
    double complex::G_coeff_XY, G_coeff_ZX, G_coeff_YZ
    double complex::G_coeff_YX, G_coeff_XZ, G_coeff_ZY
    double complex::RXY, RZX, RYZ
    double complex::RYX, RXZ, RZY
    double complex::M_X, M_Y, M_Z
    double complex::CSUM, CNORM
    double complex::mu_X, mu_Y, mu_Z
    double complex::eta_C_0, eta_C_neg_half_pi, eta_C_pos_half_pi
    double complex::N0X(N_total), N0Y(N_total), N0Z(N_total)

    ! NOTE: x,y coupled to Y while z couples to P
    ! g_xY = gx, g_yY = gy, g_zP = gz
    ! g_xP = gyP = gzY = 0
    g_xY = G_matrix(1)
    g_yY = G_matrix(2)
    g_zP = G_matrix(3)
    gxy = G_matrix(4)
    gyz = G_matrix(5)
    gzx = G_matrix(6)

    a_hat_vector = cmplx(0.d0, 0.d0)
    a_dag_vector = cmplx(0.d0, 0.d0)

    X_vector = cmplx(0.d0, 0.d0)
    Y_vector = cmplx(0.d0, 0.d0)
    Z_vector = cmplx(0.d0, 0.d0)

    N0X = cmplx(0.d0, 0.d0)
    N0Y = cmplx(0.d0, 0.d0)
    N0Z = cmplx(0.d0, 0.d0)

    imag_num = cmplx(0.d0, 1.d0)
    one = cmplx(1.d0, 0.0d0)

    ! SUSCEPTIBILITIES
    ! NOTE: using this eta_C_phi is implicitly assuming the x,y coupling to Y only, and z coupling to P
    ! eta_C_phi = e^(-i.phi) * chi_C - e^(i.phi) * chi_C_star
    eta_C_0 = chi_C - chi_C_star
    eta_C_neg_half_pi = imag_num * (chi_C + chi_C_star)
    eta_C_pos_half_pi = -imag_num * (chi_C + chi_C_star)
    mu_X = chi_X - chi_X_star
    mu_Y = chi_Y - chi_Y_star
    mu_Z = chi_Z - chi_Z_star

    ! coeff of X-Y coupling- Combines direct and indirect paths
    G_coeff_XY = gxy + imag_num * eta_C_0 * g_xY * g_yY
    G_coeff_YX = gxy + imag_num * eta_C_0 * g_xY * g_yY
    G_coeff_XZ = gzx + imag_num * eta_C_neg_half_pi * g_zP * g_xY
    G_coeff_ZX = gzx + imag_num * eta_C_pos_half_pi * g_zP * g_xY
    G_coeff_YZ = gyz + imag_num * eta_C_neg_half_pi * g_yY * g_zP
    G_coeff_ZY = gyz + imag_num * eta_C_pos_half_pi * g_yY * g_zP

    ! NORMALIZATIONS
    M_X = 1.d0 + g_xY**2 * mu_X * eta_C_0
    M_Y = 1.d0 + g_yY**2 * mu_Y * eta_C_0
    M_Z = 1.d0 + g_zP**2 * mu_Z * eta_C_0

    ! this is the coefficient of the amplitude (phase) quadrature's noise vector for x,y (z)
    beta_x = imag_num * sqrt(kappa) * mu_X * g_xY
    beta_y = imag_num * sqrt(kappa) * mu_Y * g_yY
    beta_z = imag_num * sqrt(kappa) * mu_Z * g_zP

    ! 1D positional noise vectors, i.e. qj_1D
    ! qj_1D = M_j^-1 * (sqrt(Gamma_M).qj_thermal + i.mu_j.sqrt(kappa)*(gjY.Y_in + gjP.P_in)
    ! M_j = 1 + mu_j.eta_C_0.gj.gj*
    ! g_j = g_jY + i.g_jP
    ! qj_thermal = chi_j.bj + chi_j*.bj* = thermal noise for each axis

    ! optical noises
    ! Y_in = chi_C.a + chi_C*.a*
    ! P_in = i(chi_C*.a* - chi_C.a)

    ! x_1D = zero-th order X noise vector; weights of a, a*, bx, bx*
    N0X(1) = beta_x * chi_C / M_X
    N0X(2) = beta_x * chi_C_star / M_X
    N0X(3) = sqrt(Gamma_M) * chi_X / M_X
    N0X(4) = sqrt(Gamma_M) * chi_X_star / M_X

    ! y_1D = zero-th order Y noise vector; weights of a, a*, by, by*
    N0Y(1) = beta_y * chi_C / M_Y
    N0Y(2) = beta_y * chi_C_star / M_Y
    N0Y(5) = sqrt(Gamma_M) * chi_Y / M_Y
    N0Y(6) = sqrt(Gamma_M) * chi_Y_star / M_Y

    ! z_1D = zero-th order Z noise vector; weights of a, a*, bz, bz*
    N0Z(1) = -imag_num * beta_z * chi_C / M_Z
    N0Z(2) = imag_num * beta_z * chi_C_star / M_Z
    N0Z(7) = sqrt(Gamma_M) * chi_Z / M_Z
    N0Z(8) = sqrt(Gamma_M) * chi_Z_star / M_Z

    ! Entries in matrix R from matrix equation r_3D = R.r_3D + r_1D
    ! r_3D = (x_3D, y_3D, z_3D) - column vector of position operators for each axis considered in 3D
    ! r_1D = (x_1D, y_1D, z_1D) - column vector of position operators for each axis considered in 1D

    ! i.e. R is a 3x3 matrix:
    !  0  RXY RXZ
    ! RYX  0  RYZ
    ! RZX RZY  0
    RXY = imag_num * mu_X * G_coeff_XY / M_X
    RXZ = imag_num * mu_X * G_coeff_XZ / M_X
    RYX = imag_num * mu_Y * G_coeff_YX / M_Y
    RYZ = imag_num * mu_Y * G_coeff_YZ / M_Y
    RZX = imag_num * mu_Z * G_coeff_ZX / M_Z
    RZY = imag_num * mu_Z * G_coeff_ZY / M_Z

    ! So, find inverse of final matrix A = 1 - R from matrix equation r_3D = A^-1.r_1D
    ! A^-1 = adj(A)/det(A) = CSUM/CNORM

    ! det(A)
    CNORM = 1.d0 - ((RYX * RXY) + (RYZ * RZY)+ (RXZ * RZX)) - ((RXY * RYZ * RZX) + (RXZ * RZY * RYX))

    ! ADD 3D BACK-ACTION TERMS
    DO ii=1, N_total
        ! 1st row of adj(A), applies to x_3D
        CSUM = (1.d0 - RYZ * RZY) * N0X(ii) + (RXY + RXZ * RZY) * N0Y(ii) + (RXZ + RXY * RYZ) * N0Z(ii)
        X_vector(1, ii) = X_vector(1, ii) + CSUM / CNORM

        ! 2nd row of adj(A), applies to y_3D
        CSUM = (RYX + RYZ * RZX) * N0X(ii) + (1.d0 - RXZ * RZX) * N0Y(ii) + (RYZ + RYX * RXZ) * N0Z(ii)
        Y_vector(1, ii) = Y_vector(1, ii) + CSUM / CNORM

        ! 3rd row of adj(A), applies to z_3D
        CSUM =  (RZX + RZY * RYX) * N0X(ii) + (RZY + RZX * RXY) * N0Y(ii) + (1.d0 - RXY * RYX) * N0Z(ii)
        Z_vector(1, ii) = Z_vector(1, ii) + CSUM / CNORM
    END DO

    DO ii=1, N_total
        a_hat_vector(1, ii) = imag_num * chi_C * (g_xY * X_vector(1, ii) + g_yY * Y_vector(1, ii))
        a_hat_vector(1, ii) = a_hat_vector(1, ii) - chi_C * g_zP * Z_vector(1, ii)

        a_dag_vector(1, ii) = -imag_num * chi_C_star * (g_xY * X_vector(1, ii) + g_yY * Y_vector(1, ii))
        a_dag_vector(1, ii) = a_dag_vector(1, ii) - chi_C_star * g_zP * Z_vector(1, ii)
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
        a_hat_vector(1, 3) = G_coeff_XY
        a_hat_vector(1, 4) = G_coeff_YX
        a_hat_vector(1, 5) = mu_X / M_X
        a_hat_vector(1, 6) = mu_Y / M_Y
    END IF

    RETURN
END


SUBROUTINE INTEGRATE_SPECTRUM(num_omega_samples, omega_increment, spectrum_array, area_under_spectrum)
    ! """
    ! Integrates a given array of spectrum values using the trapezium rule
    ! """
    IMPLICIT NONE

    integer::num_omega_samples
    double precision::omega_increment
    double precision::spectrum_array(num_omega_samples)
    double precision::area_under_spectrum

    integer::ii
    double precision::pi2, trapeze_base

    pi2 = 2.d0 * dacos(-1.d0)

    ! integrate the position spectrum of bead
    ! quick hack - use trapezoidal rule - improve later?
    ! each trapeze base is the omega_increment i.e. omega_range/num_omega_samples

    ! NOTE: divide by 2π to convert omega to non-angular frequency (Hz)
    ! TODO: explain the division by 2π as a result of Parseval's theorem / "Sum-rule" ?
    ! this means the area is the area under the spectrum curve plotted vs. Hz
    trapeze_base = omega_increment/pi2

    area_under_spectrum = 0.d0
    DO ii=1, num_omega_samples - 1
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
    kappa, Gamma_M, &
    bead_debug, stdout)
    ! """
    ! subroutine below is provided by user
    ! and calculates  relevant parameters
    ! """
    IMPLICIT NONE

    integer::N_period, N_total
    double precision::bead_radius, bead_density
    double precision::vacuum_permittivity, relative_permittivity
    double precision::speed_of_light, hbar, kB, bath_temperature, Gravity
    double precision::linewidth_k, beam_waist_X, beam_waist_Y
    double precision::cavity_waist, cavity_length, Finesse, air_pressure
    double precision::tweezer_input_power, DelFSR, theta_0
    double precision::detuning_Hz, detuning_Hz_Antonio_0_91, detuning_Hz_Antonio_1_825
    double precision::X0, Y0, Z0
    INCLUDE 'CSCAVITY.h'

    double precision::Rayleigh_range, bead_mass, polarisability
    double precision::E_tweezer, E_cavity
    double precision::kappa, Gamma_M
    integer::bead_debug, stdout

    double precision::pi, pi2
    double precision::omega_optical, coeff, cavity_volume, amplitude
    double precision::kappa_in, kappa_nano

    ! zero eq. initial values
    pi = dacos(-1.d0)
    pi2 = 2.d0 * pi

    Rayleigh_range = 0.5d0 * linewidth_k * beam_waist_X * beam_waist_Y
    bead_mass = bead_density * 4.*pi/3. * bead_radius**3
    polarisability = 4.* pi * vacuum_permittivity * (relative_permittivity - 1.) / (relative_permittivity + 2.) * bead_radius**3

    omega_optical = speed_of_light * linewidth_k
    ! add a factor of 4 here. Not in GALAX1-5 routines!!!
    cavity_volume = cavity_length * pi * cavity_waist**2 / 4.d0

    ! Depth of cavity field. Weak and unimportant for CS case
    amplitude = omega_optical * polarisability / 2. / cavity_volume / vacuum_permittivity

    E_tweezer = sqrt(4. * tweezer_input_power / (beam_waist_X * beam_waist_Y * pi * speed_of_light * vacuum_permittivity))
    E_cavity = sqrt(hbar * omega_optical / (2. * cavity_volume * vacuum_permittivity))

    kappa_in = pi * speed_of_light / Finesse / cavity_length
    coeff = linewidth_k * polarisability / vacuum_permittivity / omega_optical**2
    kappa_nano = 4. * coeff**2 * DelFSR * cos(linewidth_k * x0)**2
    kappa = kappa_in + kappa_nano

    ! take usual expression e.g. Levitated review by Li Geraci etc
    ! 1 bar = 10^5 pascal; air_pressure /mbar = 10^2 Pascal
    ! Gamma_M = 1600 * air_pressure / (pi * air_speed * bead_density * bead_radius)
    ! air_speed = 500 ms^-1
    Gamma_M = 1600. * air_pressure / pi / 500. / bead_density / bead_radius

    IF (bead_debug == 1) THEN
        WRITE(stdout, *) '---------------------------------------------'
        WRITE(stdout, *) 'BEAD PARAMETERS'
        WRITE(stdout, *) '---------------------------------------------'
        WRITE(stdout, *) 'Rayleigh_range /m = ', Rayleigh_range
        WRITE(stdout, *) 'bead_mass /kg = ', bead_mass
        WRITE(stdout, *) 'bead polarisability / Farad.m^2 = ', polarisability
        WRITE(stdout, *)
        WRITE(stdout, *) 'Cavity trap, amplitude/2π (zeroth shift) /Hz = '
        WRITE(stdout, *) amplitude / pi2,  amplitude * cos(linewidth_k * x0)**2
        WRITE(stdout, *)
        WRITE(stdout, *) 'E_tweezer /Volts.m^-1 = ', E_tweezer
        WRITE(stdout, *) 'E_cavity /Volts.m^-1 = ', E_cavity
        WRITE(stdout, *)
        WRITE(stdout, *) 'kappa_in /Hz = ', kappa_in / pi2
        WRITE(stdout, *) 'kappa_nano /2π*Hz = ', kappa_nano
        WRITE(stdout, *) 'Optical damping: kappa /Hz = ', kappa / pi2
        WRITE(stdout, *)
        WRITE(stdout, *) 'Mechanical damping: Gamma_M /Hz = ', Gamma_M
    END IF

    RETURN
END


SUBROUTINE CALCULATE_EQUILIBRIUM_PARAMETERS(&
    theta, X0_value, &
    G_matrix, &
    Rayleigh_range, bead_mass, polarisability, &
    E_tweezer, E_cavity, &
    detuning, kappa, Gamma_M, &
    omega_x, omega_y, omega_z, &
    num_phonons_X, num_phonons_Y, num_phonons_Z, &
    phonons_from_cooling_formula_array, &
    equilibrium_debug, stdout)
    ! """
    ! subroutine below obtains the optomechanical parameters
    ! """
    IMPLICIT NONE

    integer::N_period, N_total
    double precision::bead_radius, bead_density
    double precision::vacuum_permittivity, relative_permittivity
    double precision::speed_of_light, hbar, kB, bath_temperature, Gravity
    double precision::linewidth_k, beam_waist_X, beam_waist_Y
    double precision::cavity_waist, cavity_length, Finesse, air_pressure
    double precision::tweezer_input_power, DelFSR, theta_0
    double precision::detuning_Hz, detuning_Hz_Antonio_0_91, detuning_Hz_Antonio_1_825
    double precision::X0, Y0, Z0
    INCLUDE 'CSCAVITY.h'

    double precision::theta, X0_value
    double precision::G_matrix(6)
    double precision::Rayleigh_range, bead_mass, polarisability
    double precision::E_tweezer, E_cavity
    double precision::detuning, kappa, Gamma_M
    double precision::omega_x, omega_y, omega_z
    double precision::num_phonons_X, num_phonons_Y, num_phonons_Z
    ! analytic equilibrium phonon numbers
    double precision::phonons_from_cooling_formula_array(3)
    integer::equilibrium_debug, stdout

    double precision::pi, pi2
    double precision::kX0
    double precision::N_photon
    double precision::photon_field_real, photon_field_imag
    double precision::correctionXY, correctionZ
    double precision::phonons_from_cooling_formula_X, phonons_from_cooling_formula_Y, phonons_from_cooling_formula_Z
    double precision::g_xY, g_yY, g_zP
    double precision::gxy, gzx, gyz
    double precision::X_zpf, Y_zpf, Z_zpf
    double precision::k_X_zpf, k_Y_zpf, k_Z_zpf
    double precision::E_drive
    double precision::half_kappa
    double precision::cooling_rate_x, cooling_rate_y, cooling_rate_z

    pi = dacos(-1.d0)
    pi2 = 2.d0 * pi
    half_kappa = 0.5d0 * kappa

    omega_x = polarisability * E_tweezer**2 / bead_mass / beam_waist_X**2
    omega_y = polarisability * E_tweezer**2 / bead_mass / beam_waist_Y**2
    omega_z = 0.5d0 * polarisability * E_tweezer**2 / bead_mass / Rayleigh_range**2

    ! use negative E_drive, in units of hbar
    ! electric field of the combined cavity and tweezer light
    E_drive = -0.5d0 * polarisability * E_tweezer * E_cavity * sin(theta)

    kX0 = linewidth_k * X0_value
    ! photon number in cavity
    ! real part of photon field
    photon_field_real = detuning * (E_drive / hbar) * cos(kX0) / (half_kappa**2 + detuning**2)
    photon_field_imag = -half_kappa * (E_drive / hbar) * cos(kX0) / (half_kappa**2 + detuning**2)

    N_photon = ((E_drive / hbar) * cos(kX0))**2 / (half_kappa**2 + detuning**2)

    ! ADD CS POTENTIAL CORRECTION to frequency squared
    correctionXY = -2. * E_drive / bead_mass * photon_field_real * linewidth_k**2 * cos(kX0)
    correctionZ = -2. * E_drive / bead_mass * photon_field_real * (linewidth_k - 1.d0 / Rayleigh_range)**2 * cos(kX0)
    omega_x = sqrt(omega_x + correctionXY * sin(theta)**2)
    omega_y = sqrt(omega_y + correctionXY * cos(theta)**2)
    omega_z = sqrt(omega_z + correctionZ)

    ! zero point fluctuations
    X_zpf = sqrt(hbar / (2.d0 * bead_mass * omega_x))
    Y_zpf = sqrt(hbar / (2.d0 * bead_mass * omega_y))
    Z_zpf = sqrt(hbar / (2.d0 * bead_mass * omega_z))

    k_X_zpf = linewidth_k * X_zpf
    k_Y_zpf = linewidth_k * Y_zpf
    k_Z_zpf = (linewidth_k - 1.d0 / Rayleigh_range) * Z_zpf

    ! optomechanical couplings
    g_xY = (E_drive / hbar) * k_X_zpf * sin(kX0) * sin(theta)
    g_yY = (E_drive / hbar) * k_Y_zpf * sin(kX0) * cos(theta)
    g_zP = -(E_drive / hbar) * k_Z_zpf * cos(kX0)

    gxy = (E_drive / hbar) * k_X_zpf * k_Y_zpf * photon_field_real * cos(kX0) * sin(2. * theta)
    gyz = 2.d0 * (E_drive / hbar) * k_Y_zpf * k_Z_zpf * photon_field_imag * sin(kX0) * cos(theta)
    gzx = 2.d0 * (E_drive / hbar) * k_Z_zpf * k_X_zpf * photon_field_imag * sin(kX0) * sin(theta)

    ! assign these couplings to an array
    G_matrix(1) = g_xY
    G_matrix(2) = g_yY
    G_matrix(3) = g_zP
    G_matrix(4) = gxy
    G_matrix(5) = gyz
    G_matrix(6) = gzx

    IF (equilibrium_debug == 1) THEN
        WRITE(stdout, *)
        WRITE(stdout, *) '---------------------------------------------'
        WRITE(stdout, *) 'EQUILIBRIUM PARAMETERS'
        WRITE(stdout, *) '---------------------------------------------'
        WRITE(stdout, *) '(E_drive / hbar) / 2π = ', (E_drive / hbar) / pi2
        WRITE(stdout, *)
        WRITE(stdout, *) 'detuning /Hz = ', detuning / pi2
        WRITE(stdout, *) 'half_kappa /Hz = ', half_kappa / pi2
        WRITE(stdout, *) 'Number of photons in cavity = ', N_photon
        WRITE(stdout, *)
        WRITE(stdout, *) 'omega_x /Hz = ', omega_x / pi2
        WRITE(stdout, *) 'omega_y /Hz = ', omega_y / pi2
        WRITE(stdout, *) 'omega_z /Hz = ', omega_z / pi2
        WRITE(stdout, *)
        WRITE(stdout, *) 'g_xY /Hz', g_xY / pi2
        WRITE(stdout, *) 'g_yY /Hz', g_yY / pi2
        WRITE(stdout, *) 'g_zP /Hz', g_zP / pi2
        WRITE(stdout, *)
        WRITE(stdout, *) 'gxy /Hz', gxy / pi2
        WRITE(stdout, *) 'gyz /Hz', gyz / pi2
        WRITE(stdout, *) 'gzx /Hz', gzx / pi2
    END IF

    CALL CALCULATE_COOLING_RATES(&
        detuning, kappa, &
        g_xY, g_yY, g_zP, &
        omega_x, omega_y, omega_z, &
        cooling_rate_x, cooling_rate_y, cooling_rate_z)

        ! X
        num_phonons_X = kB * bath_temperature / hbar / omega_x
        phonons_from_cooling_formula_X = num_phonons_X * Gamma_M / (abs(cooling_rate_x) + Gamma_M)
        phonons_from_cooling_formula_array(1) = abs(phonons_from_cooling_formula_X)

        ! Y
        num_phonons_Y = kB * bath_temperature / hbar / omega_y
        phonons_from_cooling_formula_Y = num_phonons_Y * Gamma_M / (abs(cooling_rate_y) + Gamma_M)
        phonons_from_cooling_formula_array(2) = abs(phonons_from_cooling_formula_Y)

        ! Z
        num_phonons_Z = kB * bath_temperature / hbar / omega_z
        phonons_from_cooling_formula_Z = num_phonons_Z * Gamma_M / (abs(cooling_rate_z) + Gamma_M)
        phonons_from_cooling_formula_array(3) = abs(phonons_from_cooling_formula_Z)

    IF (equilibrium_debug == 1) THEN
        WRITE(stdout, *)
        WRITE(stdout, *) '---------------------------------------------'
        WRITE(stdout, *) 'COOLING RATES'
        WRITE(stdout, *) '---------------------------------------------'
        WRITE(stdout, *) 'cooling_rate_x', cooling_rate_x
        WRITE(stdout, *) 'Gamma_M', Gamma_M
        WRITE(stdout, *) 'kX0', kX0
        WRITE(stdout, *) 'X phonons: at room bath_temperature; at equilibrium'
        WRITE(stdout, 100) num_phonons_X, abs(phonons_from_cooling_formula_X)
        WRITE(stdout, *) 'X temperature: at room bath_temperature; at equilibrium'
        WRITE(stdout, *) bath_temperature, bath_temperature * Gamma_M / (Gamma_M + abs(cooling_rate_x))
        WRITE(stdout, *)
        WRITE(stdout, *) 'cooling_rate_y', cooling_rate_y
        WRITE(stdout, *) 'Gamma_M', Gamma_M
        WRITE(stdout, *) 'Y0', Y0
        WRITE(stdout, *) 'Y phonons: at room bath_temperature; at equilibrium'
        WRITE(stdout, 100) num_phonons_Y, abs(phonons_from_cooling_formula_Y)
        WRITE(stdout, *) 'Y temperature: at room bath_temperature; at equilibrium'
        WRITE(stdout, *) bath_temperature, bath_temperature * Gamma_M / (Gamma_M + abs(cooling_rate_y))
        WRITE(stdout, *)
        WRITE(stdout, *) 'cooling_rate_z', cooling_rate_z
        WRITE(stdout, *) 'Gamma_M', Gamma_M
        WRITE(stdout, *) 'Z0', Z0
        WRITE(stdout, *) 'Z phonons: at room bath_temperature; at equilibrium'
        WRITE(stdout, 100) num_phonons_Z, abs(phonons_from_cooling_formula_Z)
        WRITE(stdout, *) 'Z temperature: at room bath_temperature; at equilibrium'
        WRITE(stdout, *) bath_temperature, bath_temperature * Gamma_M / (Gamma_M + abs(cooling_rate_z))
    END IF

    100 FORMAT(4D16.8)

    RETURN
END


SUBROUTINE CALCULATE_COOLING_RATES(&
    detuning, kappa, &
    g_xY, g_yY, g_zP, &
    omega_x, omega_y, omega_z, &
    cooling_rate_x, cooling_rate_y, cooling_rate_z)
    ! """
    ! subroutine below obtains the optomechanical parameters
    ! Note that detuning effectively is detuning + AOPT
    ! It is corrected by zeroth order optical contribution in the linearised case
    ! """

    IMPLICIT NONE
    double precision::detuning, kappa
    double precision::g_xY, g_yY, g_zP
    double precision::omega_x, omega_y, omega_z
    double precision::cooling_rate_x, cooling_rate_y, cooling_rate_z

    double precision::pi, half_kappa
    double precision::C_positive, C_negative

    pi = dacos(-1.d0)
    half_kappa = 0.5d0 * kappa

    ! NOTE: here we neglect opto shift
    ! X cooling
    C_positive = 1.d0 / ((detuning + omega_x)**2 + half_kappa**2)
    C_negative = 1.d0 / ((detuning - omega_x)**2 + half_kappa**2)
    cooling_rate_x = -g_xY**2 * kappa * (C_positive - C_negative)

    ! Y cooling
    C_positive = 1.d0 / ((detuning + omega_y)**2 + half_kappa**2)
    C_negative = 1.d0 / ((detuning - omega_y)**2 + half_kappa**2)
    cooling_rate_y = -g_yY**2 * kappa * (C_positive - C_negative)

    ! Z cooling
    C_positive = 1.d0 / ((detuning + omega_z)**2 + half_kappa**2)
    C_negative = 1.d0 / ((detuning - omega_z)**2 + half_kappa**2)
    cooling_rate_z = -g_zP**2 * kappa * (C_positive - C_negative)

    RETURN
END
