! This code has 1 beam plus 3 dimensions
! Calculates quantum noise spectra, and uses them to determine the frequency shift cancellation point
! USES real mean fields
! works out position and momentum optical quadrature
! Gamma_M = mechanical damping due to gas, also noise from fluctuation dissipation.

IMPLICIT NONE

! parameter file with input values
integer::N_total
double precision::hbar, speed_of_light, kB, vacuum_permittivity
double precision::Free_Spectral_Range, bath_temperature, air_pressure, air_speed
double precision::bead_diameter, bead_density, bead_permittivity
double precision::tweezer_input_power, beam_waist_X, beam_waist_Y, theta_tweezer_over_pi
double precision::cavity_waist, cavity_length, Finesse
double precision::half_kappa_exp_kHz, tweezer_wavelength
INCLUDE 'CSCAVITY.h'

! user choices
integer::equation_choice, use_default_script_params, num_y0_samples

! loop this many samples over omega
integer::num_omega_samples, num_kD_slices
PARAMETER(num_omega_samples=14000, num_kD_slices=30)

double precision, DIMENSION(num_omega_samples)::&
    omega_array, &
    S_XX_array, S_YY_array, S_ZZ_array, &
    S_homodyne_array, S_heterodyne_array

double precision, DIMENSION(num_kD_slices**2)::&
    all_half_kappa_kHz, &
    all_detuning_kHz

integer::ii, jj, kD
double precision::pi, pi2, kHz
double precision::y0_value, lambda_coeff, node, half_node, y0_to_record_spectra_at
double precision::linewidth_k
double precision::theta_homodyne_0
double precision::theta_tweezer_radians
double precision::cavity_photons
double precision::G_matrix(6)
double precision::Rayleigh_range, bead_mass, polarisability
double precision::E_tweezer, E_cavity
double precision::detuning, detuning_kHz
double precision::kappa, half_kappa_kHz
double precision::kD_limit, kD_increment
double precision::Gamma_M
double precision::omega_0_prefactor
double precision::omega_x_0, omega_y_0, omega_z_0
double precision::omega_x_CS, omega_y_CS, omega_z_CS
double precision::omega_x_OS, omega_y_OS, omega_z_OS
double precision::omega_x_CS_plus_OS, omega_y_CS_plus_OS, omega_z_CS_plus_OS
double precision::max_omega, min_omega, omega_increment, omega_value, omega_kHz
double precision::theta_homodyne
double precision::E_drive
double precision::S_XX_value, S_YY_value, S_ZZ_value
double complex::S_XY_value, S_YX_value, S_corr_value
double precision::S_homodyne_value, S_homodyne_negative, S_heterodyne_value
double precision::num_phonons_X, num_phonons_Y, num_phonons_Z
double precision::num_photons, num_phonons
double precision::area_under_spectrum
double precision::phonons_from_cooling_formula_array(3), phonons_array(3), x_y_phonons

character(len=75)::phi_c_dir
character(len=25)::point_dir
character(len=3)::point_num_str

! integer labels for WRITE for the .dat output files
integer::QLT_spectra, kappa_delta_point
integer::omega_CS, omega_OS, omega_CS_plus_OS, omega_QLT
integer::deviation_0, deviation_CS, deviation_OS, deviation_CS_plus_OS
integer::OM_couplings, mech_couplings
integer::cavity_parameters
integer::stdout

! various debug flags :-)
integer::loop_debug, extra_debug, bead_debug, equilibrium_debug, cooling_debug
loop_debug = 1
extra_debug = 0
bead_debug = 1
equilibrium_debug = 1
cooling_debug = 0

! NOTE: make sure none of these file-specifying numbers clash
QLT_spectra = 1
kappa_delta_point = 2

omega_CS = 21
omega_OS = 22
omega_CS_plus_OS = 23
omega_QLT = 24

deviation_0 = 31 ! omega_QLT - omega_0
deviation_CS = 32 ! omega_QLT - omega_CS
deviation_OS = 33 ! omega_QLT - omega_OS
deviation_CS_plus_OS = 34 ! omega_QLT - omega_CS_plus_OS

OM_couplings = 41 ! g_xY, g_yY, g_zP
mech_couplings = 42 ! gxy, gyz, gzx

cavity_parameters = 5
stdout = 6

WRITE(stdout, *)  'Which equations to use for Sxy?'
WRITE(stdout, *) '1 = Full 3D'
WRITE(stdout, *) '2 = Simplified Sxy ~ (Sxx - Syy) * Gxy/(omega_y - omega_x)'
READ(*, *) equation_choice

WRITE(stdout, *) &
    'Use default parameters? num_y0_samples = 250' &
    //'; (min_omega, max_omega) = (25, 165)kHz' &
    //'; limit for sampling (kappa/2, Delta) = 300kHz'
WRITE(stdout, *) '1 = Yes'
WRITE(stdout, *) '2 = No'
READ(*, *) use_default_script_params

IF (use_default_script_params == 2) THEN
    WRITE(stdout, *)  'How many y0 positions to sample between antinode and node?'
    READ(*, *) num_y0_samples

    WRITE(stdout, *) 'Start frequency sampling at what value? (kHz)'
    READ(*, *) min_omega

    ! highest mechanical frequency is omega_y ~ 155kHz for these parameters
    WRITE(stdout, *) 'Stop frequency sampling at what value? (kHz)'
    READ(*, *) max_omega
ELSE
    num_y0_samples = 2500
    min_omega = 25
    max_omega = 165
    kD_limit = 300
END IF

! FORMAT can be understood from:
! https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html
! http://www.personal.psu.edu/jhm/f90/lectures/23.html

! PRINT, OPEN, WRITE, CLOSE can be understood from:
! https://www.tutorialspoint.com/fortran/fortran_file_input_output.htm
! https://meteor.geol.iastate.edu/classes/mt227/lectures/Formatting.pdf

! e.g. the below is:
! '10' = up to 8 things written to the same line
! 'ES' = scientific exponential notation (1-9E+power, instead of 0.1-0.9)
! '23' = number of positions to be used
! '15' = number of digits to the right of the decimal point
! NOTE: this is *not* the same as the precision of the number
1015 FORMAT(10ES23.15)

node = 0.25d0
half_node = 0.125d0
y0_to_record_spectra_at = half_node

pi = dacos(-1.d0)
pi2 = 2.d0 * pi
kHz = pi2 * 1.d3

linewidth_k = pi2 / tweezer_wavelength ! m^-1
theta_tweezer_radians = theta_tweezer_over_pi * pi  ! /rad

theta_homodyne_0 = 0.12d0
theta_homodyne = theta_homodyne_0 * pi

omega_increment = (max_omega - min_omega) / num_omega_samples
IF (loop_debug == 1) THEN
    WRITE(stdout, *)
    WRITE(stdout, *) 'omega range /kHz = ', min_omega, ' to ', max_omega
    WRITE(stdout, *) 'number of omega samples = ', num_omega_samples
    WRITE(stdout, *) 'omega increment /kHz = ', omega_increment
END IF

! convert to radians
min_omega = pi2 * min_omega * 1.d3
omega_increment = pi2 * omega_increment * 1.d3

CALL CALCULATE_BEAD_PARAMETERS(&
    linewidth_k, theta_tweezer_radians, &
    Rayleigh_range, bead_mass, polarisability, &
    E_tweezer, E_cavity, E_drive, &
    Gamma_M, &
    omega_0_prefactor, &
    omega_x_0, omega_y_0, omega_z_0, &
    bead_debug, stdout)

! NOTE: rename the final folder name of phi_c_dir (num_kD_pts = num_kD_slices**2 currently)
phi_c_dir = '/Users/julian/Code/qopto/3D/phi-c-freq-output/num_kD_pts-y0-num_y0_samples/'
CALL system('mkdir '//trim(phi_c_dir))
OPEN(cavity_parameters, file=trim(phi_c_dir)//'cavity_parameters.dat', status='unknown')
WRITE(cavity_parameters, *) "equation_choice", equation_choice

WRITE(cavity_parameters, *) "hbar", hbar
WRITE(cavity_parameters, *) "speed_of_light", speed_of_light
WRITE(cavity_parameters, *) "vacuum_permittivity", vacuum_permittivity

WRITE(cavity_parameters, *) "Free_Spectral_Range", Free_Spectral_Range
WRITE(cavity_parameters, *) "air_pressure", air_pressure
WRITE(cavity_parameters, *) "air_speed", air_speed

WRITE(cavity_parameters, *) "bead_diameter", bead_diameter
WRITE(cavity_parameters, *) "bead_density", bead_density
WRITE(cavity_parameters, *) "bead_permittivity", bead_permittivity

WRITE(cavity_parameters, *) "cavity_waist", cavity_waist
WRITE(cavity_parameters, *) "cavity_length", cavity_length
WRITE(cavity_parameters, *) "Finesse", Finesse

WRITE(cavity_parameters, *) "tweezer_wavelength", tweezer_wavelength

WRITE(cavity_parameters, *) "tweezer_input_power", tweezer_input_power
WRITE(cavity_parameters, *) "beam_waist_X", beam_waist_X
WRITE(cavity_parameters, *) "beam_waist_Y", beam_waist_Y
WRITE(cavity_parameters, *) "theta_tweezer_radians", theta_tweezer_radians

WRITE(cavity_parameters, *) "omega_x_0", omega_x_0 / kHz
WRITE(cavity_parameters, *) "omega_y_0", omega_y_0 / kHz
WRITE(cavity_parameters, *) "omega_z_0", omega_z_0 / kHz

! SOME DESCRIPTION OF HOW TO THINK ABOUT WHICH (kappa/2,Delta) POINTS TO MEASURE:
! for omega_y ~ 136kHz, define grid of (kappa/2, detuning) points to find phi_c for
! need (kappa/2, detuning) to be > circle of radius of omega_y
! all_half_kappa_kHz = (&
!     ! right half grid (300, -300) to (150, -10)
!     300, 300, 300, 300, 300, 300, 300 &
!     250, 250, 250, 250, 250, 250, 250 &
!     200, 200, 200, 200, 200, 200, 200 &
!     150, 150, 150, 150, 150, 150, 150 &
!     ! bottom left grid (10, -150) to (125, -300)
!     10, 10, 10, 10, &
!     25, 25, 25, 25, &
!     50, 50, 50, 50, &
!     75, 75, 75, 75, &
!     100, 100, 100, 100, &
!     125, 125, 125, 125, &
!     ! points closer to the omega_x circle
!     140, 120, 100, & ! toward the radius
!     140, 140, 140, & ! kappa/2=140kHz line
!     10, 50, 100, & ! detuning=-140kHz line
!     120, 80, & ! final 2 points
!     ! looking into the antinode cancellation?
!     5)

! all_detunings_kHz = (&
!     ! right half grid
!     10, 50, 100, 150, 200, 250, 300 & ! 300
!     10, 50, 100, 150, 200, 250, 300 & ! 250
!     10, 50, 100, 150, 200, 250, 300 & ! 200
!     10, 50, 100, 150, 200, 250, 300 & ! 150
!     ! bottom left grid
!     150, 200, 250, 300, &
!     150, 200, 250, 300, &
!     150, 200, 250, 300, &
!     150, 200, 250, 300, &
!     150, 200, 250, 300, &
!     150, 200, 250, 300, &
!     ! circle
!     140, 120, 100, & ! toward the radius
!     10, 50, 100, & ! kappa/2=140kHz line
!     140, 140, 140, & ! detuning=-140kHz line
!     80, 120, & ! final 2 points
!     ! looking into the antinode cancellation?
!     140)

! 30 CUSTOM POINTS USED FOR phi-c-freq-output/30pts-y0-2500/:
! ! initial 8, going toward the omega_x radius on kappa/2 = -detuning line
! all_half_kappa_kHz(1) = 250
! all_detuning_kHz(1) = -250

! all_half_kappa_kHz(2) = 200
! all_detuning_kHz(2) = -200

! all_half_kappa_kHz(3) = 150
! all_detuning_kHz(3) = -150

! all_half_kappa_kHz(4) = 140
! all_detuning_kHz(4) = -140

! all_half_kappa_kHz(5) = 130
! all_detuning_kHz(5) = -130

! all_half_kappa_kHz(6) = 120
! all_detuning_kHz(6) = -120

! all_half_kappa_kHz(7) = 110
! all_detuning_kHz(7) = -110

! all_half_kappa_kHz(8) = 100
! all_detuning_kHz(8) = -100

! ! 5 more, surrounding the bottom of the circle more closely
! all_half_kappa_kHz(9) = 10
! all_detuning_kHz(9) = -140

! all_half_kappa_kHz(10) = 50
! all_detuning_kHz(10) = -140

! all_half_kappa_kHz(11) = 100
! all_detuning_kHz(11) = -140

! all_half_kappa_kHz(12) = 90
! all_detuning_kHz(12) = -130

! all_half_kappa_kHz(13) = 80
! all_detuning_kHz(13) = -120

! ! 5 more, surrounding the right side of the circle more closely
! all_half_kappa_kHz(14) = 140
! all_detuning_kHz(14) = -10

! all_half_kappa_kHz(15) = 140
! all_detuning_kHz(15) = -50

! all_half_kappa_kHz(16) = 140
! all_detuning_kHz(16) = -100

! all_half_kappa_kHz(17) = 130
! all_detuning_kHz(17) = -90

! all_half_kappa_kHz(18) = 120
! all_detuning_kHz(18) = -80

! ! 3 more, bottom left
! all_half_kappa_kHz(19) = 50
! all_detuning_kHz(19) = -200

! all_half_kappa_kHz(20) = 100
! all_detuning_kHz(20) = -200

! all_half_kappa_kHz(21) = 150
! all_detuning_kHz(21) = -200

! ! 3 more, top right
! all_half_kappa_kHz(22) = 200
! all_detuning_kHz(22) = -50

! all_half_kappa_kHz(23) = 200
! all_detuning_kHz(23) = -100

! all_half_kappa_kHz(24) = 200
! all_detuning_kHz(24) = -150

! ! antinode cancellation?
! all_half_kappa_kHz(25) = 5
! all_detuning_kHz(25) = -140

! all_half_kappa_kHz(26) = 5
! all_detuning_kHz(26) = -150

! all_half_kappa_kHz(27) = 5
! all_detuning_kHz(27) = -160

! all_half_kappa_kHz(28) = 5
! all_detuning_kHz(28) = -170

! ! 2 more, the experimental values
! all_half_kappa_kHz(29) = 198
! all_detuning_kHz(29) = -176

! all_half_kappa_kHz(30) = 198
! all_detuning_kHz(30) = -360
! used for phi-c-output-4 --------------------------------------------------------------------------

! fill the (kappa/2,Delta) space of points
kD_increment = kD_limit / num_kD_slices
kD = 1
DO ii=1, num_kD_slices
    DO jj=1, num_kD_slices
        ! this will label the grid of points left-to-right, top-to-bottom in (kappa,detuning) space
        ! e.g. if num_kD_slices = 9
        ! 1, 2, 3, 4, 5, 6, 7, 8, 9
        ! 10,11,12,13,14,15,16,17,18 etc.
        all_half_kappa_kHz(kD) = jj * kD_increment
        ! crucial minus sign here for red-detuning when detuning = laser - cavity frequency
        all_detuning_kHz(kD) = -ii * kD_increment
        kD = kD + 1
    END DO
END DO

DO kD=1, num_kD_slices**2
    ! detuning of trap beam i.e. laser frequency minus cavity frequency a.k.a. Delta /kHz
    detuning = all_detuning_kHz(kD) * kHz ! /rads^-1
    kappa = 2.d0 * all_half_kappa_kHz(kD) * kHz ! /rads^-1

    ! converts integer to string
    WRITE(point_num_str, '(I0)') kD
    ! WRITE(halfk_str, '(I0)') int(all_half_kappa_kHz(kD))
    ! WRITE(delta_str, '(I0)') int(all_detuning_kHz(kD))

    ! point_dir = 'point_halfk'//trim(halfk_str)//'_D'//trim(delta_str)//'/'
    point_dir = 'point_'//trim(point_num_str)//'/'
    CALL system('mkdir '//trim(phi_c_dir)//trim(point_dir))
    OPEN(kappa_delta_point, file=trim(phi_c_dir)//trim(point_dir)//'kappa_delta_point.dat', status='unknown')
    WRITE(kappa_delta_point, *) kappa, detuning

    ! OPEN(kD * 10, file=trim(phi_c_dir)//trim(point_dir)//'QLT_spectra.dat', status='unknown')
    ! WRITE(kD * 10, *) &
    !     'kHz ', &
    !     ' S_XX S_YY S_ZZ', &
    !     ' Re(S_XY) Im(S_XY)', &
    !     ' Re(S_YX) Im(S_YX)', &
    !     ' Re(S_corr) Im(S_corr)'

    ! trim(experimental_positions_dir) as it could be '' at this point
    OPEN(omega_QLT, file=trim(phi_c_dir)//trim(point_dir)//'omega_QLT.dat', status='unknown')
    OPEN(omega_CS, file=trim(phi_c_dir)//trim(point_dir)//'omega_CS.dat', status='unknown')
    OPEN(omega_OS, file=trim(phi_c_dir)//trim(point_dir)//'omega_OS.dat', status='unknown')
    OPEN(omega_CS_plus_OS, file=trim(phi_c_dir)//trim(point_dir)//'omega_CS_plus_OS.dat', status='unknown')

    OPEN(deviation_0, file=trim(phi_c_dir)//trim(point_dir)//'deviation_0.dat', status='unknown')
    OPEN(deviation_CS, file=trim(phi_c_dir)//trim(point_dir)//'deviation_CS.dat', status='unknown')
    OPEN(deviation_OS, file=trim(phi_c_dir)//trim(point_dir)//'deviation_OS.dat', status='unknown')
    OPEN(deviation_CS_plus_OS, file=trim(phi_c_dir)//trim(point_dir)//'deviation_CS_plus_OS.dat', status='unknown')

    OPEN(OM_couplings, file=trim(phi_c_dir)//trim(point_dir)//'OM_couplings.dat', status='unknown')
    OPEN(mech_couplings, file=trim(phi_c_dir)//trim(point_dir)//'mech_couplings.dat', status='unknown')

    WRITE(omega_QLT, *) 'y0/lambda omega_x_QLT omega_y_QLT omega_z_QLT'
    WRITE(omega_CS, *) 'y0/lambda omega_x_CS omega_y_CS omega_z_CS'
    WRITE(omega_OS, *) 'y0/lambda omega_x_OS omega_y_OS omega_z_OS'
    WRITE(omega_CS_plus_OS, *) 'y0/lambda omega_x_CS_plus_OS omega_y_CS_plus_OS omega_z_CS_plus_OS'

    WRITE(deviation_0, *) 'y0/lambda deviation_0_X deviation_0_Y deviation_0_Z'
    WRITE(deviation_CS, *) 'y0/lambda deviation_CS_X deviation_CS_Y deviation_CS_Z'
    WRITE(deviation_OS, *) 'y0/lambda deviation_OS_X deviation_OS_Y deviation_OS_Z'
    WRITE(deviation_CS_plus_OS, *) 'y0/lambda deviation_CS_plus_OS_X deviation_CS_plus_OS_Y deviation_CS_plus_OS_Z'

    WRITE(OM_couplings, *) 'y0/lambda g_xY g_yY g_zP'
    WRITE(mech_couplings, *) 'y0/lambda gxy gyz gzx'

    ! loop over the y0 equilibrium position
    DO ii=1, num_y0_samples
        ! increase y0 in increments of num_y0_samples segments of a quarter wavelength (0.25 * lambda)
        lambda_coeff = ii * node / num_y0_samples

        IF (num_y0_samples == 1) THEN
            ! if ignoring y0 sampling, just run the loop once
            lambda_coeff = y0_to_record_spectra_at
        END IF

        y0_value = lambda_coeff * tweezer_wavelength
        cavity_photons = (E_drive / hbar)**2 * cos(linewidth_k * y0_value)**2 / ((0.5d0 * kappa)**2 + detuning**2)
        IF (loop_debug == 1) THEN
            WRITE(stdout, *)
            WRITE(stdout, *) '#########################################################################'
            WRITE(stdout, *) '(kappa/2, detuning) loop', kD, 'of', num_kD_slices**2
            WRITE(stdout, *) '(kappa/2, detuning) = (', 0.5d0 * kappa / kHz, detuning / kHz, ') kHz'
            WRITE(stdout, *) 'y0 loop', ii, 'of', num_y0_samples
            WRITE(stdout, *) '-------------------------------------------------------------------------'
            WRITE(stdout, *) 'y0 /lambda = ', lambda_coeff
            WRITE(stdout, *) 'Number of photons in cavity', cavity_photons
        END IF

        CALL CALCULATE_EQUILIBRIUM_PARAMETERS(&
            y0_value, &
            linewidth_k, theta_tweezer_radians, &
            G_matrix, &
            Rayleigh_range, bead_mass, &
            E_drive, &
            ! use kappa_exp here, not kappa_calc
            detuning, kappa, Gamma_M, &
            omega_x_0, omega_y_0, omega_z_0, &
            omega_x_CS, omega_y_CS, omega_z_CS, &
            omega_x_OS, omega_y_OS, omega_z_OS, &
            omega_x_CS_plus_OS, omega_y_CS_plus_OS, omega_z_CS_plus_OS, &
            num_phonons_X, num_phonons_Y, num_phonons_Z, &
            phonons_from_cooling_formula_array, &
            equilibrium_debug, cooling_debug, stdout)

        ! shot noise
        num_photons = 0.d0

        ! open loop over frequency for noise spectra
        ! e.g. S_XX(omega) = FT(autocorrelation<X(t).X^dag(t + tau)>)
        S_homodyne_array = 0.d0

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
                ! OPT corrrection is a dynamical effect which calculating the PSDs via QLT should produce
                ! so just use the CS correction as the input for the QLT PSDs
                omega_x_CS, omega_y_CS, omega_z_CS, &
                S_XX_value, S_YY_value, S_ZZ_value, &
                ! cross-correlation spectra (can be complex)
                S_XY_value, S_YX_value, S_corr_value, &
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
                omega_x_CS, omega_y_CS, omega_z_CS, &
                S_XX_value, S_YY_value, S_ZZ_value, &
                S_XY_value, S_YX_value, S_corr_value, &
                S_homodyne_value, S_heterodyne_value)

            ! update
            S_homodyne_negative = 0.5d0 * (S_homodyne_negative + S_homodyne_value)
            omega_kHz = omega_value / kHz

            ! WRITE(kD * 10, 1015) &
            !     omega_kHz, &
            !     S_XX_value, S_YY_value, S_ZZ_value, &
            !     real(S_XY_value), aimag(S_XY_value), &
            !     real(S_YX_value), aimag(S_YX_value), &
            !     real(S_corr_value), aimag(S_corr_value)

            S_XX_array(jj) = S_XX_value
            S_YY_array(jj) = S_YY_value
            S_ZZ_array(jj) = S_ZZ_value
            S_heterodyne_array(jj) = S_heterodyne_value
            ! to find optimal squeezing
            S_homodyne_array(jj) = S_homodyne_negative

        ! end of loop over omega samples
        END DO

        WRITE(omega_CS, 1015) &
            lambda_coeff, &
            omega_x_CS / kHz, &
            omega_y_CS / kHz, &
            omega_z_CS / kHz

        WRITE(omega_OS, 1015) &
            lambda_coeff, &
            omega_x_OS / kHz, &
            omega_y_OS / kHz, &
            omega_z_OS / kHz

        WRITE(omega_CS_plus_OS, 1015) &
            lambda_coeff, &
            omega_x_CS_plus_OS / kHz, &
            omega_y_CS_plus_OS / kHz, &
            omega_z_CS_plus_OS / kHz

        WRITE(omega_QLT, 1015) &
            lambda_coeff, &
            omega_array(maxloc(S_XX_array)) / kHz, &
            omega_array(maxloc(S_YY_array)) / kHz, &
            omega_array(maxloc(S_ZZ_array)) / kHz

        WRITE(deviation_0, 1015) &
            lambda_coeff, &
            (omega_array(maxloc(S_XX_array)) - omega_x_0) / kHz , &
            (omega_array(maxloc(S_YY_array)) - omega_y_0) / kHz , &
            (omega_array(maxloc(S_ZZ_array)) - omega_z_0) / kHz

        WRITE(deviation_CS, 1015) &
            lambda_coeff, &
            (omega_array(maxloc(S_XX_array)) - omega_x_CS) / kHz , &
            (omega_array(maxloc(S_YY_array)) - omega_y_CS) / kHz , &
            (omega_array(maxloc(S_ZZ_array)) - omega_z_CS) / kHz

        WRITE(deviation_OS, 1015) &
            lambda_coeff, &
            (omega_array(maxloc(S_XX_array)) - omega_x_OS) / kHz , &
            (omega_array(maxloc(S_YY_array)) - omega_y_OS) / kHz , &
            (omega_array(maxloc(S_ZZ_array)) - omega_z_OS) / kHz

        ! NOTE: this should be ~zero, as the PSD peaks should simulate the dynamical effect of the optical spring
        WRITE(deviation_CS_plus_OS, 1015) &
            lambda_coeff, &
            (omega_array(maxloc(S_XX_array)) - omega_x_CS_plus_OS) / kHz , &
            (omega_array(maxloc(S_YY_array)) - omega_y_CS_plus_OS) / kHz , &
            (omega_array(maxloc(S_ZZ_array)) - omega_z_CS_plus_OS) / kHz

        WRITE(OM_couplings, 1015) &
            lambda_coeff, &
            G_matrix(1) / kHz, &
            G_matrix(2) / kHz, &
            G_matrix(3) / kHz

        WRITE(mech_couplings, 1015) &
            lambda_coeff, &
            G_matrix(4) / kHz, &
            G_matrix(5) / kHz, &
            G_matrix(6) / kHz

    ! end of loop over y0 samples
    END DO

    ! NOTE: not sure why this was needed now, but it addressed some issue to do with two different files using the same number
    ! CLOSE(kD * 10)
    CLOSE(kappa_delta_point)
    CLOSE(omega_CS)
    CLOSE(omega_OS)
    CLOSE(omega_CS_plus_OS)
    CLOSE(omega_QLT)
    CLOSE(deviation_0)
    CLOSE(deviation_CS)
    CLOSE(deviation_OS)
    CLOSE(deviation_CS_plus_OS)
    CLOSE(OM_couplings)
    CLOSE(mech_couplings)

! end of loop over (kappa,detuning) points
END DO

WRITE(stdout, *) '#########################################################################'
WRITE(stdout, *) 'KEY SIMULATION PARAMETERS'
WRITE(stdout, *) '-------------------------------------------------------------------------'
WRITE(stdout, *) 'omega/2π range', min_omega, 'to', max_omega, 'kHz'
WRITE(stdout, *) 'omega/2π resolution', omega_increment, 'kHz'
WRITE(stdout, *) '(kappa/2, detuning)/2π resolution', kD_limit/num_kD_slices, 'kHz'
WRITE(stdout, *) 'phi/2π resolution', 0.25d0/num_y0_samples


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
    S_XY_value, S_YX_value, S_corr_value, &
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
    double complex::S_XY_value, S_YX_value, S_corr_value
    double precision::S_homodyne_value, S_heterodyne_value

    double precision::pi
    double precision::G_average
    double precision::omega_heterodyne, omega_addition
    double complex::bigGxy
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

    S_corr_value = 0.5 * (S_XY_value + S_YX_value)

    IF (equation_choice == 2) THEN
        bigGxy = a_hat_vector(1, 4)

        ! bigGxy = cavity-mediated + direct
        ! bigGxy = imag_num * eta_C * g_xY * g_yY + gxy
        ! gxy/(g_xY*g_yY) = B cot(phi)**2
        ! B = 2*detuning/(detuning**2 + half_kappa**2)

        ! So:
        ! bigGxy = g_xY * g_yY * (imag_num * eta_C + B * cot(phi)))
        ! simplified version of cross-correlation
        S_corr_value = (bigGxy / (omega_y - omega_x)) * (S_XX_value - S_YY_value)
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

    ! coeff of X-Y coupling - Combines direct and indirect paths
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
    linewidth_k, theta_tweezer_radians, &
    Rayleigh_range, bead_mass, polarisability, &
    E_tweezer, E_cavity, E_drive, &
    Gamma_M, &
    omega_0_prefactor, &
    omega_x_0, omega_y_0, omega_z_0, &
    bead_debug, stdout)
    ! """
    ! subroutine below is provided by user
    ! and calculates relevant parameters
    ! """
    IMPLICIT NONE

    integer::N_total
    double precision::hbar, speed_of_light, kB, vacuum_permittivity
    double precision::Free_Spectral_Range, bath_temperature, air_pressure, air_speed
    double precision::bead_diameter, bead_density, bead_permittivity
    double precision::tweezer_input_power, beam_waist_X, beam_waist_Y, theta_tweezer_over_pi
    double precision::cavity_waist, cavity_length, Finesse
    double precision::half_kappa_exp_kHz, tweezer_wavelength
    INCLUDE 'CSCAVITY.h'

    double precision::linewidth_k, theta_tweezer_radians
    double precision::Rayleigh_range, bead_mass, polarisability
    double precision::E_tweezer, E_cavity, E_drive
    double precision::Gamma_M
    double precision::omega_0_prefactor
    double precision::omega_x_0, omega_y_0, omega_z_0
    integer::bead_debug, stdout

    double precision::pi, pi2, kHz
    double precision::bead_radius
    double precision::omega_cavity, cavity_volume, amplitude

    ! zero eq. initial values
    pi = dacos(-1.d0)
    pi2 = 2.d0 * pi
    kHz = pi2 * 1.d3
    bead_radius = 0.5d0 * bead_diameter

    Rayleigh_range = 0.5d0 * linewidth_k * beam_waist_X * beam_waist_Y
    bead_mass = bead_density * 4.*pi/3. * bead_radius**3
    polarisability = 4.* pi * vacuum_permittivity * (bead_permittivity - 1.) / (bead_permittivity + 2.) * bead_radius**3

    omega_cavity = speed_of_light * linewidth_k
    ! add a factor of 4 here. Not in GALAX1-5 routines!!!
    cavity_volume = cavity_length * pi * cavity_waist**2 / 4.d0

    ! Depth of cavity field. Weak and unimportant for CS case
    amplitude = omega_cavity * polarisability / 2. / cavity_volume / vacuum_permittivity

    E_tweezer = sqrt(4. * tweezer_input_power / (beam_waist_X * beam_waist_Y * pi * speed_of_light * vacuum_permittivity))
    E_cavity = sqrt(hbar * omega_cavity / (2. * cavity_volume * vacuum_permittivity))
    ! use negative E_drive, in units of hbar
    ! electric field of the combined cavity and tweezer light
    E_drive = -0.5d0 * polarisability * E_tweezer * E_cavity * sin(theta_tweezer_radians)

    ! take usual expression e.g. Levitated review by Li Geraci etc
    ! Appendix E of https://arxiv.org/pdf/1308.4503.pdf
    ! 1 bar = 10^5 pascal; air_pressure /mbar = 10^2 Pascal
    Gamma_M = (16. / pi) * ((air_pressure * 1.d2) / (bead_radius * bead_density * air_speed))

    omega_0_prefactor = polarisability * E_tweezer**2 / bead_mass
    omega_x_0 = sqrt(omega_0_prefactor / beam_waist_X**2)
    omega_y_0 = sqrt(omega_0_prefactor / beam_waist_Y**2)
    omega_z_0 = sqrt(omega_0_prefactor / (2.d0 * Rayleigh_range**2))

    IF (bead_debug == 1) THEN
        WRITE(stdout, *)
        WRITE(stdout, *) '---------------------------------------------'
        WRITE(stdout, *) 'BEAD PARAMETERS'
        WRITE(stdout, *) '---------------------------------------------'
        WRITE(stdout, *) 'Rayleigh_range /m = ', Rayleigh_range
        WRITE(stdout, *) 'bead_mass /kg = ', bead_mass
        WRITE(stdout, *) 'bead polarisability / Farad.m^2 = ', polarisability
        WRITE(stdout, *)
        WRITE(stdout, *) 'Cavity trap, amplitude/2π (zeroth shift) - both at node /Hz = '
        WRITE(stdout, *) amplitude / pi2,  amplitude * cos(0.25 * pi2)**2
        WRITE(stdout, *)
        WRITE(stdout, *) 'E_tweezer /Vm^-1 = ', E_tweezer
        WRITE(stdout, *) 'E_cavity /Vm^-1 = ', E_cavity
        WRITE(stdout, *) 'E_drive /Fm^-2.(Vm^-1)^2 = ', E_drive
        WRITE(stdout, *) 'E_drive/hbar = ', E_drive / hbar
        WRITE(stdout, *)
        WRITE(stdout, *) 'Mechanical damping: Gamma_M /Hz = ', Gamma_M
        WRITE(stdout, *)
        WRITE(stdout, *) 'MECHANICAL FREQUENCIES IN THE TWEEZER FRAME'
        WRITE(stdout, *) 'omega_x_0 /kHz = ', omega_x_0 / kHz
        WRITE(stdout, *) 'omega_y_0 /kHz = ', omega_y_0 / kHz
        WRITE(stdout, *) 'omega_z_0 /kHz = ', omega_z_0 / kHz
    END IF

    RETURN
END


SUBROUTINE CALCULATE_EQUILIBRIUM_PARAMETERS(&
    y0_value, &
    linewidth_k, theta_tweezer_radians, &
    G_matrix, &
    Rayleigh_range, bead_mass, &
    E_drive, &
    detuning, kappa, Gamma_M, &
    omega_x_0, omega_y_0, omega_z_0, &
    omega_x_CS, omega_y_CS, omega_z_CS, &
    omega_x_OS, omega_y_OS, omega_z_OS, &
    omega_x_CS_plus_OS, omega_y_CS_plus_OS, omega_z_CS_plus_OS, &
    num_phonons_X, num_phonons_Y, num_phonons_Z, &
    phonons_from_cooling_formula_array, &
    equilibrium_debug, cooling_debug, stdout)
    ! """
    ! subroutine below obtains the optomechanical parameters
    ! """
    IMPLICIT NONE

    integer::N_total
    double precision::hbar, speed_of_light, kB, vacuum_permittivity
    double precision::Free_Spectral_Range, bath_temperature, air_pressure, air_speed
    double precision::bead_diameter, bead_density, bead_permittivity
    double precision::tweezer_input_power, beam_waist_X, beam_waist_Y, theta_tweezer_over_pi
    double precision::cavity_waist, cavity_length, Finesse
    double precision::half_kappa_exp_kHz, tweezer_wavelength
    INCLUDE 'CSCAVITY.h'

    double precision::y0_value
    double precision::linewidth_k, theta_tweezer_radians
    double precision::G_matrix(6)
    double precision::Rayleigh_range, bead_mass
    double precision::E_drive
    double precision::detuning, kappa, Gamma_M
    double precision::omega_x_0, omega_y_0, omega_z_0
    double precision::omega_x_CS, omega_y_CS, omega_z_CS
    double precision::omega_x_OS, omega_y_OS, omega_z_OS
    double precision::omega_x_CS_plus_OS, omega_y_CS_plus_OS, omega_z_CS_plus_OS
    double precision::num_phonons_X, num_phonons_Y, num_phonons_Z
    ! analytic equilibrium phonon numbers
    double precision::phonons_from_cooling_formula_array(3)
    integer::equilibrium_debug, cooling_debug, stdout

    double precision::pi, pi2, kHz
    double precision::phi, half_kappa
    double precision::mean_cavity_field_prefactor, mean_cavity_field_real, mean_cavity_field_imag
    double precision::omega_CS_prefactor, dx_CS, dy_CS, dz_CS
    double precision::X_zpf, Y_zpf, Z_zpf
    double precision::k_X_zpf, k_Y_zpf, k_Z_zpf
    double precision::g_xY, g_yY, g_zP
    double precision::gxy, gzx, gyz
    double precision::sum_x, diff_x, suscept_x, dx_OS
    double precision::sum_y, diff_y, suscept_y, dy_OS
    double precision::sum_z, diff_z, suscept_z, dz_OS
    double precision::cooling_rate_x, cooling_rate_y, cooling_rate_z
    double precision::phonons_from_cooling_formula_X, phonons_from_cooling_formula_Y, phonons_from_cooling_formula_Z

    pi = dacos(-1.d0)
    pi2 = 2.d0 * pi
    kHz = pi2 * 1.d3
    half_kappa = 0.5d0 * kappa

    phi = linewidth_k * y0_value
    ! the mean cavity field a.k.a. alpha_bar
    mean_cavity_field_prefactor = (E_drive / hbar) * cos(phi) / (half_kappa**2 + detuning**2)
    mean_cavity_field_real = mean_cavity_field_prefactor * detuning
    mean_cavity_field_imag = mean_cavity_field_prefactor * (-half_kappa)

    ! Coherent scattering correction to frequency squared
    omega_CS_prefactor = -2. * E_drive / bead_mass * mean_cavity_field_real * cos(phi)
    dx_CS = omega_CS_prefactor * (linewidth_k * sin(theta_tweezer_radians))**2
    dy_CS = omega_CS_prefactor * (linewidth_k * cos(theta_tweezer_radians))**2
    dz_CS = omega_CS_prefactor * (linewidth_k - 1.d0 / Rayleigh_range)**2
    omega_x_CS = sqrt(omega_x_0**2 + dx_CS)
    omega_y_CS = sqrt(omega_y_0**2 + dy_CS)
    omega_z_CS = sqrt(omega_z_0**2 + dz_CS)

    ! zero point fluctuations
    X_zpf = sqrt(hbar / (2.d0 * bead_mass * omega_x_CS))
    Y_zpf = sqrt(hbar / (2.d0 * bead_mass * omega_y_CS))
    Z_zpf = sqrt(hbar / (2.d0 * bead_mass * omega_z_CS))

    ! define these variables for shorthand convenience
    k_X_zpf = linewidth_k * X_zpf
    k_Y_zpf = linewidth_k * Y_zpf
    k_Z_zpf = (linewidth_k - 1.d0 / Rayleigh_range) * Z_zpf

    ! optomechanical couplings
    g_xY = (E_drive / hbar) * k_X_zpf * sin(phi) * sin(theta_tweezer_radians)
    g_yY = (E_drive / hbar) * k_Y_zpf * sin(phi) * cos(theta_tweezer_radians)
    g_zP = -(E_drive / hbar) * k_Z_zpf * cos(phi)

    gxy = (E_drive / hbar) * k_X_zpf * k_Y_zpf * mean_cavity_field_real * cos(phi) * sin(2. * theta_tweezer_radians)
    gyz = 2.d0 * (E_drive / hbar) * k_Y_zpf * k_Z_zpf * mean_cavity_field_imag * sin(phi) * cos(theta_tweezer_radians)
    gzx = 2.d0 * (E_drive / hbar) * k_Z_zpf * k_X_zpf * mean_cavity_field_imag * sin(phi) * sin(theta_tweezer_radians)

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
        WRITE(stdout, *) 'OPTOMECHANICAL COUPLINGS'
        WRITE(stdout, *) 'g_xY /kHz', g_xY / kHz
        WRITE(stdout, *) 'g_yY /kHz', g_yY / kHz
        WRITE(stdout, *) 'g_zP /kHz', g_zP / kHz
        WRITE(stdout, *)
        WRITE(stdout, *) 'MECHANICAL COUPLINGS'
        WRITE(stdout, *) 'gxy /kHz', gxy / kHz
        WRITE(stdout, *) 'gyz /kHz', gyz / kHz
        WRITE(stdout, *) 'gzx /kHz', gzx / kHz
        WRITE(stdout, *)
        WRITE(stdout, *) 'COHERENT SCATTERING CORRECTED'
        WRITE(stdout, *) 'omega_x_CS /kHz = ', omega_x_CS / kHz
        WRITE(stdout, *) 'omega_y_CS /kHz = ', omega_y_CS / kHz
        WRITE(stdout, *) 'omega_z_CS /kHz = ', omega_z_CS / kHz
    END IF

    ! Add optical spring correction to the frequency squared
    ! compare this value to the position of the PSD peaks, as this is a dynamical effect
    sum_x = omega_x_CS + detuning
    diff_x = omega_x_CS - detuning
    suscept_x = (sum_x / (sum_x**2 + half_kappa**2)) - (diff_x / (diff_x**2 + half_kappa**2))
    dx_OS = 2. * omega_x_CS * g_xY**2 * suscept_x

    sum_y = omega_y_CS + detuning
    diff_y = omega_y_CS - detuning
    suscept_y = (sum_y / (sum_y**2 + half_kappa**2)) - (diff_y / (diff_y**2 + half_kappa**2))
    dy_OS = 2. * omega_y_CS * g_yY**2 * suscept_y

    sum_z = omega_z_CS + detuning
    diff_z = omega_z_CS - detuning
    suscept_z = (sum_z / (sum_z**2 + half_kappa**2)) - (diff_z / (diff_z**2 + half_kappa**2))
    dz_OS = 2. * omega_z_CS * g_zP**2 * suscept_z

    omega_x_OS = sqrt(omega_x_0**2 + dx_OS)
    omega_y_OS = sqrt(omega_y_0**2 + dy_OS)
    omega_z_OS = sqrt(omega_z_0**2 + dz_OS)

    omega_x_CS_plus_OS = sqrt(omega_x_0**2 + dx_CS + dx_OS)
    omega_y_CS_plus_OS = sqrt(omega_y_0**2 + dy_CS + dy_OS)
    omega_z_CS_plus_OS = sqrt(omega_z_0**2 + dz_CS + dz_OS)

    IF (equilibrium_debug == 1) THEN
        WRITE(stdout, *)
        WRITE(stdout, *) 'OPTICAL SPRING CORRECTED'
        WRITE(stdout, *) 'omega_x_OS /kHz = ', omega_x_OS / kHz
        WRITE(stdout, *) 'omega_y_OS /kHz = ', omega_y_OS / kHz
        WRITE(stdout, *) 'omega_z_OS /kHz = ', omega_z_OS / kHz
        WRITE(stdout, *)
        WRITE(stdout, *) 'COHERENT SCATTERING + OPTICAL SPRING CORRECTED'
        WRITE(stdout, *) 'omega_x_CS_plus_OS /kHz = ', omega_x_CS_plus_OS / kHz
        WRITE(stdout, *) 'omega_y_CS_plus_OS /kHz = ', omega_y_CS_plus_OS / kHz
        WRITE(stdout, *) 'omega_z_CS_plus_OS /kHz = ', omega_z_CS_plus_OS / kHz
    END IF

    CALL CALCULATE_COOLING_RATES(&
        detuning, kappa, &
        g_xY, g_yY, g_zP, &
        omega_x_CS, omega_y_CS, omega_z_CS, &
        cooling_rate_x, cooling_rate_y, cooling_rate_z)

        ! X
        num_phonons_X = kB * bath_temperature / hbar / omega_x_CS
        phonons_from_cooling_formula_X = num_phonons_X * Gamma_M / (abs(cooling_rate_x) + Gamma_M)
        phonons_from_cooling_formula_array(1) = abs(phonons_from_cooling_formula_X)

        ! Y
        num_phonons_Y = kB * bath_temperature / hbar / omega_y_CS
        phonons_from_cooling_formula_Y = num_phonons_Y * Gamma_M / (abs(cooling_rate_y) + Gamma_M)
        phonons_from_cooling_formula_array(2) = abs(phonons_from_cooling_formula_Y)

        ! Z
        num_phonons_Z = kB * bath_temperature / hbar / omega_z_CS
        phonons_from_cooling_formula_Z = num_phonons_Z * Gamma_M / (abs(cooling_rate_z) + Gamma_M)
        phonons_from_cooling_formula_array(3) = abs(phonons_from_cooling_formula_Z)

    IF (cooling_debug == 1) THEN
        WRITE(stdout, *)
        WRITE(stdout, *) '---------------------------------------------'
        WRITE(stdout, *) 'COOLING RATES'
        WRITE(stdout, *) '---------------------------------------------'
        WRITE(stdout, *) 'cooling_rate_x = ', cooling_rate_x
        WRITE(stdout, *) 'X phonons: at room bath_temperature; at equilibrium'
        WRITE(stdout, 100) num_phonons_X, abs(phonons_from_cooling_formula_X)
        WRITE(stdout, *) 'X temperature: at room bath_temperature; at equilibrium'
        WRITE(stdout, *) bath_temperature, bath_temperature * Gamma_M / (Gamma_M + abs(cooling_rate_x))
        WRITE(stdout, *)
        WRITE(stdout, *) 'cooling_rate_y = ', cooling_rate_y
        WRITE(stdout, *) 'Y phonons: at room bath_temperature; at equilibrium'
        WRITE(stdout, 100) num_phonons_Y, abs(phonons_from_cooling_formula_Y)
        WRITE(stdout, *) 'Y temperature: at room bath_temperature; at equilibrium'
        WRITE(stdout, *) bath_temperature, bath_temperature * Gamma_M / (Gamma_M + abs(cooling_rate_y))
        WRITE(stdout, *)
        WRITE(stdout, *) 'cooling_rate_z = ', cooling_rate_z
        WRITE(stdout, *) 'Gamma_M', Gamma_M
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
    ! This is the Imaginary part of the self energy!
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


SUBROUTINE CALCULATE_CANCELLATION_POINT(detuning, kappa, omega, phi_c)
    ! """
    ! calculates simplified approximation for "phi_c"
    ! i.e. the trapping postion at which the cancellation of the orientation and frequency shifts cancel to zero
    ! """

    IMPLICIT NONE
    double precision::detuning, kappa, omega, phi_c
    double precision::half_kappa, a, b

    ! phi_c can be related to the cavity and optical susceptibilities, chi_c and eta_c
    ! chi_c = (-i(omega + detuning) + kappa/2)^-1
    ! eta_c = chi_c(omega) - chi_c*(-omega)
    ! eta_c = 2i.detuning / z_c*
    ! eta_c = 2i.detuning * (z_c / |z_c|^2)
    ! where
    ! z_c = a + ib
    ! z_c* = a - ib

    half_kappa = 0.5d0 * kappa
    a = half_kappa**2 + detuning**2 - omega**2
    b = kappa * omega

    phi_c = atan(sqrt(&
        (a**2 + b**2) / (a*(a + omega**2)) &
    ))

    RETURN
END
