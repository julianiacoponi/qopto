! This code has 1 beam and 1 dimension
! ****************************************************************************************************
! Calculates quantum noise spectra.
! Uses real mean fields
! Works out position and momentum optical quadrature
! ****************************************************************************************************


IMPLICIT NONE
integer::ii
integer::N_total, NPERIOD, num_points
double precision::R0, RHO, EPSR, EPSI0, C, hbar, BOLTZ, TEMP, Gravity
double precision::WK, waist_radius, XL, Finesse, air_pressure, WX, WY, DelFSR
double precision::tweezer_input_power1, detuning1, X0, Y0, Z0, Theta0, theta_homodyne

! parameter file with input values
INCLUDE 'CSCAVITY-1D.h'
PARAMETER(theta_homodyne=0.d0, num_points=10000)
double precision::omega_x

! calculated  parameters
! half_Gamma_M = mechanical damping due to gas,  also noise from fluctuation dissipation
double precision::Delta, half_kappa, kappa, theta_homodyne_pi, half_Gamma_M
double precision::max_omega, step, g_coupling

double precision::pi, pi2, S_XX_value, S_homodyne, A_homodyne, S_heterodyne
double precision::omega_sweep  ! loop variable "swept" over the num_points for data collection
double precision::nbar  ! thermal bath phonon occupancy
double precision::num_phonons, num_photons, phonons_from_optical_noise, area_under_S_XX, oscillator_temperature

! for OPTOCOOL formula
double precision::TBATH, C_plus, C_minus, x_cooling, x_phonons_from_cooling_formula

! read parameter file with input values
double precision, DIMENSION(num_points)::S_XX_array, omega_array, S_homodyne_array

! which expressions to calculate spectra with
! either 1 = original, 2 = via back actions, 3 = manually derived
integer::method
method = 3
WRITE(6, *) 'Chosen method', method

! plot optical field amplitudes here
! TODO: understand how to fill these with data
OPEN(8, file="FTanX.dat", status="unknown")
OPEN(10, file="FTanOPT.dat", status="unknown")

pi = dacos(-1.d0)
pi2 = 2.d0 * pi

! homodyne angle
theta_homodyne_pi = theta_homodyne * pi

TBATH = 300.d0

! parameters needed for final S_XX calculation
Delta = -200000.0 * pi2 
kappa = 438032.94007709384 * pi2
omega_x = 128062.66787 * pi2
half_Gamma_M = 0.324068D-02
g_coupling = 53061.566105502839 * pi2
WRITE(6, *) 'Delta, kappa, omega_x, Gamma_M, g_coupling, kB*T/hbar'
WRITE(6, 500) Delta, kappa, omega_x, half_Gamma_M * 2, g_coupling, (BOLTZ * TBATH)/hbar

half_kappa = kappa * 0.5d0

! thermal bath occupancy
nbar = (BOLTZ * TBATH)/(hbar * omega_x)

WRITE(6, *) 'nbar (thermal bath phonon occupancy) ='
WRITE(6, 200) nbar

! shot noise
num_photons = 0.d0

! estimate x-axis phonons using optomechanical cooling formula for comparison
C_plus = ((Delta + omega_x)**2 + half_kappa**2)**(-1)
C_minus = ((Delta - omega_x)**2 + half_kappa**2)**(-1)
x_cooling = -(g_coupling**2 * kappa) * (C_plus - C_minus)
x_phonons_from_cooling_formula = nbar * 2 * half_Gamma_M/(abs(x_cooling) + 2. * half_Gamma_M)

! open loop over frequency sweep for noise spectra
! e.g. S_XX(omega) = FT(autocorrelation <X(t)X^T(t+tau)>)
max_omega = 2 * omega_x * 1.001
step = max_omega/num_points

DO 11 ii = 1, num_points

      omega_sweep = -max_omega + 2 * (ii - 1) * step

      ! store frequency for integration
      omega_array(ii) = omega_sweep

      ! THE BELOW IS NOT USED FOR S_XX
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ! ! work out Power Spectral Density (PSD) of homodyne
      ! A_homodyne = 0.d0
      
      ! ! work out PSDs for negative frequency for symmetrisation of homodyne (not important for S_XX)
      ! CALL CALCULATE_SPECTRA(&
      !       g_coupling, nbar, num_photons, &
      !       theta_homodyne_pi, Delta, half_kappa, half_Gamma_M, &
      !       omega_x, -omega_sweep, S_XX_value, &
      !       S_homodyne, S_heterodyne)

      ! ! update homodyne
      ! A_homodyne = A_homodyne + S_homodyne
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      ! work out same PSDs for positive and symmetrise homodyne
      CALL CALCULATE_SPECTRA(&
            method, &
            g_coupling, nbar, num_photons, &
            theta_homodyne_pi, Delta, half_kappa, half_Gamma_M, &
            omega_x, omega_sweep, S_XX_value, &
            S_homodyne, S_heterodyne)

      S_XX_array(ii) = S_XX_value

      ! THE BELOW IS NOT USED FOR S_XX
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ! ! update homodyne
      ! A_homodyne = 0.5 * (A_homodyne + S_homodyne)

      ! ! to find optimal squeezing
      ! S_homodyne_array(ii) = A_homodyne
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


11 END DO

! back action limit - at zero pressure, this is the residual phonon number from optical noise
! this should be the theoretical limit if Gamma_M is set to zero
phonons_from_optical_noise = ((omega_x + Delta)**2 + half_kappa**2)/(-4 * omega_x * Delta)
WRITE(6, *) 'x-axis back action limited phonons ='
WRITE(6, 200) phonons_from_optical_noise

! integrate the PSD for the (rescaled) position operator X^hat = (b + b^dagg)
! uses trapezoidal rule
CALL INTEGRATE_SXX(num_points, omega_array, S_XX_array, area_under_S_XX)
WRITE(6, *) 'Area under S_XX ='
WRITE(6, 200) area_under_S_XX

! Area phonons_from_optical_noise corresponds to 2n+1 so convert to get n               
num_phonons = 0.5d0 * (area_under_S_XX - 1.d0)

! oscillator_temperature = area_under_S_XX * hbar * omega_x/BOLTZ
! original code does the above i.e. goes straight from the area XRE (2n+1) to get TEMP
! do more precisely with num_phonons (n) here
oscillator_temperature = num_phonons * hbar * omega_x/BOLTZ

WRITE(6, *) 'x-axis phonons from optomechancial cooling formula:'
WRITE(6, 200) x_phonons_from_cooling_formula
WRITE(6, *) 'x-axis phonons from S_XX sideband asymmetry:'
WRITE(6, 200) num_phonons
WRITE(6, *) 'Mechanical oscillator`s temperature /K'
WRITE(6, 200) oscillator_temperature

! 200 FORMAT(7E14.3)
200 FORMAT(ES14.2)
500 FORMAT(ES14.5)
STOP
END

! Formerly 'ANALYT'
SUBROUTINE CALCULATE_SPECTRA(&
      method, &
      g_coupling, nbar, num_photons, &
      theta_homodyne_pi, Delta, half_kappa, half_Gamma_M, &
      omega_x, omega, S_XX_value, &
      S_homodyne, S_heterodyne)   
      ! """
      ! Generic routine for calculating the noise spectra of the trap and probe beams
      ! """
      IMPLICIT NONE
      integer::N_total, method
      PARAMETER(N_total=4)
      double precision::g_coupling, nbar, num_photons, theta_homodyne_pi, Gav
      double precision::S_XX_value, S_homodyne, S_heterodyne
      double precision::Delta, omega, omega_x, half_kappa, half_Gamma_M, omega_heterodyne, omega_addition
      double complex::chi_C, chi_C_star, chi_M, chi_M_star
      double complex::BAX(1, N_total), A1(1, N_total), A1dagg(1, N_total)

      ! My alternatives with just the 1D arrays needed
      double complex::X_hat_vector(N_total), a_vector(N_total), a_dagg_vector(N_total)

      ! first get the susceptibilities                    
      CALL CALCULATE_SUSCEPTIBILITIES(&
            omega, omega_x, &
            Delta, half_kappa, half_Gamma_M, &
            chi_C, chi_C_star, chi_M, chi_M_star)

      IF (method == 1) THEN
            ! work out noise vector for X^hat, a and a^dagg
            CALL CALCULATE_NOISE_VECTORS(&
                  g_coupling, half_Gamma_M, half_kappa, &
                  chi_C, chi_C_star, chi_M, chi_M_star, &
                  BAX, A1, A1dagg)

            ! Note: this isn't halved, implying X^hat = (b + b^dagg) i.e. with no sqrt(0.5) factor
            S_XX_value = (abs(BAX(1, 1)))**2 + (nbar + 1) * (abs(BAX(1, 3)))**2 + nbar * (abs(BAX(1, 4)))**2

      ELSE IF (method == 2) THEN
            ! My alternative using the back-action expressions
            ! NOTE: Putting these expressions in Mathematica match exactly with my manually derived ones
            CALL CALCULATE_NOISE_VECTORS_VIA_BACK_ACTIONS(&
                  g_coupling, half_Gamma_M, half_kappa, &
                  chi_C, chi_C_star, chi_M, chi_M_star, &
                  X_hat_vector, a_vector, a_dagg_vector)
                  
            S_XX_value = (abs(X_hat_vector(1)))**2 + (nbar + 1) * (abs(X_hat_vector(3)))**2 + nbar * (abs(X_hat_vector(4)))**2

      ELSE IF (method == 3) THEN
            ! My alternative using manually derived expressions
            CALL CALCULATE_NOISE_VECTORS_VIA_MANUAL_DERIVATIONS(&
                  g_coupling, half_Gamma_M, half_kappa, &
                  chi_C, chi_C_star, chi_M, chi_M_star, &
                  X_hat_vector, a_vector, a_dagg_vector)

            S_XX_value = (abs(X_hat_vector(1)))**2 + (nbar + 1) * (abs(X_hat_vector(3)))**2 + nbar * (abs(X_hat_vector(4)))**2

      END IF

      ! COMPARE ALTERNATIVES - this verified that N0X was equal to X_hat_vector
      ! IF (ALL(N0X /= X_hat_vector)) THEN
      !       WRITE(6, *) 'THIS IS AN ERROR????'
      !       WRITE(6, 500) N0X
      !       WRITE(6, *) '-------'
      !       WRITE(6, 500) X_hat_vector
      ! ELSE
      !       WRITE(6, *) 'THIS IS OKAY'
      ! END IF
      ! 500 FORMAT(ES14.5)

      ! THE BELOW IS NOT USED FOR S_XX
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ! ! work out homodyne spectra using same vectors
      ! CALL HOMODYNE(N_total, nbar, num_photons, theta_homodyne_pi, A1, A1dagg, S_homodyne)

      ! ! ! work out homodyne spectra using same vectors
      ! ! CALL HOMODYNE_2(N_total, nbar, num_photons, theta_homodyne_pi, A1, A1dagg, S_homodyne)

      ! S_homodyne = (S_homodyne - 1.d0)/abs(chi_C - chi_C_star)**2/g_coupling**2/half_kappa/2.

      ! ! work out heterodyne for positive frequency branch.
      ! !  Shift frequency by heterodyne beat.

      ! ! FOR THIS EXAMPLE JUST QUADRUPLE MECHANICAL FREQUENCY.
      ! omega_addition = omega_x * 4.d0
      ! omega_addition = 0.d0
      ! omega_heterodyne = (omega + omega_addition)

      ! ! work out susceptibilities shifted in frequency
      ! CALL CALCULATE_SUSCEPTIBILITIES(&
      !       omega_heterodyne, omega_x, &
      !       Delta, half_kappa, half_Gamma_M, &
      !       chi_C, chi_C_star, chi_M, chi_M_star)

      ! ! work out noise vector again
      ! CALL CALCULATE_NOISE_VECTORS(&
      !       g_coupling, half_Gamma_M, half_kappa, &
      !       chi_C, chi_C_star, chi_M, chi_M_star, &
      !       BAX, A1, A1dagg)

      ! CALL HETERODYNE(N_total, nbar, num_photons, theta_homodyne_pi, A1, A1dagg, S_heterodyne)
      ! ! WRITE(6, *) omega_heterodyne, S_heterodyne, (S_heterodyne)/abs(chi_C)**2/GMAT(1)**2
      ! S_heterodyne=(S_heterodyne - 1.d0)/abs(chi_C_star)**2/Gav**2/half_kappa/2.
      ! ! S_heterodyne = S_heterodyne - 1.d0
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      RETURN
END

! THE BELOW IS NOT USED FOR S_XX
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! SUBROUTINE HOMODYNE(N_total, nbar, num_photons, theta_homodyne_pi, A1, A1dagg, S_homodyne)
!       ! """
!       ! Work out homodyne spectra
!       ! """
!       IMPLICIT NONE
!       integer::ii, N_total
!       double precision:: theta_homodyne_pi, S_homodyne, nbar, num_photons
!       double complex::imag_n
!       double complex::homodyne_noise_vector(1, N_total)
!       double complex::PY_hat_vector(1, N_total), Y_hat_vector(1, N_total)
!       double complex::A1(1, N_total)
!       double complex::A1dagg(1, N_total)

!       ! i
!       imag_n = cmplx(0.d0, 1.d0)
!       DO ii = 1, N_total
!             ! optical field quadrature
!             Y_hat_vector(1, ii) = A1(1, ii) + A1dagg(1, ii)
!             ! photon momentum operator
!             PY_hat_vector(1, ii) = imag_n * (A1(1, ii) - A1dagg(1, ii))
!             ! since theta_homodyne_pi is set to 0, homodyne_noise_vector = Y_hat_vector
!             ! Therefore S_@@ = S_YY (@=theta)
!             homodyne_noise_vector(1, ii) = PY_hat_vector(1, ii) * sin(theta_homodyne_pi) + Y_hat_vector(1, ii) * cos(theta_homodyne_pi)
!       END DO

!       S_homodyne = 0.d0
!       S_homodyne = S_homodyne + abs(homodyne_noise_vector(1, 1))**2 * num_photons + abs(homodyne_noise_vector(1, 2))**2 * (num_photons + 1.d0)
!       S_homodyne = S_homodyne + abs(homodyne_noise_vector(1, 3))**2 * nbar + abs(homodyne_noise_vector(1, 4))**2 * (nbar + 1.d0)
 
!       RETURN
! END


! My own version with 1D vectors
! SUBROUTINE HOMODYNE_2(N_total, nbar, num_photons, theta_homodyne_pi, a_vector, a_dagg_vector, S_homodyne)
!       ! """
!       ! Work out homodyne spectra
!       ! """
!       IMPLICIT NONE
!       integer::ii, N_total
!       double precision:: theta_homodyne_pi, S_homodyne, nbar, num_photons
!       double complex::imag_n
!       double complex::homodyne_noise_vector(N_total), PY_hat_vector(N_total), Y_hat_vector(N_total)
!       double complex::a_vector(N_total), a_dagg_vector(N_total)

!       ! i
!       imag_n = cmplx(0.d0, 1.d0)
!       DO ii = 1, N_total
!             ! optical field quadrature
!             Y_hat_vector(ii) = a_vector(ii) + a_dagg_vector(ii)
!             ! photon momentum operator
!             PY_hat_vector(ii) = imag_n * (a_vector(ii) - a_dagg_vector(ii))
!             ! since theta_homodyne_pi is set to 0, homodyne_noise_vector = Y_hat_vector
!             ! Therefore S_@@ = S_YY (@=theta)
!             homodyne_noise_vector(ii) = PY_hat_vector(ii) * sin(theta_homodyne_pi) + Y_hat_vector(ii) * cos(theta_homodyne_pi)
!       END DO

!       S_homodyne = 0.d0
!       S_homodyne = S_homodyne + abs(homodyne_noise_vector(1, 1))**2 * num_photons + abs(homodyne_noise_vector(1, 2))**2 * (num_photons + 1.d0)
!       S_homodyne = S_homodyne + abs(homodyne_noise_vector(1, 3))**2 * nbar + abs(homodyne_noise_vector(1, 4))**2 * (nbar + 1.d0)
 
!       RETURN
! END


! SUBROUTINE HETERODYNE(N_total, nbar, num_photons, theta_homodyne_pi, A1, A1dagg, S_heterodyne)
!       IMPLICIT NONE
!       integer::ii, N_total
!       double precision::THETA, S_heterodyne, nbar, num_photons
      
!       double complex::imag_n
!       !double complex::XTHET1(1,N_total)
!       double complex:: PY_hat_vector(1, N_total), Y_hat_vector(1, N_total)
!       double complex:: A1(1, N_total)
!       double complex:: A1dagg(1, N_total), heterodyne_noise_vector(1, N_total)
      
!       ! (0 + i)
!       imag_n = cmplx(0.d0, 1.d0)
!       DO ii = 1, N_total
!             !        Y_hat_vector(1, ii) = A1(1, ii) + A1dagg(1, ii)
!             !        PY_hat_vector(1, ii) = imag_n * (A1(1, ii) - A1dagg(1, ii))
!             !        heterodyne_noise_vector(1, ii) = PY_hat_vector(1, ii) * sin(theta_homodyne_pi) + Y_hat_vector(1, ii) * cos(theta_homodyne_pi)
!             heterodyne_noise_vector(1, ii) = A1dagg(1, ii)
!       END DO
!       S_heterodyne = 0.d0
      
!       S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(1, 1))**2 * num_photons + abs(heterodyne_noise_vector(1, 2))**2 * (num_photons + 1.d0)
!       S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(1, 3))**2 * nbar + abs(heterodyne_noise_vector(1, 4))**2 * (nbar + 1.d0)
      
!       RETURN
! END

! ! My own version with 1D vectors
! SUBROUTINE HETERODYNE_2(N_total, nbar, num_photons, theta_homodyne_pi, A1, A1dagg, S_heterodyne)
!       IMPLICIT NONE
!       integer::ii, N_total
!       double precision::theta_homodyne_pi, S_heterodyne, nbar, num_photons      
!       double complex::imag_n
!       double complex::PY_hat_vector(N_total), Y_hat_vector(N_total)
!       double complex::a_vector(N_total), a_dagg_vector(N_total), heterodyne_noise_vector(N_total)
      
!       ! (0 + i)
!       imag_n = cmplx(0.d0, 1.d0)
!       DO ii = 1, N_total
!             !        Y_hat_vector(ii) = A1(ii) + A1dagg(ii)
!             !        PY_hat_vector(ii) = imag_n * (A1(ii) - A1dagg(ii))
!             !        heterodyne_noise_vector(ii) = PY_hat_vector(ii) * sin(theta_homodyne_pi) + Y_hat_vector(ii) * cos(theta_homodyne_pi)
!             heterodyne_noise_vector(ii) = A1dagg(ii)
!       END DO
!       S_heterodyne = 0.d0
      
!       S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(1))**2 * num_photons + abs(heterodyne_noise_vector(2))**2 * (num_photons + 1.d0)
!       S_heterodyne = S_heterodyne + abs(heterodyne_noise_vector(3))**2 * nbar + abs(heterodyne_noise_vector(4))**2 * (nbar + 1.d0)
      
!       RETURN
! END


! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! Formerly 'SUSCEPT'
SUBROUTINE CALCULATE_SUSCEPTIBILITIES(&
      omega, omega_x, &
      Delta, half_kappa, half_Gamma_M, &
      chi_C, chi_C_star, chi_M, chi_M_star)
      ! """
      ! Calculates the optical (chi_C, chi_C_star) and mechanical (chi_M, chi_M_star) susceptibilities
      ! """
      IMPLICIT NONE 
      double precision::omega, omega_x
      double precision::Delta, half_Gamma_M, half_kappa
      double complex::chi_C, chi_C_star, chi_M, chi_M_star
      double precision::chi_real, chi_imag, chi_denominator

      ! Optical susceptibilities
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

      ! x-axis mechanical susceptibilities
      ! chi_M
      chi_denominator = (half_Gamma_M)**2 + (omega - omega_x)**2
      chi_real =  (half_Gamma_M)/chi_denominator
      chi_imag = (omega - omega_x)/chi_denominator
      chi_M = cmplx(chi_real, chi_imag)

      ! chi_M*(-omega) - i.e. turn omega into -omega, then chi_imag into -chi_imag.
      chi_denominator = (half_Gamma_M)**2 + (-omega - omega_x)**2
      chi_real =  (half_Gamma_M)/chi_denominator
      chi_imag = (-omega - omega_x)/chi_denominator
      chi_M_star = cmplx(chi_real, -chi_imag)

      RETURN
END

! Formerly 'Avect'
SUBROUTINE CALCULATE_NOISE_VECTORS(&
      g_coupling, half_Gamma_M, half_kappa, &
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
      double precision::g_coupling, half_Gamma_M, half_kappa, sqrt_Gamma_M
      PARAMETER(N_total=4)
      double complex::chi_C, chi_C_star, chi_CBA
      double complex::chi_M, chi_M_star, chi_MBA
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

      ! THE BELOW IS NOT USED FOR S_XX
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ! ! now work out the optical trap field, a
      ! CA1 = imag_n * chi_C
      ! ! now work out the photon field, a^dagg
      ! CA1dagg = -imag_n * chi_C_star

      ! DO ii = 1, N_total
      !       A1(1, ii) = CA1 * g_coupling * BAX(1, ii)
      !       A1dagg(1, ii) = CA1dagg * g_coupling * BAX(1, ii)
      ! END DO

      ! ! add shot or incoming noise
      ! ! trap beam: add cavity-filtered contribution
      ! A1(1, 1) = A1(1, 1) + Sqrtkapp * chi_C
      ! A1dagg(1, 2) = A1dagg(1, 2) + Sqrtkapp * chi_C_star
      
      ! ! cavity output : add incoming imprecision
      ! ! work out a_out = a_in - Sqrtkapp(a)
      ! DO ii = 1, N_total
      !       A1(1, ii) = -A1(1, ii) * Sqrtkapp
      !       A1dagg(1, ii) = -A1dagg(1, ii) * Sqrtkapp
      ! END DO
      
      ! A1(1, 1) = one + A1(1, 1)
      ! A1dagg(1, 2) = one + A1dagg(1, 2)
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      
      RETURN
END

SUBROUTINE CALCULATE_NOISE_VECTORS_VIA_BACK_ACTIONS(&
      g_coupling, half_Gamma_M, half_kappa, &
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
      double precision::g_coupling, half_Gamma_M, half_kappa, sqrt_kappa, sqrt_Gamma_M
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

      ! THE BELOW IS NOT USED FOR S_XX
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ! DO ii = 1, N_total
      !       ! now work out the optical trap field, a
      !       a_vector(ii) = ig * chi_C * X_hat_vector(ii)
      !       ! now work out the photon field, a^dagg
      !       a_dagg_vector(ii) = -ig * chi_C_star * X_hat_vector(ii)
      ! END DO

      ! ! add shot or incoming noise
      ! ! trap beam: add cavity-filtered contribution
      ! a_vector(1) = a_vector(1) + sqrt_kappa * chi_C
      ! a_dagg_vector(2) = a_dagg_vector(2) + sqrt_kappa * chi_C_star
      
      ! ! cavity output: add incoming imprecision
      ! ! work out a_out = a_in - sqrt_kappa * a
      ! DO ii = 1, N_total
      !       a_vector(ii) = -a_vector(ii) * sqrt_kappa
      !       a_dagg_vector(ii) = -a_dagg_vector(ii) * sqrt_kappa
      ! END DO
      
      ! a_vector(1) = one + a_vector(1)
      ! a_dagg_vector(2) = one + a_dagg_vector(2)
      ! ! TODO: write out expressions implied by the above for a and a^dagg
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      RETURN
END


SUBROUTINE CALCULATE_NOISE_VECTORS_VIA_MANUAL_DERIVATIONS(&
      g_coupling, half_Gamma_M, half_kappa, &
      chi_C, chi_C_star, chi_M, chi_M_star, &
      X_hat_vector, a_vector, a_dagg_vector)
      ! """
      ! Work out vectors (in terms of noise operators) of:
      ! - the (rescaled) position operator X^hat and,
      ! - the optical field operators (a, a^dagg)
      ! (alternative but equivalent expressions I derived manually)
      ! """
      IMPLICIT NONE
      integer::N_total
      double precision::g_coupling, half_Gamma_M, half_kappa, sqrt_kappa, sqrt_Gamma_M
      PARAMETER(N_total=4)
      double complex::one, imag_n, ig
      double complex::chi_C, chi_C_star, chi_CBA
      double complex::chi_M, chi_M_star, chi_MBA
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

      ! optical back-action
      ! susceptibility of optical field to back-action thermal pressure fluctuations
      chi_CBA = ig * (chi_C - chi_C_star)

      ! mechanical back-action
      ! susceptibility of mechanical oscillator to back-action radiation pressure fluctuations
      chi_MBA = ig * (chi_M - chi_M_star)

      A_M = 1.d0/(1.d0 - ig * chi_M * chi_CBA)
      A_M_star = 1.d0/(1.d0 + ig * chi_M_star * chi_CBA)

      P_M = A_M * chi_M
      P_M_star = A_M_star * chi_M_star

      Q_M = ig * P_M * P_M_star * chi_CBA
      ! Q_M_star = -Q_M

      B_M = 1.d0/(1.d0 + ig * Q_M * chi_CBA)
      ! B_M_star = B_M

      sqrt_kappa = sqrt(2.d0 * half_kappa)
      sqrt_Gamma_M = sqrt(2.d0 * half_Gamma_M)

      k_M1 = B_M * P_M * sqrt_Gamma_M
      k_M2 = B_M * Q_M * sqrt_Gamma_M
      k_M3 = ig * B_M * (P_M - Q_M) * chi_C * sqrt_kappa
      k_M4 = ig * B_M * (P_M - Q_M) * chi_C_star * sqrt_kappa

      k_M1_star = B_M * P_M_star * sqrt_Gamma_M
      k_M2_star = -B_M * Q_M * sqrt_Gamma_M
      k_M3_star = -ig * B_M * (P_M_star + Q_M) * chi_C_star * sqrt_kappa
      k_M4_star = -ig * B_M * (P_M_star + Q_M) * chi_C * sqrt_kappa
      
      ! X^hat noise vector; weights of a_in, a_in^dagg, b_in, b_in^dagg
      ! f_3 = k_M3 + k_M4_star i.e. for a_in
      X_hat_vector(1) = k_M3 + k_M4_star
      ! f_4 = k_M4 + k_M3_star i.e. for a_in^dagg
      X_hat_vector(2) = k_M4 + k_M3_star
      ! f_1 = k_M1 + k_M2_star i.e. for b_in
      X_hat_vector(3) = k_M1 + k_M2_star
      ! f_2 = k_M2 + k_M1_star i.e. for b_in^dagg
      X_hat_vector(4) = k_M2 + k_M1_star

      ! THE BELOW IS NOT USED FOR S_XX
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ! DO ii = 1, N_total
      !       ! now work out the optical trap field, a
      !       a_vector(ii) = ig * chi_C * X_hat_vector(ii)
      !       ! now work out the photon field, a^dagg
      !       a_dagg_vector(ii) = -ig * chi_C_star * X_hat_vector(ii)
      ! END DO

      ! ! add shot or incoming noise
      ! ! trap beam: add cavity-filtered contribution
      ! a_vector(1) = a_vector(1) + sqrt_kappa * chi_C
      ! a_dagg_vector(2) = a_dagg_vector(2) + sqrt_kappa * chi_C_star
      
      ! ! cavity output: add incoming imprecision
      ! ! work out a_out = a_in - sqrt_kappa * a
      ! DO ii = 1, N_total
      !       a_vector(ii) = -a_vector(ii) * sqrt_kappa
      !       a_dagg_vector(ii) = -a_dagg_vector(ii) * sqrt_kappa
      ! END DO
      
      ! a_vector(1) = one + a_vector(1)
      ! a_dagg_vector(2) = one + a_dagg_vector(2)
      ! ! TODO: write out expressions implied by the above for a and a^dagg
      ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      RETURN
END


! Formerly 'NORM1'
SUBROUTINE INTEGRATE_SXX(num_points, omega_array, S_XX_array, area_under_S_XX)
      IMPLICIT NONE 
      integer::ii, num_points
      double precision::pi2, area_under_S_XX, trapeze_base
      double precision, DIMENSION(num_points)::omega_array
      double precision::S_XX_array(num_points)
     
      pi2 = 2.d0 * dacos(-1.d0)

      ! integrate the position spectrum of bead
      ! quick hack - use trapezoidal rule - improve later

      ! each trapeze base is the step size in the omega_array i.e. max_omega/num_points 
      trapeze_base = abs(omega_array(2) - omega_array(1))
      DO ii = 1, num_points - 1
            area_under_S_XX = area_under_S_XX + 0.5d0 * (S_XX_array(ii) + S_XX_array(ii + 1))
      END DO

      ! this is equal to 2n+1
      area_under_S_XX = area_under_S_XX * trapeze_base/pi2
      RETURN
END