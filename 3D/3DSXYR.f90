! This code has 1 beam plus 3 dimensions
! Calculates quantum noise spectra.
! USES real mean fields
! works out position and momentum optical quadrature

! GammaM=mechanical damping due to gas,  also noise from fluctuation dissipation.
! NB no intent(in) or INTENT(OUT) used- add these later!!


IMPLICIT NONE
integer::ll, kk, mm, it, im, ii, jj, NT, ithet, Nthet, Nx0, ix0, Ndet

integer::NTOT, NPERIOD, NPTS
double precision::R0, RHO, EPSR, EPSI0, C, hbar, BOLTZ, TEMP, Gravity
double precision::WK, waist, XL, Finesse, Press, WX, WY, DelFSR
double precision::PIN1, DET1, X0, Y0, Z0, Theta0, Thetahom, Delx0, XX0, Deldet, DET2PIini, test
double precision::Delthet, cavphotn
! parameter file with input values

INCLUDE 'CSCAVITY.h'
! PARAMETER(Nthet=1, Delthet=0.01d0, Ndet=4, Deldet=50.d3, Thetahom=0.12d0, NPTS=10000)
PARAMETER(Nthet=1, Delthet=0.01d0, Ndet=1, Deldet=100.d3, Thetahom=0.12d0, NPTS=10000)
double precision::OMx, OMy, OMz

! calculated  parameters
double precision::Det2pi, kappnano, kappin, kapp2, kappa, Polaris, thet, THEhom
double precision::XM, GAMMAM, Wkx0, epsTW, epsCAV, ZR, Ed
double precision::XMAX, DEL

double precision::pi, pi2, OMsweep, XXQM, YYQM, ZZQM, SHOM1, AHOM1, SHET1, XXQMs, YYQMs, ZZQMs
double precision::TPERIOD, AVNX, AVNY, AVNZ, AVPHOT, AVRE, OPTOCOOL, PHONONS
double precision::COOLX, COOLY, COOLZ, TEMPX, TEMPY, TEMPZ, TBATH

! read parameter file with input values
double precision::GMAT(6)
double precision::PHON(3), AV(3)
double precision, DIMENSION(NPTS)::SXXQM, SYYQM, SZZQM, OMSTOR, Shom
double complex::XI, ZOM, SXYQM, SYXQM, XYQMs, YXQMs, XYQM, YXQM

! Here plot optical field amplitudes
OPEN(8, file="FTanXYZ.dat", status="unknown")
OPEN(9, file="FTSXY.dat", status="unknown")
OPEN(10, file="FTanOPT.dat", status="unknown")
OPEN(12, file="GCOUPLE.dat", status="unknown")
OPEN(14, file="PHONS.dat", status="unknown")

pi = dacos(-1.d0)
pi2 = 2.d0 * pi
DET2PIini = DET1 * pi2
TBATH = TEMP
NT = NTOT

! OPEN LOOP OVER DEtun
! DO  10 ithet=1, Nthet
DO  10 ix0 = 1, Ndet
        DET2PI = DET2PIini + Deldet * pi * 2 * (ix0-1)
        WRITE(6, *) 'detuning', Det2pi / 2 / pi
        ! thet = theta0 + Delthet * (ithet - 1)
        thet = theta0
        ! XX0 = Wk * x0 + (ix0 - 1) * Delx0
        xx0 = Wk * x0
        Wkx0 = xx0
        thet = thet * pi
        THEhom = Thetahom * pi

        WRITE(6, *) 'theta/pi', thet/pi
        ! note BEAD routine uses x0 but not theta
        CALL BEAD_PARAMETERS(Wkx0, Polaris, epsTW, epsCAV, ZR, XM, kappin, kappnano, GAMMAM)

        CALL EQUILIBRIUM_PARAMETERS(thet, Det2pi, Wkx0, Polaris, epsTW, epsCAV, XM, ZR, kappnano, kappin, GammaM, OMX, OMY, OMZ, GMAT, PHON)
        WRITE(12, 120) thet, (abs(Gmat(ii)), ii = 1, 6)
        120 FORMAT(7E14.6)

        Ed = 0.5d0 * Polaris * epsTW * epsCAV * sin(thet)
        ! Ed = Ed / hbar
        ! cavphotn = Ed * Ed * cos(Wkx0) * cos(Wkx0) / hbar / hbar

        kappa = (kappnano + kappin)
        kapp2 = kappa * 0.5d0
        WRITE(6, *) 'kappa in main', kappa/2/pi, kapp2/2/pi
        ! cavphotn = Cavphotn / (kapp2**2 + Det2pi**2)

        ! thermal bath occupancy
        AVNX = BOLTZ * TBATH / hbar / OMX
        AVNY = BOLTZ * TBATH / hbar / OMY
        AVNZ = BOLTZ * TBATH / hbar / OMZ
        ! shot noise
        AVPHOT = 0.d0

        ! open loop over frequency for noise spectra eg S_xx(\omega)=FT(autocorrelation <X(t)X^T(t+tau)>)
        XMAX = 2 * OMX * 1.001
        DEL = XMAX / NPTS
        SXXQM = 0.d0
        SYYQM = 0.d0
        SZZQM = 0.d0
        SXYQM = 0.d0
        SYXQM = 0.d0
        SHOM = 0.d0

        DO 11 mm = 1, NPTS
                OMsweep = -XMAX + 2*(mm - 1) * DEL
                ! store frequency for integration
                OMSTOR(mm) = OMsweep
                ! work out PSD of homodyne
                AHOM1 = 0.d0
                ! symmetrised
                XXQMs = 0.d0
                YYQMs = 0.d0
                ZZQMs = 0.d0
                XYQMs = 0.d0
                YXQMs = 0.d0! ! work out  for negative frequency for symmetrisation
                CALL ANALYT(GMAT, AVNX, AVNY, AVNZ, AVPHOT, THEhom, DET2pi, Kapp2, GAMMAM, OMX, OMY, OMZ, -OMsweep, XXQM, YYQM, ZZQM, XYQM, YXQM, SHOM1, SHET1)

                ! update homodyne and symm disp.
                AHOM1 = AHOM1 + SHOM1

                ! work out same for positive and symmetrise homodyne
                CALL ANALYT(GMAT, AVNX, AVNY, AVNZ, AVPHOT, THEhom, DET2pi, Kapp2, GAMMAM, OMX, OMY, OMZ, OMsweep, XXQM, YYQM, ZZQM, XYQM, YXQM, SHOM1, SHET1)
                ! update
                AHOM1 = 0.5 * (AHOM1 + SHOM1)
                XXQMs = XXQMs + XXQM
                YYQMs = YYQMs + YYQM
                ZZQMs = ZZQMs + ZZQM
                ! cross-correlation spectra (can be complex)
                XYQMs = XYQMs + XYQM
                YXQMs = YXQMs + YXQM
                omsweep = omsweep / 2 / pi*1.d-3
                WRITE(8, 101) OMsweep, XXQMs, YYQMs, ZZQMs
                WRITE(9, 101) OMsweep, XYQMs, YXQMs
                WRITE(10, 101) OMsweep, SHET1, AHom1

                101 FORMAT(7D14.6)
                SXXQM(mm) = XXQM
                SYYQM(mm) = YYQM
                SZZQM(mm) = SHET1
                ! to find optimal squeezing
                SHOM(mm) = AHOM1
        11 END DO

        AVRE = (omx + Det2pi)**2 + kapp2**2
        AVRE = AVRE / 4 / omx / (-det2pi)

        WRITE(6, *) 'X: back action limited phonons', AVRE
        AVRE = (omy + Det2pi)**2 + kapp2**2
        AVRE = AVRE / 4 / omY / (-det2pi)
        WRITE(6, *) 'Y: back action limited phonons', AVRE
        AVRE = (omz + Det2pi)**2 + kapp2**2
        AVRE = AVRE / 4 / omz / (-det2pi)
        WRITE(6, *) 'Z: back action limited phonons', AVRE
        ! integrate and normalise the quantum  noise spectra.

        CALL NORM1(NPTS, OMX, TBATH, GAMMAM, TEMPX, SXXQM, OMSTOR, AVRE)
        ! Area AVRE corresponds to 2n+1 so convert to get n
        PHONONS = PHON(1)
        AVRE = 0.5d0 * (AVRE - 1.d0)
        AV(1) = AVRE
        WRITE(6, *) 'X phonons from formula,   from SXX FT'
        WRITE(6, 200) PHONONS, AVRE
        CALL NORM1(NPTS, OMY, TBATH, GAMMAM, TEMPY, SYYQM, OMSTOR, AVRE)
        ! Area AVRE corresponds to 2n+1 so convert to get n
        PHONONS = PHON(2)
        AVRE = 0.5d0 * (AVRE - 1.d0)
        AV(2) = AVRE
        WRITE(6, *) 'Y phonons from formula,   from SXX FT'
        WRITE(6, 200) PHONONS, AVRE

        ! PHONONS FOR HETERODYNE!!!!NOT zz!!!!
        CALL NORM1(NPTS, OMZ, TBATH, GAMMAM, TEMPZ, SZZQM, OMSTOR, AVRE)
        ! Area AVRE corresponds to 2(nx+ny)+2 so DIFFERENT CONVERSION  to get nx+ny
        PHONONS = PHON(3)
        AVRE = 0.5d0 * (AVRE - 2.d0)
        AV(3) = AVRE
        WRITE(6, *) 'Z phonons from formula,   from SXX FT'
        WRITE(6, 200) PHONONS, AVRE
        WRITE(14, 120) Det2pi/2/pi, AV(1), AV(2), AV(3), PHON(1), PHON(2), PHON(3)

! loop over theta
10 END DO

100 FORMAT(I3, 3E14.6, 1x, 2(E14.6))
200 FORMAT(7E14.6)

STOP
END


FUNCTION OPTOCOOL(G1, Det1X, KAPP2, OMEGAM, GAMMAM)
        ! """
        !  Function to evaluate optomechanical cooling formula
        ! """
        IMPLICIT NONE
        double precision::G1, Det1X, KAPP2, OMEGAM
        double precision::C1, C2, C3, C4
        double precision::OPTOCOOL, COOL1, COOL2, GAMMAM
        double precision:: hbar, BOLTZ, TBATH
        PARAMETER(BOLTZ=1.4d-23, hbar=1.05d-34, TBATH=300.0)
        COOL1=0.d0

        ! now work out opto cooling expression
        ! Trap beam
        C1 = 2. * KAPP2 * G1 * G1
        C2 = DET1X


        C3 = (C2 + omegam)**2 + kapp2**2
        C3 = 1.d0 / C3
        C4 = (C2 - omegam)**2 + kapp2**2
        C4 = 1.d0 / C4
        COOL1 = -C1 * (C3 - C4)

        OPTOCOOL = COOL1
        100 FORMAT(4D16.6)
        RETURN
END


SUBROUTINE ANALYT(GMAT, AVNX, AVNY, AVNZ, AVPHOT, THETA, DET, Kapp2, GAMMAM, OMX, OMY, OMZ, OMEGA, XXF, YYF, ZZF, XYF, YXF, SHOM1, SHET1)
        ! """
        !  Generic routine for noise spectra of trap and probe beams
        ! """
        IMPLICIT NONE
        integer::ii, m, jj, NT
        PARAMETER(NT=8)
        double precision::DET
        double precision::gamm
        double precision::pi, pi2, XXF, YYF, ZZF, PHN, SHOM1, SHET1
        double precision::OMEGA, KAPP2, GAMMAM, OMX, OMY, OMZ, omhet, oma
        double precision::CHIMSQ, XNORMSQ
        double precision::AVNX, AVNY, AVNZ, AVPHOT, THETA, Gav
        double precision::GMAT(6)
        double complex::XYF, YXF
        double complex::CHIR1, CHISTMOM1, GXY, XMUX, XMUY
        double complex::CHIMX, CHIMMOMX, CHIMY, CHIMMOMY, CHIMZ, CHIMMOMZ, RXY, RYX
        double complex::BAX(1, NT), BAY(1, NT), BAZ(1, NT), A1(1, NT), A1dagg(1, NT)


        pi = dacos(-1.d0)

        SHOM1 = 0.d0
        BAX = 0.d0
        BAY = 0.d0
        BAZ = 0.d0

        ! WORK OUT NOISE SPECTRA
        ! First do susceptibilities
        CALL SUSCEPT(OMEGA, DET, Kapp2, gammam, OMX, OMY, OMZ, CHIR1, CHISTMOM1, CHIMX, CHIMMOMX, CHIMY, CHIMMOMY, CHIMZ, CHIMMOMZ)
        ! G2NORM = G2 * G2 * (abs(CHIR2 - cos(2 * theta) * CHISTMOM2))**2

        ! work out noise vector for x,  a1 and a2
        CALL AVECT(GMAT, Gammam, kapp2, CHIR1, CHISTMOM1, CHIMX, CHIMMOMX, CHIMY, CHIMMOMY, CHIMZ, CHIMMOMZ, A1, A1dagg, BAX, BAY, BAZ)
        ! XX = sqrt(0.5) (b+b^dagg) so halve XX
        XXF = (abs(BAX(1, 1)))**2
        XXF = XXF + (AVNX + 1)*(abs(BAX(1, 3)))**2 + AVNX*(abs(BAX(1, 4)))**2 + (AVNY + 1)*abs(BAX(1, 5))**2 + AVNY*abs(BAX(1, 6))**2
        XXF = XXF + (AVNZ + 1)*(abs(BAX(1, 7)))**2 + AVNZ*(abs(BAX(1, 8)))**2

        YYF = (abs(BAY(1, 1)))**2
        YYF = YYF + (AVNX + 1)*(abs(BAY(1, 3)))**2 + AVNX*(abs(BAY(1, 4)))**2 + (AVNY + 1)*abs(BAY(1, 5))**2 + AVNY*abs(BAY(1, 6))**2
        YYF = YYF + (AVNZ + 1)*(abs(BAY(1, 7)))**2 + AVNZ*(abs(BAY(1, 8)))**2

        ZZF = (abs(BAZ(1, 1)))**2
        ZZF = ZZF + (AVNX + 1)*(abs(BAZ(1, 3)))**2 + AVNX*(abs(BAZ(1, 4)))**2 + (AVNY + 1)*abs(BAZ(1, 5))**2 + AVNY*abs(BAZ(1, 6))**2
        ZZF = ZZF + (AVNZ + 1)*(abs(BAZ(1, 7)))**2 + AVNZ*(abs(BAZ(1, 8)))**2

        ! WORK OUT CROSS SPECTRA PSD S_XY and SYX
        XYF = BAX(1, 1) * conjg(BAY(1, 1))
        XYF = XYF + (AVNX + 1)*BAX(1, 3)*conjg(BAY(1, 3)) + AVNX*BAX(1, 4) * conjg(BAY(1, 4)) + (AVNY + 1)*BAX(1, 5)*conjg(BAY(1, 5))
        XYF = XYF + AVNY*BAX(1, 6)*conjg(BAY(1, 6)) + (AVNZ + 1)*BAX(1, 7)*conjg(BAY(1, 7)) + AVNZ*BAX(1, 8)*conjg(BAY(1, 8))

        YXF = BAY(1, 1)*conjg(BAX(1, 1))
        YXF = YXF + (AVNX + 1)*BAY(1, 3)*conjg(BAX(1, 3)) + AVNX*BAY(1, 4)*conjg(BAX(1, 4)) + (AVNY + 1)*BAY(1, 5)*conjg(BAX(1, 5))
        YXF = YXF + AVNY*BAY(1, 6)*conjg(BAX(1, 6)) + (AVNZ + 1)*BAY(1, 7)*conjg(BAX(1, 7)) + AVNZ*BAY(1, 8)*conjg(BAX(1, 8))

        XYF = 0.5 * (XYF + YXF)
        ! RXY = A1(1, 1)
        ! RYX = A1(1, 2)
        ! YXF = XXF * real(RYX) + YYF * real(RXY)

        GXY = A1(1, 4)
        ! imaginary part of mu_x / M_x
        ! XMUX = -AIMAG(A1(1, 5))
        ! simplfy above
        ! XMUX = 1. / (omx - omy)
        ! imaginary part of mu_x / M_x
        ! XMUY = -AIMAG(A1(1, 6))
        ! XMUY = -1. / (omx - omy)

        ! IN THIS VERSION WE WORK OUT Eq. 18,  simplified version of cross-correlation
        YXF = (XXF - YYF) * GXY / (omy - omx)

        ! work out homodyne spectra using same vectors
        CALL HOMODYNE(NT, AVNX, AVNY, AVNZ, AVPHOT, THETA, A1, A1dagg, SHOM1)
        ! rescaling of heterodyne for Gx / Gy
        Gav = (GMAT(1) + GMAT(2)) / 2.

        SHOM1 = (SHOM1 - 1.d0) / abs(CHIR1 - CHISTMOM1)**2 / Gav**2 / kapp2 / 2.
        ! work out heterodyne for positive frequency branch.
        ! Shift frequency by heterodyne beat.
        OMA = OMx * 4.d0
        OMA = 0.d0
        OMHET = (OMEGA + OMa)

        ! work out susceptibilities shifted in frequency
        CALL SUSCEPT(OMhet, DET, Kapp2, gammam, OMX, OMY, OMZ, CHIR1, CHISTMOM1, CHIMX, CHIMMOMX, CHIMY, CHIMMOMY, CHIMZ, CHIMMOMZ)

        ! work out noise vector again
        CALL AVECT(GMAT, Gammam, kapp2, CHIR1, CHISTMOM1, CHIMX, CHIMMOMX, CHIMY, CHIMMOMY, CHIMZ, CHIMMOMZ, A1, A1dagg, BAX, BAY, BAZ)

        CALL HETERODYNE(NT, AVNX, AVNY, AVNZ, AVPHOT, THETA, A1, A1dagg, SHET1)
        ! WRITE(6, *) omhet, SHET1, (SHET1) / abs(CHIR1)**2 / GMAT(1)**2
        SHET1= (SHET1 - 1.d0) / abs(CHISTMOM1)**2 / Gav**2 / kapp2 / 2.
        ! SHET1 = SHET1 - 1.d0
        100 FORMAT(6D16.8)
        RETURN
END


SUBROUTINE HOMODYNE(NTOT, AVNX, AVNY, AVNZ, AVPHOT, THETA, A1, A1dagg, SHOM1)
        ! """
        ! Work out homodyne spectra
        ! """
        IMPLICIT NONE
        integer  ::ii, m, jj, NTOT
        double precision:: THETA, SHOM1, AVNX, AVNY, AVNZ, AVPHOT

        double complex::XI
        double complex::XTHET1(1, NTOT)
        double complex::XPM1(1, NTOT), XAM1(1, NTOT)
        double complex::A1(1, NTOT)
        double complex::A1dagg(1, NTOT)
        !  zero arrays
        XPM1 = (0.d0, 0.d0)
        XAM1 = (0.d0, 0.d0)

        XTHET1 = (0.d0, 0.d0)

        ! i
        XI = cmplx(0.d0, 1.d0)
        DO ii = 1, NTOT
                XAM1(1, ii) = A1(1, ii) + A1dagg(1, ii)
                XPM1(1, ii) = XI * (A1(1, ii) - A1dagg(1, ii))

                XTHET1(1, ii) = XPM1(1, ii) * sin(theta) + XAM1(1, ii) * cos(theta)
        END DO
        SHOM1=0.d0

        SHOM1 = SHOM1 + abs(XTHET1(1, 1))**2*AVPHOT + abs(XTHET1(1, 2))**2*(AVPHOT + 1.d0)
        SHOM1 = SHOM1 + abs(XTHET1(1, 3))**2*AVNX + abs(XTHET1(1, 4))**2*(AVNX + 1.d0)
        SHOM1 = SHOM1 + abs(XTHET1(1, 5))**2*AVNY + abs(XTHET1(1, 6))**2*(AVNY + 1.d0)
        SHOM1 = SHOM1 + abs(XTHET1(1, 7))**2*AVNZ + abs(XTHET1(1, 8))**2*(AVNZ + 1.d0)
        100   FORMAT(6D16.8)
        RETURN
END


SUBROUTINE SUSCEPT(OMEGA, DET1x, Kapp2, gamm, OMX, OMY, OMZ, CHIR1, CHISTMOM1, CHIMX, CHIMMOMX, CHIMY, CHIMMOMY, CHIMZ, CHIMMOMZ)
        ! """
        ! work out susceptibilities for noise spectra
        ! """
        IMPLICIT NONE
        double precision::DET1X
        double precision::gamm, gamm2, Kapp2
        double precision::OMEGA, OMX, OMY, OMZ
        double precision::t1, t2, ANORM, BNORM
        double complex::CHIR1, CHISTMOM1, CHIMX, CHIMMOMX, CHIMY, CHIMMOMY, CHIMZ, CHIMMOMZ

        ! gamm is actually gamma/2
        Gamm2=gamm
        ! FIELD1 susceptibilities

        ! Chi_R
        ANORM = KAPP2**2 + (omega + DET1X)**2
        t1 = kapp2 / ANORM
        t2 = (omega + DET1X) / ANORM
        CHIR1 = cmplx(t1, t2)
        ! chi_r^*(-omega)
        ANORM = KAPP2**2 + (-omega + DET1X)**2
        t1 = kapp2 / ANORM
        t2 = (-omega + DET1X) / ANORM
        CHISTMOM1 = cmplx(t1, -t2)

        ! X MECHANICAL susceptibilities
        ! chi_M X
        BNORM = (Gamm2)**2 + (omega - OMX)**2
        t1 = (Gamm2) / BNORM
        t2 = (omega - OMX) / BNORM
        CHIMX = cmplx(t1, t2)
        ! CHI_M*(-om) X
        T1 = (Gamm2)**2 + (-omega - OMX)**2
        CHIMMOMX = cmplx(Gamm2 / T1,  ( omega + OMX)/T1)

        ! Y MECHANICAL susceptibilities
        ! chi_M Y
        BNORM = (Gamm2)**2 + (omega - OMY)**2
        t1 = (Gamm2) / BNORM
        t2 = (omega - OMY) / BNORM
        CHIMY = cmplx(t1, t2)
        ! CHI_M*(-om) Y
        T1 = (Gamm2)**2 + (-omega - OMY)**2
        CHIMMOMY = cmplx(Gamm2 / T1, (omega + OMY) / T1)

        ! Z MECHANICAL susceptibilities
        ! chi_M Z
        BNORM = (Gamm2)**2 + (omega - OMZ)**2
        t1 = (Gamm2) / BNORM
        t2 = (omega - OMZ) / BNORM
        CHIMZ = cmplx(t1, t2)
        ! CHI_M*(-om) Z
        T1 = (Gamm2)**2 + (-omega - OMZ)**2
        CHIMMOMZ = cmplx(Gamm2 / T1,  (omega + OMZ) / T1)

        100 FORMAT(6D16.8)
        RETURN
END


SUBROUTINE AVECT(GMAT, Gamm, kapp2, CHIR1, CHISTMOM1, CHIMX, CHIMMOMX, CHIMY, CHIMMOMY, CHIMZ, CHIMMOMZ, A1, A1dagg, BAX, BAY, BAZ)
        ! """
        ! Work our vectors of X and optical fields but in terms of noise operators
        ! """
        IMPLICIT NONE
        integer::ii, m, jj, NTOT, NPERIOD
        double precision::ANORM, BNORM
        double precision::t1, t2, Sqrtkapp, Sqrtgamm, Gamm, Gamm2, kapp2
        double precision::GX, GY, GZ, GXY, GZX, GYZ
        DOUBLE PRECISION::GMAT(6)
        double precision::ASQ1, alp1
        PARAMETER(NTOT=8)
        double complex::CHIR1, CHISTMOM1
        double complex::CHIMX, CHIMMOMX, CHIMY, CHIMMOMY, CHIMZ, CHIMMOMZ
        double complex:: XI, ONE, Betx, BetY, BetZ, GcoefXY, GcoefZX, GcoefYZ, GcoefYX, GcoefXZ, GcoefZY
        double complex::RXY, RYX, RXZ, RZX, RYZ, RZY
        double complex::CMX, CMY, CMZ, CSUM, CNORM, eta0c, etaMpi2c, etaPpi2c, etaX, etaY, etaZ, CA1, CA1dagg
        double complex, INTENT(OUT)::BAX(1, NTOT), BAY(1, NTOT), BAZ(1, NTOT)
        double complex, INTENT(OUT)::A1(1, NTOT)
        double complex, INTENT(OUT)::A1dagg(1, NTOT)
        double complex::N0X(NTOT), N0Y(NTOT), N0Z(NTOT)
        GX = GMAT(1)
        GY = GMAT(2)
        GZ = GMAT(3)
        GXY = GMAT(4)
        GYZ = GMAT(5)
        GZX = GMAT(6)

        ! 2D
        ! actually Gamm is gamma/2
        Gamm2 = Gamm
        ! zero arrays
        BAX = cmplx(0.d0, 0.d0)
        BAY = cmplx(0.d0, 0.d0)
        BAZ = cmplx(0.d0, 0.d0)
        N0X = cmplx(0.d0, 0.d0)
        N0Y = cmplx(0.d0, 0.d0)
        N0Z = cmplx(0.d0, 0.d0)
        A1 = cmplx(0.d0, 0.d0)
        A1dagg = cmplx(0.d0, 0.d0)
        ! i
        XI = cmplx(0.d0, 1.d0)
        ONE = cmplx(1.d0, 0.0d0)

        !      SUSCEPTIBILITIES
        eta0c = CHIR1 - CHISTMOM1
        etaMpi2c = XI * (CHIR1 + CHISTMOM1)
        etaPpi2c = -XI * (CHIR1 + CHISTMOM1)
        etaX = CHIMX - CHIMMOMX
        etaY = CHIMY - CHIMMOMY
        etaZ = CHIMZ - CHIMMOMZ

        ! coeff of X-Y coupling- Combines direct and indirect paths
        GcoefXY = (GXY + Xi * eta0c * G X * GY)
        GcoefYX = (GXY + Xi * eta0c * GX * GY)
        GcoefXZ = (GZX + Xi * etaMpi2c * GZ * GX)
        GcoefZX = (GZX + Xi * etaPpi2c * GZ * GX)
        GcoefYZ = (GYZ + Xi * etaMpi2c * GY * GZ)
        GcoefZY = (GYZ + Xi * etaPpi2c * GY * GZ)

        ! NORMALIZATIONS
        CMX = 1.d0 + GX**2 * etax * eta0c
        CMY = 1.d0 + GY**2 * etay * eta0c
        CMZ = 1.d0 + GZ**2 * etaz * eta0c

        Sqrtkapp = sqrt(2.d0 * KAPP2)
        Sqrtgamm = sqrt(2.d0 * Gamm2)
        BETx = XI * Sqrtkapp * etax * GX
        BETy = XI * Sqrtkapp * etay * GY
        BETz = XI * Sqrtkapp * etaz * GZ

        ! zero-th order X noise vector; weights of a1, a1*, bx, bx*, by, by*, bz, bz*
        N0X(1) = BETX * CHIR1 / CMX
        N0X(2) = BETX * CHISTMOM1 / CMX
        N0X(3) = Sqrtgamm * CHIMX / CMX
        N0X(4) = Sqrtgamm * CHIMMOMX / CMX

        ! Zero-th order Y noise vector;weights of a1, a1*, bx, bx*, by, by*, bz, bz*
        N0Y(1) = BETY * CHIR1 / CMY
        N0Y(2) = BETY * CHISTMOM1 / CMY
        N0Y(5) = Sqrtgamm * CHIMY / CMY
        N0Y(6) = Sqrtgamm * CHIMMOMY / CMY
        ! Zero-th Z noise vector;weights of a1, a1*, bx, bx*, by, by*, bz, bz*
        N0Z(1) = -XI * BETz * CHIR1 / CMZ
        N0Z(2) = XI * BETz * CHISTMOM1 / CMZ
        N0Z(7) = Sqrtgamm * CHIMZ / CMZ
        N0Z(8) = Sqrtgamm * CHIMMOMZ / CMZ

        ! Higher order
        RXY = xi * etax * GcoefXY / CMX
        RYX = xi * etay * GcoefYX / CMY
        RXZ = xi * etax * GcoefXZ / CMX
        RZX = xi * etaz * GcoefZX / CMZ
        RYZ = xi * etay * GcoefYZ / CMY
        RZY = xi * etaz * GcoefZY / CMZ

        !
        CNORM = 1.d0 - RZX * RXZ - RZY * RYZ - RYX * RXY - RZX * RXY * RYZ - RYX * RXZ * RZY

        ! ADD 3D BACK-ACTION TERMS
        DO ii=1, NTOT
                CSUM = (1.d0 - RZY * RYZ) * N0X(ii) + (RXY + RXZ * RZY) * N0Y(ii) + (RXZ + RXY * RYZ) * N0Z(ii)
                BAX(1, ii) = BAX(1, ii) + CSUM / CNORM
                CSUM = (1.d0 - RZX * RXZ) * N0Y(ii) + (RYX + RYZ * RZX) * N0X(ii) + (RYZ + RYX * RXZ) * N0Z(ii)
                BAY(1, ii) = BAY(1, ii) + CSUM / CNORM
                CSUM = (1.d0 - RYX * RXY) * N0Z(ii) + (RZX + RZY * RYX) * N0X(ii) + (RZY + RZX * RXY) * N0Y(ii)
                ! put coefficients in BAZ
                BAZ(1, ii) = BAZ(1, ii) + CSUM / CNORM
        END DO

        ! now work out the optical trap field =a1
        CA1 = XI * CHIR1
        ! now work out the photon field a1dagger
        CA1dagg = -XI * CHISTMOM1

        DO ii=1, NTOT
                !! BUG FIX OF 29/7/2020
                A1(1, ii)=CA1*(GX*BAX(1, ii)+GY*BAY(1, ii))-GZ*CHIR1*BAZ(1, ii)
                A1dagg(1, ii)=CA1dagg*(GX*BAX(1, ii)+GY*BAY(1, ii))-GZ*CHISTMOM1*BAZ(1, ii)
        END DO

        ! add shot or incoming noise
        ! trap beam: add cavity-filtered contribution
        A1(1, 1) = A1(1, 1) + Sqrtkapp * CHIR1
        A1dagg(1, 2) = A1dagg(1, 2) + Sqrtkapp * CHISTMOM1

        ! cavity output : add incoming imprecision
        ! work out a_out=a_in-Sqrtkapp(a)
        DO ii=1, NTOT
                A1(1, ii)=-A1(1, ii)*Sqrtkapp
                A1dagg(1, ii)=-A1dagg(1, ii)*Sqrtkapp
        END DO

        A1(1, 1) = ONE + A1(1, 1)
        A1dagg(1, 2) = ONE + A1dagg(1, 2)
        A1(1, 1) = RXY
        A1(1, 2) = RYX
        A1(1, 3) = GcoefXY
        A1(1, 4) = GcoefYX
        A1(1, 5) = etax/CMX
        A1(1, 6) = etay/CMY
        RETURN
END


SUBROUTINE NORM1(NP, OMEGAM, TBATH, Gamm, TEMP, SXXW, OMSTOR, XRE)
        IMPLICIT NONE
        integer::ii, jj, NP
        double precision::omegam

        double precision::PI2, XRE, XIM, COOL, gamm, TEM
        double precision::TEMP2, TEMP, C1, C2, DEL

        double precision::TBATH, hbar, BOLTZ
        PARAMETER(hbar=1.05d-34, BOLTZ=1.4d-23)
        double precision, DIMENSION(NP)::OMSTOR
        double precision::SXXW(NP)


        pi2 = 2.d0 * dacos(-1.d0)
        ! integrate the position spectrum of bead
        ! quick hack - use trapezoidal rule- improve later
        XRE = 0.d0
        XIM = 0.d0
        DEL = abs(OMSTOR(2) - OMSTOR(1))
        DO ii=1, NP-1
                Tem = 0.5d0 * (SXXW(ii) + SXXW(ii + 1))
                XRE = XRE + TEm
        END DO
        XRE = XRE * DEL / pi2
        TEMP= XRE * hbar * OMEGAM / BOLTZ

        100 FORMAT(6E14.6)
        RETURN
END

SUBROUTINE HETERODYNE(NTOT, AVNX, AVNY, AVNZ, AVPHOT, THETA, A1, A1dagg, SHET1)
        IMPLICIT NONE
        integer::ii, m, jj, NTOT
        double precision:: THETA, SHET1, AVNX, AVNY, AVNZ, AVPHOT
        double complex::XI
        !double complex::XTHET1(1, NTOT)
        double complex:: XPM1(1, NTOT), XAM1(1, NTOT)
        double complex:: A1(1, NTOT)
        double complex:: A1dagg(1, NTOT), XHET1(1, NTOT)
        ! zero arrays
        ! XPM1=(0.d0, 0.d0)
        ! XAM1=(0.d0, 0.d0)

        XHET1=(0.d0, 0.d0)

        ! i
        XI=cmplx(0.d0, 1.d0)
        DO ii=1, NTOT
                ! XAM1(1, ii)=A1(1, ii)+A1dagg(1, ii)
                ! XPM1(1, ii)=XI*(A1(1, ii)-A1dagg(1, ii))
                ! XTHET1(1, ii)=XPM1(1, ii)*sin(theta)+XAM1(1, ii)*cos(theta)
                XHET1(1, ii)=A1dagg(1, ii)
        END DO
        SHET1=0.d0

        SHET1 = SHET1 + abs(XHET1(1, 1))**2 * AVPHOT + abs(XHET1(1, 2))**2 * (AVPHOT+1.d0)
        SHET1 = SHET1 + abs(XHET1(1, 3))**2 * AVNX + abs(XHET1(1, 4))**2 * (AVNX+1.d0)
        SHET1 = SHET1 + abs(XHET1(1, 5))**2 * AVNY + abs(XHET1(1, 6))**2 * (AVNY+1.d0)
        SHET1 = SHET1 + abs(XHET1(1, 7))**2 * AVNZ + abs(XHET1(1, 8))**2 * (AVNZ+1.d0)

        100 FORMAT(6D16.8)
        RETURN
END
