SUBROUTINE BEAD_PARAMETERS(thet, Polaris, epsTW, epsCAV, ZR, XM, kappin, kappnano, GAMMAM)
    ! """
    ! subroutine below is provided by user
    ! and calculates  relevant parameters
    ! """
IMPLICIT NONE 
    integer::ii, m, jj
    integer::NPERIOD, NTOT
    double precision::R0, RHO, EPSR, EPSI0, C, hbar, BOLTZ, TEMP, Gravity
    double precision::WK, waist, XL, Finesse, Press, WX, WY, DelFSR
    double precision::PIN1, DET1, X0, Y0, Z0, Theta0

    INCLUDE 'CSCAVITY.h'
    double precision::epsTW, epsCAV, ZR, A, thet
    double precision::W2, W2M, OMOPT, coeff
    double precision::KAPPA, XM, kappnano, kappin
    double precision::pi, pi2
    double precision::Polaris, VOL
    double precision::Gammam

    ! zero eq. initial values
    PI = dacos(-1.d0)
    XM = RHO * 4. * pi/3. * R0**3
    WRITE(6, *) 'Mass of bead = (Kg)'
    WRITE(6, 100) XM
    Polaris=4.*pi * EPSI0 * (EPSR - 1.)/(EPSR + 2.) * R0**3
    WRITE(6, *) 'Polarisability of bead = '
    WRITE(6, 100) Polaris

    OMOPT = C * WK
    ! note waist is waist radius
    W2 = waist**2

    ! add a factor of 4 here. Not in GALAX1-5 routines!!!
    VOL = XL * Pi * W2 / 4.d0

    ! Depth of cavity field. Weak and unimportant for CS case
    A = OMOPT * POLARIS/2. / VOL / EPSI0
    WRITE(6, *) 'Cavity trap A/(2pi) =  (Hz); zero-th shift'
    WRITE(6, 100) A/pi/2.,  A * cos(WK * x0)**2

    KAPPin = pi * c / finesse / XL
    WRITE(6, *) 'kappaIn =  (Hz)'
    WRITE(6, 100) Kappin / 2 / pi

    epsTW = 4.*PIN1 / (Wx * Wy * pi * c * EPSI0)
    epsTW = sqrt(epsTW)
    epsCAV = hbar * OMOPT/(2. * VOL * Epsi0)
    epsCAV = sqrt(epsCAV)
    WRITE(6, *) 'epsTW, epsCAV', epsTW, epsCAV

    coeff = WK * Polaris / EPSI0 / OMOPT**2 / pi
    kappnano = 4.*coeff**2 * DelFSR * cos(WK * x0) * cos(WK * x0)
    WRITE(6, *) 'kappnano', kappnano
    kappa = kappnano + kappin

    ZR = WX * WY * WK / 2.d0
    WRITE(6, *) 'ZR = ', ZR

    ! take usual expression eg Levitated review by Li Geraci etc
    ! 1 bar= 10^ 5 pascal; Press is in mbar = 10^ 2 pascal
    ! gamma=16 P/(pi*v*rho*R)
    ! v=speed of air=500 /s
    GAMMAM = 1600.* press / pi
    GAMMAM = GAMMAM / 500 / RHO / R0
    ! Fix of Feb.2016 our GAMMAM => GAMMAM/2!!
    ! gamma/2=8 P/(pi*v*rho*R)
    GAMMAM = GAMMAM / 2.d0
    WRITE(6, *) 'mechanical damping* 2pi'
    WRITE(6, 100) GAMMAM

    100 FORMAT(4D16.8)

    RETURN
END


SUBROUTINE EQUILIBRIUM_PARAMETERS(&
    thet, Det, Wkx0, Polaris, &
    epsTW, epsCAV, XM, ZR, &
    kappnano, kappin, GammaM, &
    OMX, OMY, OMZ, GMAT, PHON)
    ! """
    ! subroutine below obtains the optomechanical parameters
    ! """
    IMPLICIT NONE
    integer::ii, m, jj, NT, iwrite

    integer::NPERIOD, NTOT
    double precision::R0, RHO, EPSR, EPSI0, C, hbar, BOLTZ, TEMP, Gravity
    double precision::WK, waist, XL, Finesse, Press, WX, WY, DelFSR, Det
    double precision::PIN1, DET1, X0, Y0, Z0, Theta0

    INCLUDE 'CSCAVITY.h'
    double precision::GX, GY, Gz, GXY, GZX, GYZ, XZPF, YZPF, ZZPF, thet
    double precision::OMx, OMy, OMz,  Edip, Ediph, ZR, Nphoton, ALPre, ALPim
    double precision::COOLx, COOLy, COOLz, epsCAV, epsTW, Polaris

    double precision::DET2pi, GammaM
    double precision::KAPPA, kapp2, kappnano, kappin, A, XM
    double precision::pi, pi2
    double precision::WKX0, Wkxx0, Wky0, Wkz0
    double precision::C1, C2, C3, C4
    ! analytic equil phonon numbers
    double precision::PHON(3)
    double precision::GMAT(6)

    iwrite = 6
    PI = dacos(-1.d0)
    ! Det2pi=Det1*2*pi
    Det2pi = Det
    kappa = kappnano + kappin
    kapp2 = 0.5d0 * kappa
    WRITE(6, *) 'kappa / 2 / pi (kHz) =', kappa / 2 / pi

    ! note that detunings include zeroth order correction for linearised versions
    ! as a first pass work out frequencies with Wk0 = 0
    WkX0 = WK * X0

    OmX = Polaris * epsTw**2 / XM / WX**2
    OmY = Polaris * epsTw**2 / XM / Wy**2
    OMz = 0.5d0 * Polaris * epsTw**2 / XM / ZR**2
    ! OmX = (310d3 * 2 * pi)**2
    ! Omy = (275d3 * 2 * pi)**2
    ! Omz = (80d3 * 2 * pi)**2

    WRITE(6, *) 'Wkxeq / pi', Wkx0 / pi
    ! READ(20, *) Wkxx0, Wky0, Wkz0, Alpre, Alpim
    ! Wkx0 = Wkx0 + Wkxx0 * sin(thet) + Wky0 * cos(thet)
    WRITE(6, *) 'Wkxeq / pi corrected', Wkx0 / pi

    ! optomechanical drive frequency
    WRITE(6, *) 'epsTW, epsCAV'
    WRITE(6, 100)  epsTW, epsCAV

    ! Sept 5 we will use negative Edip
    Edip = -0.5d0 * Polaris * epsTW * epsCAV * sin(thet)
    Ediph = Edip / hbar
    WRITE(6, *) 'Edip / 2 / pi / hbar = '
    WRITE(6, 100) Ediph / 2 / pi

    ! photon number in cavity
    ! real part of photon field
    ALPRe = Det2pi * Ediph * cos(Wkx0) / (kapp2**2 + Det2pi**2)
    ALPim = -kapp2 * Ediph * cos(Wkx0) / (kapp2**2 + Det2pi**2)

    Nphoton = Ediph * Ediph * cos(Wkx0) * cos(Wkx0)
    Nphoton = Nphoton / (kapp2**2 + Det2pi**2)
    WRITE(iwrite, *) 'delta, kappa/2, number of photons in cavity'
    WRITE(6, 100) Det2pi / 2 / pi, kapp2/ 2 / pi, Nphoton

    ! ADD CS POTENTIAL CORRECTION to frequency squared
    C1 = -Edip / XM * 2. * ALPRe * WK**2 * cos(Wkx0)
    omx = omx + C1 * sin(thet) * sin(thet)
    omy = omy + C1 * cos(thet) * cos(thet)
    omz = omz - 2. * Edip / XM * ALPre * (WK - 1.d0/ZR)**2 *cos(Wkx0)

    omx = sqrt(omx)
    omy = sqrt(omy)
    omz = sqrt(omz)
    WRITE(iwrite, *) 'mech freq Omx/2pi=', omx/2/pi
    WRITE(iwrite, *) 'mech freq Omy/2pi=', omy/2/pi
    WRITE(iwrite, *) 'mech freq Omz/2pi=', omz/2/pi

    ! optomechanical couplings
    XZPF = sqrt(hbar / (2.d0 * XM * OMX))
    YZPF = sqrt(hbar / (2.d0 * XM * OMY))
    ZZPF = sqrt(hbar / (2.d0 * XM * OMZ))

    GX = Ediph * WK * XZPF * sin(thet) * sin(Wkx0)
    GY = Ediph * WK * YZPF * cos(thet) * sin(Wkx0)
    GZ= -Ediph * (WK - 1.d0 / ZR) * ZZPF * cos(Wkx0)
    WRITE(iwrite, *) 'GX ,  GY,  GZ in Hz', GX/2/pi, GY/2/pi, Gz/2/pi
    ! corrected sign on 29/8/2019
    GXY = Ediph * WK * XZPF * WK * YZPF * ALPre * sin(2. * thet) * cos(Wkx0)
    GZX = 2.d0 * Ediph * (WK - 1.d0 / ZR) * ZZPF * WK * XZPF * ALPim * sin(Wkx0) * sin(thet)
    GYZ=2.d0*Ediph*(WK-1.d0/ZR)*ZZPF*WK*YZPF*ALPim*sin(Wkx0)*cos(thet)
    WRITE (iwrite, *) 'GXY ,  GYZ,  GZX in Hz', GXY/2/pi, GYZ/2/pi, GZX/2/pi
    GMAT(1) = GX
    GMAT(2) = GY
    GMAT(3) = GZ
    GMAT(4) = GXY
    GMAT(5) = GYZ
    GMAT(6) = GZX

    ! calculate cooling rates
    WRITE(iwrite, *) 'COOLING RATES'
    CALL COOLING_RATES(WKX0, DET2pi, Kappa, GX, GY, GZ, OMx, OMy, OMz, COOLx, COOLy, COOLz, Gammam)

    WRITE(iwrite, *) 'X optomechanical cooling rate,  mechanical damping, kx0'
    ! we compare phonon damping with Gamma_tot = 2 * Gammam
    WRITE(iwrite, 100) COOLx, 2. * GAMMAM, Wkx0
    WRITE(iwrite, *) 'Number of phonons: at room temp;at equil'
    C1 = boltz * temp / hbar / omX
    C2 = 2.d0 * C1 * GAMMAM / (abs(COOLx) + 2. * GammaM)
    WRITE(iwrite, 100) C1, abs(C2)
    PHON(1) = abs(C2)
    WRITE(iwrite, *) 'X temperature: at room temp;at equil'
    WRITE(iwrite, *) Temp, Temp * 2. * gammam / (2. * gammam + abs(coolx))

    !  Y
    WRITE(iwrite, *) 'y optomechanical cooling rate,  mechanical damping, Y0'
    ! we compare phonon damping with Gamma_tot = 2 * Gammam
    WRITE(iwrite, 100) COOLy, 2.*GAMMAM, y0
    WRITE(iwrite, *) 'Number of phonons: at room temp;at equil'
    C1 = boltz * temp / hbar / omy
    C2 = 2.d0 * C1 * GAMMAM / (abs(COOLy) + 2. * gammam)
    WRITE(iwrite, 100)C1, abs(C2)
    WRITE(iwrite, *) 'Y temperature: at room temp;at equil'
    WRITE(iwrite, *) Temp, Temp * 2. * gammam / (2. * gammam + abs(Cooly))
    PHON(2) = abs(C2)

    !  Z
    WRITE(iwrite, *) 'z optomechanical cooling rate,  mechanical damping,  Z0'
    ! we compare phonon damping with Gamma_tot = 2 * Gammam (gammam is act gamm/2)
    WRITE(iwrite, 100) COOLz, 2. * GAMMAM, z0
    WRITE(iwrite, *) 'Number of phonons: at room temp;at equil'
    C1 = boltz * temp / hbar / omz
    C2 = 2.d0 * C1 * GAMMAM / (abs(COOLz) + 2. * gammam)
    WRITE(iwrite, 100) C1, abs(C2)
    WRITE(iwrite, *) 'Z temperature: at room temp;at equil'
    WRITE(iwrite, *) Temp, Temp * 2. * gammam / (2. * gammam + abs(Coolz))
    PHON(3) = abs(C2)

    100 FORMAT(4D16.8)

    RETURN
END


SUBROUTINE COOLING_RATES(WKX0, DETUN, Kappa, GX, GY, GZ, OMx, OMy, OMz, COOLx, COOLy, COOLz, Gammam)
    ! """
    ! subroutine below obtains the optomechanical parameters
    ! Note that DETUN effectively is
    ! DEtun+AOPT. It is corrected by zeroth order optical contribution
    ! in the linearised case
    ! """

    IMPLICIT NONE
    integer::ii, m, jj, iwrite
    integer::NPERIOD, NTOT
    double precision::R0, RHO, EPSR, EPSI0, C, hbar, BOLTZ, TEMP, Gravity
    double precision::WK, waist, XL, Finesse, Press, WX, WY, DelFSR
    double precision::PIN1, DET1, X0, Y0, Z0, Theta0

    INCLUDE 'CSCAVITY.h'
    double precision::DETUN, GammaM, GX, GY, Gz, OMx, OMy, OMz
    double precision::KAPPA, KAPP2
    double precision::pi, pi2
    double precision::WKX0
    double precision::C1, C2, C3, C4, COOLx, COOLy, COOLz

    PI=dacos(-1.d0)
    Kapp2=0.5d0*Kappa
    WRITE(iwrite, *) 'GX ,  GY,  GZ in Hz', GX/2/pi, GY/2/pi, Gz/2/pi

    ! X cooling
    ! for maximum cooling sin2kx_0=1
    ! C1 = sin(2.d0 * WK0)**2 * 2.* Kapp2 * WK**2 * A**2 * ASQ * hbar / XM / 2.d0 / OMX
    C1 = GX**2 * kappa
    C2 = DETUN
    ! WRITE(iwrite, *) 'shifted detuning = det + A'
    ! WRITE(iwrite, 100) C2 / 2.d0 / pi

    C3=(C2+omX)**2+kapp2**2
    C3=1.d0/C3
    C4=(C2-omx)**2+kapp2**2
    C4=1.d0/C4
    COOLx=-C1*(C3-C4)

    ! Y cooling
    C1 = Gy**2 * kappa
    C2 = DETUN
    ! WRITE(iwrite, *)'shifted detuning=det+A'
    ! WRITE(iwrite, 100)C2/2.d0/pi

    C3 = (C2 + omY)**2 + kapp2**2
    C3 = 1.d0 / C3
    C4 = (C2 - omY)**2 + kapp2**2
    C4 = 1.d0 / C4
    COOLy = -C1 * (C3 - C4)

    ! Z cooling
    C1 = Gz**2 * kappa
    C2 = DETUN
    ! WRITE(iwrite, *) 'shifted detuning = det + A'
    ! WRITE(iwrite, 100) C2/ 2.d0 / pi
    ! here we neglect opto shift!
    C3 = (C2 + omz)**2 + kapp2**2
    C3 = 1.d0/C3
    C4 = (C2-omz)**2 + kapp2**2
    C4 = 1.d0/C4
    COOLz = -C1 * (C3 - C4)

    100 FORMAT(4D16.8)
    RETURN
END


