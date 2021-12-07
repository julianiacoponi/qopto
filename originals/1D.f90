!This code has 1 beam plus 3 dimensions
!**************************************
! Calculates quantum noise spectra.
! USES real mean fields
! works out position and momentum optical quadrature
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! GammaM=mechanical damping due to gas, also noise from fluctuation dissipation.
!
!*********************************************
! NB no intent(in) or intent(out) used- add these later!!
!*********************************************************
!*********************************************
     IMPLICIT NONE
     integer  ::ll,kk,mm,it,im,ii,jj,NT,ithet,Nthet,Nx0,ix0,Ndet

        integer  ::NTOT,NPERIOD,NPTS
     double precision::R0,RHO,EPSR,EPSI0,C,hbar,BOLTZ,TEMP,Gravity
     double precision::WK,waist,XL,Finesse,Press,WX,WY,DelFSR
      double precision::PIN1,DET1,X0,Y0,Z0,Theta0,Thetahom,Delx0,XX0,Deldet,DET2PIini,test
      double precision::Delthet,cavphotn
! parameter file with input values
     include 'CSCAVITY-1D.h'
       Parameter(Nthet=1,Delthet=0.01d0,Ndet=1,Deldet=20.d3,Thetahom=0.d0,NPTS=10000)
        double precision::OMx,OMy,OMz
!calculated  parameters
!
double precision::Det2pi,kappnano,kappin,kapp2,kappa,Polaris,thet,THEhom
double precision::XM,GAMMAM,Wkx0,epsTW,epsCAV,ZR,Ed
       double precision:: XMAX,DEL,GX

        double precision::pi,pi2,OMsweep,XXQM,SHOM1,AHOM1,SHET1,XXQMs
        double precision:: TPERIOD,AVNX,AVPHOT,AVRE,OPTOCOOL,PHONONS
       double precision::COOLX,TEMPX,TBATH
! read parameter file with input values
       DOUBLE PRECISION::GMAT(6)
      DOUBLE PRECISION::PHON(3),AV(3)
        DOUBLE PRECISION, DIMENSION(NPTS):: SXXQM,OMSTOR,Shom
         DOUBLE COMPLEX::XI,ZOM

!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!***************************************************************
! Here plot optical field amplitudes
     OPEN(8,file="FTanX.dat",status="unknown")
     OPEN(10,file="FTanOPT.dat",status="unknown")
    open(14,file="PHONS.dat",status="unknown")
!
   pi=dacos(-1.d0)
   pi2=2.d0*pi
   DET2PIini=DET1*pi2
    TBATH=TEMP
    NT=NTOT

!  OPEN LOOP OVER DEtuning or homodyne angle.   &&&&&&&&&&&&&&&&&&&&&&&
!     do  10 ithet=1,Nthet
    do  10 ix0=1,Ndet
     DET2PI=DET2PIini+Deldet*pi*2*(ix0-1)
      write(6,*)'detuning',Det2pi/2/pi
!     thet=theta0+Delthet*(ithet-1)
      thet=theta0
!     XX0=Wk*x0+(ix0-1)*Delx0
       xx0=Wk*x0
      Wkx0=xx0
     thet=thet*pi
      THEhom=Thetahom*pi
      write(6,*)'theta/pi', thet/pi
! note BEAD routine uses x0 but not theta

       
      call BEAD_parameters(Wkx0,Polaris,epsTW,epsCAV,ZR,XM,kappin,kappnano,GAMMAM)
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  call EQUIL_PAR(thet,Det2pi,Wkx0,Polaris,epsTW,epsCAV,XM,ZR,kappnano,kappin,GammaM,OMX,OMY,OMZ,GMAT,PHON)
    write(12,120)thet,(abs(Gmat(ii)),ii=1,6)

!-----------------------------------------------------
     GX=GMAT(1)
120 Format(7E14.6)
     Ed=0.5d0*Polaris*epsTW*epsCAV*sin(thet)

      kappa=(kappnano+kappin)
      kapp2=kappa*0.5d0
      write(6,*)'kappa in main',kappa/2/pi,kapp2/2/pi
     write(6,*)'GammaM=',GammaM
     write(6,*)'OMX in Hz=',OMX/pi2! thermal bath occupancy
      AVNX=BOLTZ*TBATH/hbar/OMX
     
! shot noise
      AVPHOT=0.d0
!***************************************************************
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! open loop over frequency for noise spectra eg S_xx(\omega)=FT(autocorrelation <X(t)X^T(t+tau)>)
            XMAX=2*OMX*1.001
            DEL=XMAX/NPTS
             SXXQM=0.d0
   
             SHOM=0.d0
             do 11 mm=1,NPTS
             OMsweep=-XMAX+2*(mm-1)*DEL
! store frequency for integration
             OMSTOR(mm)=OMsweep
! work out PSD of homodyne
        AHOM1=0.d0
! symmetrised 
       XXQMs =0.d0
      
! ! work out  for negative frequency for symmetrisation
     
    CALL ANALYT(GX,AVNX,AVPHOT,THEhom,DET2pi,Kapp2,GAMMAM,OMX,-OMsweep,XXQM,SHOM1,SHET1)

! update homodyne and symm disp.
      AHOM1=AHOM1+SHOM1

! work out same for positive and symmetrise homodyne
    CALL ANALYT(GX,AVNX,AVPHOT,THEhom,DET2pi,Kapp2,GAMMAM,OMX,OMsweep,XXQM,SHOM1,SHET1)
! update
       AHOM1=0.5*(AHOM1+SHOM1)
      XXQMs=XXQMs+XXQM
    


   write(8,101)OMsweep/2/pi,XXQMs
 write(10,101)OMsweep/2/pi,SHET1,AHom1
101   format(7D14.6)
             SXXQM(mm)=XXQM
             
          
! to find optimal squeezing
              SHOM(mm)=AHOM1
11            enddo

AVRE=(omx+Det2pi)**2+kapp2**2
AVRE=AVRE/4/omx/(-det2pi)
    
write(6,*)'X: back action limited phonons',AVRE
  
  
! integrate and normalise the quantum  noise spectra

      CALL NORM1(NPTS,OMX,TBATH,GAMMAM,TEMPX,SXXQM,OMSTOR,AVRE)
! Area AVRE corresponds to 2n+1 so convert to get n
                PHONONS=PHON(1)
     AVRE=0.5d0*(AVRE-1.d0)
     AV(1)=AVRE
      write(6,*)'X phonons from formula,  from SXX FT'
      write(6,200)PHONONS,AVRE

     write(14,120)Det2pi/2/pi,AV(1),PHON(1)
! loop over theta
10    enddo

!***********************************************************************
100   FORMAT(I3,3E14.6,1x,2(E14.6))
200  FORMAT(7E14.6)

    STOP
    END
!********************************************************************
!  Function to evaluate optomechanical cooling formula
!********************************************************************

     FUNCTION optocool(G1,Det1X,KAPP2,OMEGAM,GAMMAM)
       Implicit None
      double precision::G1,Det1X,KAPP2,OMEGAM

       double precision::C1,C2,C3,C4
       double precision::OPTOCOOL,COOL1,COOL2,GAMMAM
       double precision:: hbar,BOLTZ,TBATH
       PARAMETER(BOLTZ=1.4d-23,hbar=1.05d-34,TBATH=300.0)
        COOL1=0.d0


! now work out opto cooling expression 
! Trap beam
       C1=2.*KAPP2*G1*G1
       C2=DET1X


      C3=(C2+omegam)**2+kapp2**2
       C3=1.d0/C3
       C4=(C2-omegam)**2+kapp2**2
       C4=1.d0/C4
       COOL1=-C1*(C3-C4)

        optocool=COOL1
100   format(4D16.6)
   return
    end
!********************************************************************
!  Generic routine for noise spectra of trap and probe beams
!********************************************************************

     SUBROUTINE ANALYT(GX,AVNX,AVPHOT,THETA,DET,Kapp2,GAMMAM,OMX,OMEGA,XXF,SHOM1,SHET1)

!************************************************************************
!************************************************************************
                
	IMPLICIT NONE
        integer  ::ii,m,jj,NT
         PARAMETER(NT=4)
        double precision:: DET
        double precision::gamm
        double precision::pi,pi2,XXF,PHN,SHOM1,SHET1
        double precision:: OMEGA,KAPP2,GAMMAM,OMX,OMY,OMZ,omhet,oma
        double precision:: CHIMSQ,XNORMSQ
       DOUBLE PRECISION:: AVNX,AVPHOT,THETA,Gav
       DOUBLE PRECISION::GX
        double complex::CHIR1,CHISTMOM1
       double complex::CHIMX,CHIMMOMX
      double complex:: BAX(1,NT),A1(1,NT),A1dagg(1,NT)


       pi=dacos(-1.d0)

       SHOM1=0.d0
        BAX=0.d0

! *******WORK OUT NOISE SPECTRA
! First do susceptibilities
                    
        call SUSCEPT(OMEGA,DET,Kapp2,gammam,OMX,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX)


! work out noise vector for x, a1 and a2
       Call Avect(GX,Gammam,kapp2,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX,A1,A1dagg,BAX)
!     XX= sqrt(0.5) (b+b^dagg) so halve XX
      XXF=(ABS(BAX(1,1)))**2
     XXF=XXF+ (AVNX+1)*(ABS(BAX(1,3)))**2+AVNX*(ABS(BAX(1,4)))**2

! work out homodyne spectra using same vectors
      call Homodyne(NT,AVNX,AVPHOT,THETA,A1,A1dagg,SHOM1)


     SHOM1=(SHOM1-1.d0)/abs(CHIR1-CHISTMOM1)**2/GX**2/kapp2/2.
! work out heterodyne for positive frequency branch.
!  Shift frequency by heterodyne beat.

!--------------------------------------------------------------------------
! FOR THIS EXAMPLE JUST QUADRUPLE MECHANICAL FREQUENCY.
       OMA=OMx*4.d0
       OMA=0.d0
       OMHET=(OMEGA+OMa)
!-----------------------------------
! work out susceptibilities shifted in frequency

        call SUSCEPT(OMhet,DET,Kapp2,gammam,OMX,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX)

! work out noise vector again
      Call Avect(GX,Gammam,kapp2,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX,A1,A1dagg,BAX)

    call Heterodyne(NT,AVNX,AVPHOT,THETA,A1,A1dagg,SHET1)
!     write(6,*)omhet,SHET1,(SHET1)/ABS(CHIR1)**2/GMAT(1)**2
     SHET1=(SHET1-1.d0)/ABS(CHISTMOM1)**2/Gav**2/kapp2/2.
 !     SHET1=SHET1-1.d0
100   format(6D16.8)
      RETURN
       END
!*************************************************
 ! Work out homodyne spectra
!**************************************************
      Subroutine Homodyne(NTOT,AVNX,AVPHOT,THETA,A1,A1dagg,SHOM1)
        IMPLICIT NONE
         integer  ::ii,m,jj,NTOT
         double precision:: THETA,SHOM1,AVNX,AVPHOT
  
         double complex::XI
        double complex::XTHET1(1,NTOT)
         double complex:: XPM1(1,NTOT),XAM1(1,NTOT)
      double complex:: A1(1,NTOT)
      double complex:: A1dagg(1,NTOT)
!  zero arrays
           XPM1=(0.d0,0.d0)
           XAM1=(0.d0,0.d0)

           XTHET1=(0.d0,0.d0)

! i
        XI=CMPLX(0.d0,1.d0)
       do ii=1,NTOT
        XAM1(1,ii)=A1(1,ii)+A1dagg(1,ii)
        XPM1(1,ii)=XI*(A1(1,ii)-A1dagg(1,ii))

        XTHET1(1,ii)=XPM1(1,ii)*sin(theta)+XAM1(1,ii)*cos(theta)
           enddo
        SHOM1=0.d0

       SHOM1=SHOM1+ABS(XTHET1(1,1))**2*AVPHOT+ABS(XTHET1(1,2))**2*(AVPHOT+1.d0)
       SHOM1=SHOM1+ABS(XTHET1(1,3))**2*AVNX+ABS(XTHET1(1,4))**2*(AVNX+1.d0)
 
100   format(6D16.8)
RETURN
END


!*************************************************************************
! work out susceptibilities for noise spectra
!*************************************************************************
 SUBROUTINE SUSCEPT(OMEGA,DET1x,Kapp2,gamm,OMX,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX)
	IMPLICIT NONE 
  
        double precision:: DET1X
        double precision::gamm,gamm2,Kapp2
        double precision:: OMEGA,OMX
        double precision:: t1,t2,ANORM,BNORM
        double complex::CHIR1,CHISTMOM1,CHIMX,CHIMMOMX
! gamm is actually gamma/2
        Gamm2=gamm
! FIELD1 susceptibilities

!*****************************************
!Chi_R
       ANORM=KAPP2**2+ (omega+DET1X)**2
       t1=kapp2/ANORM
       t2= (omega+DET1X)/ANORM
       CHIR1=CMPLX(t1,t2)
! chi_r^*(-omega)
      ANORM=KAPP2**2+ (-omega+DET1X)**2
      t1=kapp2/ANORM
      t2= (-omega+DET1X)/ANORM
      CHISTMOM1=CMPLX(t1,-t2)
!******************************************
! X MECHANICAL susceptibilities
! chi_M X
        BNORM=(Gamm2)**2 + (omega-OMX)**2
      t1= (Gamm2)/BNORM
      t2=(omega-OMX)/BNORM
      CHIMX=CMPLX(t1,t2)
! CHI_M*(-om) X
      T1=(Gamm2)**2 + (-omega-OMX)**2
     CHIMMOMX=CMPLX(Gamm2/T1, (omega+OMX)/T1)
!****************************************************

!****************************************************
100   format(6D16.8)
      RETURN
       END



!*****************************************************************
! Work our vectors of X and optical fields but in terms of noise operators
!*******************************************************************
  Subroutine Avect(GX,Gamm,kapp2,CHIR1,CHISTMOM1,CHIMX,CHIMMOMX,A1,A1dagg,BAX)

       IMPLICIT NONE
        integer  ::ii,m,jj,NTOT,NPERIOD
        double precision:: ANORM,BNORM
        double precision:: t1,t2,Sqrtkapp,Sqrtgamm,Gamm,Gamm2,kapp2
        double precision::GX            
       double precision::ASQ1,alp1
        PARAMETER(NTOT=4)
        double complex::CHIR1,CHISTMOM1
       double complex::CHIMX,CHIMMOMX
        double complex:: XI,ONE,Betx
  
       double complex:: CMX,CSUM,CNORM,eta0c,etaMpi2c,etaPpi2c,etaX,CA1,CA1dagg
      double complex,intent(out):: BAX(1,NTOT)
      double complex, intent(out):: A1(1,NTOT)
      double complex, intent(out):: A1dagg(1,NTOT)
      double complex:: N0X(NTOT)    

  
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
! actually Gamm is gamma/2
       Gamm2=Gamm
! zero arrays
        BAX=CMPLX(0.d0,0.d0)
        N0X=CMPLX(0.d0,0.d0) 
        A1=CMPLX(0.d0,0.d0)
        A1dagg=CMPLX(0.d0,0.d0)
! i
        XI=CMPLX(0.d0,1.d0)
        ONE=CMPLX(1.d0,0.0d0)

!      SUSCEPTIBILITIES
        eta0c=CHIR1-CHISTMOM1
        etaMpi2c=XI*(CHIR1+CHISTMOM1)
        etaPpi2c=-XI*(CHIR1+CHISTMOM1)
        etaX=CHIMX-CHIMMOMX


! NORMALIZATIONS
        CMX=1.d0+GX**2*etax*eta0c
         Sqrtkapp=sqrt(2.d0*KAPP2)
         Sqrtgamm=sqrt(2.d0*Gamm2)
        BETx=XI*Sqrtkapp*etax*GX
  
!&&&&&&&&&&&&&&&&&&&&&&&&
!
! 1D X noise vector; weights of a1,a1*,bx,bx*
         N0X(1)=BETX*CHIR1/CMX
         N0X(2)=BETX*CHISTMOM1/CMX
         N0X(3)=Sqrtgamm*CHIMX/CMX
         N0X(4)=Sqrtgamm*CHIMMOMX/CMX
! FOR 1D
      do ii=1,NTOT
      BAX(1,ii)=N0X(ii) 
      enddo
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! now work out the optical trap field =a1
         CA1=XI*CHIR1
! now work out the photon field a1dagger
         CA1dagg=-XI*CHISTMOM1

         do ii=1,NTOT
  
   A1(1,ii)=CA1*GX*BAX(1,ii)
   A1dagg(1,ii)=CA1dagg*GX*BAX(1,ii)
        enddo
! add shot or incoming noise
! trap beam: add cavity-filtered contribution
        A1(1,1)=A1(1,1)+Sqrtkapp*CHIR1
       A1dagg(1,2)=A1dagg(1,2)+Sqrtkapp*CHISTMOM1

! cavity output : add incoming imprecision
! work out a_out=a_in-Sqrtkapp(a)
      do ii=1,NTOT
      A1(1,ii)=-A1(1,ii)*Sqrtkapp
      A1dagg(1,ii)=-A1dagg(1,ii)*Sqrtkapp
      enddo
      A1(1,1)=ONE+A1(1,1)
     A1dagg(1,2)=ONE+A1dagg(1,2)


         return
         end


!*************************************************************************
!************************************************************************

       SUBROUTINE NORM1(NP,OMEGAM,TBATH,Gamm,TEMP,SXXW,OMSTOR,XRE)

!************************************************************************
!************************************************************************
                
	IMPLICIT NONE 
        integer  ::ii,jj,NP
         double precision::omegam

       double precision::PI2,XRE,XIM,COOL,gamm,TEM
      double precision::TEMP2,TEMP,C1,C2,DEL
 
        double precision::TBATH,hbar,BOLTZ
       PARAMETER(hbar=1.05d-34,BOLTZ=1.4d-23)
         DOUBLE PRECISION, DIMENSION(NP):: OMSTOR
        DOUBLE PRECISION:: SXXW(NP)
     

         pi2=2.d0*dacos(-1.d0)
!  integrate the position spectrum of bead
! quick hack - use trapezoidal rule- improve later
        XRE=0.d0
        XIM=0.d0
        DEL=ABS(OMSTOR(2)-OMSTOR(1))
        do ii=1,NP-1
        Tem=0.5d0*(SXXW(ii)+SXXW(ii+1))
         XRE=XRE+TEm
         enddo
         XRE=XRE*DEL/pi2
        TEMP=XRE*hbar*OMEGAM/BOLTZ

 100       format(6E14.6) 
        return
         end 
!************************
!**************************************************

Subroutine Heterodyne(NTOT,AVNX,AVPHOT,THETA,A1,A1dagg,SHET1)
IMPLICIT NONE
integer  ::ii,m,jj,NTOT
double precision:: THETA,SHET1,AVNX,AVPHOT

double complex::XI
!double complex::XTHET1(1,NTOT)
double complex:: XPM1(1,NTOT),XAM1(1,NTOT)
double complex:: A1(1,NTOT)
double complex:: A1dagg(1,NTOT),XHET1(1,NTOT)
!  zero arrays
!         XPM1=(0.d0,0.d0)
!         XAM1=(0.d0,0.d0)

XHET1=(0.d0,0.d0)

! i
XI=CMPLX(0.d0,1.d0)
do ii=1,NTOT
!        XAM1(1,ii)=A1(1,ii)+A1dagg(1,ii)
!        XPM1(1,ii)=XI*(A1(1,ii)-A1dagg(1,ii))
!        XTHET1(1,ii)=XPM1(1,ii)*sin(theta)+XAM1(1,ii)*cos(theta)
XHET1(1,ii)=A1dagg(1,ii)
enddo
SHET1=0.d0

SHET1=SHET1+ABS(XHET1(1,1))**2*AVPHOT+ABS(XHET1(1,2))**2*(AVPHOT+1.d0)
SHET1=SHET1+ABS(XHET1(1,3))**2*AVNX+ABS(XHET1(1,4))**2*(AVNX+1.d0)



100   format(6D16.8)
RETURN
END
