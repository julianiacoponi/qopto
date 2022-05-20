! input parameters
! NPERIOD=number of periods in prop. (of fastest oscillation: kappa^-1 or om_M)
! NTOT=number of equations so 8=>1 optical mode +3D
! R0= sphere radius
! XL=cavity length 
! Pin1=input power in Watts tweezer beam
!
! Press = air pressure in millibars
! XL=cavity length
! rho=sphere density
! WK=2*pi/lambda=k
! DelFSR = 1 FSR= 14 GHz
! DET1=detuning in KHz trap beam
! Theta0= angle between tweezer polarization and cavity axis. Given as as
! FRACTION of pi so pi/4  is theta0=0.25
! waist is  waist radius
!
! RR0 and V0 = ion trap parameters
 PARAMETER(NPERIOD=80000,NTOT=4,R0=71.5d-9, &
 !   PARAMETER(NPERIOD=160000,NTOT=4,R0=71.5d-9, &
            rho=2198.d0,EPSR=2.1d0,Epsi0=8.854d-12, &
!              rho=2198.d0,EPSR=1.45d0**2,Epsi0=8.854d-12, &
    c=3.d8,hbar=1.05d-34,BOLTZ=1.4d-23,TEMP=300.d0,Gravity=9.8d0,&
    WK=5.9d6,waist=61.d-6,WX=1.02d-6,Wy=0.879d-6,XL=1.223d-2,Finesse=2.8d4, &
    Press=1.e-6,Pin1=0.3886d0,DET1=-200.d3,DelFSR=14.d9,theta0=0.5d0,&
! Equilibrium positions 
	    !  X0=0.125*lambda ,Y0=waist/sqrt(2)
	    Y0=0.0d-6,X0=0.25d0*1.064d-6,Z0=0.d-6)
