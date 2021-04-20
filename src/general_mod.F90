!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!      Updatewfn, Get_KinCoul, Get_FCKinCoul, Get_Nuclearpotential,
!         ORTHONORMALIZE, ADJUSTSIGN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE general_mod

  USE io_tools
  USE atomdata
  USE numerov_mod
  USE globalmath
  USE gridmod
  USE radialDirac
  USE radialked
  USE radialsr

  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  SUBROUTINE Updatewfn(Grid,Pot,Orbit,rvin,success)
  !      Given new potential rvin, generate new Orbit%wfn,Orbit%eig,Pot%den
  !
  !       note that rvin=rvh+rvx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Updatewfn(Grid,Pot,Orbit,rvin,success)
    TYPE (GridInfo), INTENT(INOUT) :: Grid
    TYPE (PotentialInfo), INTENT(INOUT) :: Pot
    TYPE (OrbitInfo), INTENT(INOUT) :: Orbit
    REAL(8), INTENT(IN) ::  rvin(:)
    LOGICAL :: success

    !  program to calculate wavefunctions given potential rvin

    INTEGER :: icount,i,j,k,n,it,start,np,ierr,nroot,s1,s2
    INTEGER :: is,ip,id,jf,ig,io,l,nfix,ir,nzeff,jierr,nz,kappa
    REAL(8) :: h,emin,zz
    REAL(8), ALLOCATABLE :: dum(:)
    LOGICAL :: OK

    n=Grid%n; h=Grid%h;    nz=Pot%nz;   zz=Pot%zz
    success=.TRUE.

    allocate(dum(n))

    !!! testing
    write (STD_OUT,*) ' in Updatewfn   needvtau ', needvtau

    Pot%rv=rvin+Pot%rvn(1:n)
    dum=rvin
    if (needvtau) dum(1)=dum(1)-Pot%rvx(1)
    CALL zeropot(Grid,dum,Pot%v0,Pot%v0p)
    IF (ABS(Pot%v0)> 1.d6) Pot%v0=0
    IF (ABS(Pot%v0p)> 1.d6) Pot%v0p=0
    
    IF (finitenucleus) then
            Pot%v0=Pot%v0+Pot%Nv0
            Pot%v0p=Pot%v0p+Pot%Nv0p
    Endif        

    !  solve for bound states of Schroedinger equation
    !
    icount=0
    jierr=0
    it=0
    !  s states :
    IF (Orbit%nps.GT.0) THEN
       it=it+1
       emin=-nz*nz-0.1d0
       l=0
       nroot=Orbit%nps
       start=1;s1=start;s2=start+nroot-1
       IF (scalarrelativistic) THEN
          Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE IF (diracrelativistic) THEN
          kappa=-1     
          Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
       ELSE IF (needvtau) THEN
         write(std_out,*) 'about to call boundked ', nz,emin      
         Call boundked(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE
          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       ENDIF
          IF (.NOT.OK) THEN
             success=.FALSE.
          ENDIF
    ENDIF
    !  p states :
    IF (Orbit%npp.GT.1) THEN
       it=it+1
       emin=-nz*nz/4.d0-0.5d0
       l=1
       nroot=Orbit%npp-1
       s1=s2+1;s2=s1+nroot-1
       IF (scalarrelativistic) THEN
          Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE IF (diracrelativistic) THEN
          kappa=1     
          Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
          s1=s2+1;s2=s1+nroot-1
          kappa=-2     
          emin=-nz*nz/4.d0-0.5d0
          Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
       ELSE IF (needvtau) THEN
         Call boundked(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE
          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       ENDIF
          IF (.NOT.OK) THEN
             success=.FALSE.
          ENDIF
    ENDIF
    !  d states :
    IF (Orbit%npd.GT.2) THEN
       it=it+1
       emin=-nz*nz/9.d0-0.5d0
       l=2
       nroot=Orbit%npd-2
       s1=s2+1;s2=s1+nroot-1
       IF (scalarrelativistic) THEN
          Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE IF (diracrelativistic) THEN
          kappa=2     
          Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
          kappa=-3     
       s1=s2+1;s2=s1+nroot-1
       emin=-nz*nz/9.d0-0.5d0
          Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
       ELSE IF (needvtau) THEN
         Call boundked(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE
          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       ENDIF
          IF (.NOT.OK) THEN
             success=.FALSE.
          ENDIF
    ENDIF
    !  f states :
    IF (Orbit%npf.GT.3) THEN
       it=it+1
       emin=-nz*nz/16.d0-0.5d0
       l=3
       nroot=Orbit%npf-3
       s1=s2+1;s2=s1+nroot-1
       IF (scalarrelativistic) THEN
          Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE IF (diracrelativistic) THEN
          kappa=3     
          Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
          kappa=-4     
       s1=s2+1;s2=s1+nroot-1
       emin=-nz*nz/16.d0-0.5d0
          Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
       ELSE IF (needvtau) THEN
         Call boundked(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE
          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       ENDIF
          IF (.NOT.OK) THEN
             success=.FALSE.
          ENDIF
    ENDIF
    !  g states :
    IF (Orbit%npg.GT.4) THEN
       it=it+1
       emin=-nz*nz/25.d0-0.5d0
       l=4
       nroot=Orbit%npg-4
       s1=s2+1;s2=s1+nroot-1
       IF (scalarrelativistic) THEN
          Call Boundsr(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE IF (diracrelativistic) THEN
          kappa=4     
          Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
          kappa=-5     
       s1=s2+1;s2=s1+nroot-1
       emin=-nz*nz/25.d0-0.5d0
          Call BoundD(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             Orbit%lwfn(:,s1:s2),kappa,nroot,emin,ierr,OK)
       ELSE IF (needvtau) THEN
         Call boundked(Grid,Pot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),&
&             l,nroot,emin,ierr,OK)
       ELSE
          CALL BoundNumerov(Grid,Pot%rv,Pot%v0,Pot%v0p,Pot%nz,&
&              l,nroot,Orbit%eig(s1:s2),Orbit%wfn(:,s1:s2),BDsolve,OK)
       ENDIF
       IF (.NOT.OK) THEN
          success=.FALSE.
       ENDIF
    ENDIF

    deallocate(dum)

  END SUBROUTINE Updatewfn

  SUBROUTINE Get_KinCoul(Grid,Pot,Orbit,SCF,noalt)
    !  program to calculate Kinetic energy and Coulomb Energies from Orbit%wfn
    !   also update Pot%rvh
    TYPE(GridInfo), INTENT(INOUT) :: Grid
    TYPE(PotentialInfo), INTENT(INOUT) :: Pot
    TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
    TYPE(SCFInfo), INTENT(INOUT) :: SCF
    LOGICAL, OPTIONAL :: noalt

    REAL(8) :: ecoul,eex,etot,ekin,eone,h,x,v0,qcal,etxc,small,rescale,fpi
    REAL(8) :: electrons,fac,xocc
    INTEGER :: icount,i,j,k,l,m,n,io
    REAL(8), ALLOCATABLE :: dum(:)
    REAL(8) :: small0=1.d-6
    INTEGER :: counter=1


    n=Grid%n; h=Grid%h
    small=small0

    ALLOCATE(dum(n),STAT=k)
    IF (k /= 0) THEN
       WRITE(STD_OUT,*) 'Error in Get_KinCoul allocation ', n,k
       STOP
    ENDIF

    !update density
!   Note that kinetic energy density (tau) is in Rydberg units
    Orbit%den=0.d0;Orbit%tau=0.d0

    DO io=1,Orbit%norbit
       IF (Orbit%occ(io).GT.small) THEN
         DO i=1,Grid%n      
            IF (ABS(Orbit%wfn(i,io))<machine_zero)Orbit%wfn(i,io)=0
            IF (diracrelativistic) then
              IF (ABS(Orbit%lwfn(i,io))<machine_zero)Orbit%lwfn(i,io)=0
            ENDIF
         ENDDO   
         CALL taufromwfn(Grid,Orbit%wfn(:,io),Orbit%l(io),Orbit%otau(:,io))
         xocc=Orbit%occ(io)
         DO i=1,Grid%n
            Orbit%tau(i)=Orbit%tau(i)+xocc*Orbit%otau(i,io)
            Orbit%den(i)=Orbit%den(i)+xocc*(Orbit%wfn(i,io)**2)
            IF (diracrelativistic) then
              Orbit%den(i)=Orbit%den(i)+xocc*((Orbit%lwfn(i,io))**2)
            ENDIF
         ENDDO
       ENDIF
    ENDDO
    !Kinetic energy density is in Ry (no need of 1/2 factor)
    !Orbit%tau=0.5d0*Orbit%tau

    qcal=integrator(Grid,Orbit%den)
    WRITE(STD_OUT,*) 'qcal = ', qcal

    !For FC, Big change here , electrons are total , not valence anymore!
    electrons=Pot%q
    rescale=electrons/qcal
    Orbit%den(1:n)=Orbit%den(1:n)*rescale
    Orbit%tau(1:n)=Orbit%tau(1:n)*rescale
    WRITE(STD_OUT,*) 'rescaled qcal = ', integrator(Grid,Orbit%den), Pot%q

!   Determine difference with tauW (Weizsaker)       
    fpi=4*pi
    dum(2:Grid%n)=Orbit%den(2:Grid%n)/(fpi*Grid%r(2:Grid%n)**2)
    CALL extrapolate(Grid,dum)
    CALL derivative(Grid,dum,Orbit%deltatau)
    Do i=1,Grid%n
      if (dum(i)>machine_zero) then
        Orbit%deltatau(i)=0.25d0*(Orbit%deltatau(i)**2)/dum(i)
      else
        Orbit%deltatau(i)=0.d0
      endif
    enddo        
    dum(2:Grid%n)=Orbit%tau(2:Grid%n)/(fpi*Grid%r(2:Grid%n)**2)
    call extrapolate(Grid,dum)
    Orbit%deltatau=dum-Orbit%deltatau

    CALL poisson(Grid,Pot%q,Orbit%den,Pot%rvh,ecoul,v00=v0)

    dum=0
    dum(2:n)=Pot%rvn(2:n)*Orbit%den(2:n)/Grid%r(2:n)
    SCF%estatic=integrator(Grid,dum)+ecoul

    !WRITE(STD_OUT,*) ' n  l     occupancy       energy'
    ekin=0.0d0 ; if (frozencorecalculation) ekin=SCF%corekin
    eone=0.0d0
    DO io=1,Orbit%norbit
       if(.not.frozencorecalculation &
&          .or.frozencorecalculation.and.(.not.Orbit%iscore(io))) then
         !WRITE(STD_OUT,'(i2,1x,i2,4x,1p,2e15.7)') &
         !&  Orbit%np(io),Orbit%l(io),&
         !&  Orbit%occ(io),Orbit%eig(io)
         eone=eone+Orbit%occ(io)*Orbit%eig(io)
         IF (counter>1.and..not.present(noalt).and..not.needvtau) THEN
            CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io),Pot%rv,x)
         ELSE
            !CALL kinetic(Grid,Orbit%wfn(:,io),Orbit%l(io),x)
            x=integrator(Grid,Orbit%otau(:,io))
         ENDIF
         ekin=ekin+Orbit%occ(io)*x
       endif
    ENDDO

    write(std_out,*) 'KinCoul check ekin ', ekin,integrator(Grid,Orbit%tau)
!    if (needvtau) then
!       do i=2,Grid%n
!           write(500+counter,'(1p,3e20.9)') Grid%r(i),Orbit%tau(i) &
!&                        , Orbit%den(i)/(Grid%r(i)**2)                   
!       enddo
!    endif
    SCF%eone=eone
    SCF%ekin=ekin
    SCF%ecoul=ecoul
    counter=counter+1

    DEALLOCATE(dum)
  END SUBROUTINE Get_KinCoul

  SUBROUTINE Get_FCKinCoul(Grid,Pot,Orbit,FC,SCF,noalt)
    !  program to calculate Frozencore part of
    !    Kinetic energy and Coulomb Energies from Orbit%wfn
    !  Valence energy defined in paper --
    !    U. von Barth and C. D. Gelatt, Phys. Rev. B 21, 2222(1980)
    !  Exchange terms treated later
    TYPE(GridInfo), INTENT(INOUT) :: Grid
    TYPE(PotentialInfo), INTENT(INOUT) :: Pot
    TYPE(OrbitInfo), INTENT(INOUT) :: Orbit
    TYPE(FCInfo), INTENT(INOUT) :: FC
    TYPE(SCFInfo), INTENT(INOUT) :: SCF
    LOGICAL, OPTIONAL :: noalt

    REAL(8) :: tv,tc,x,y,electrons,vcoul,ccoul,vnucl,rescale,fac
    INTEGER :: icount,i,j,k,l,m,n,io
    REAL(8), ALLOCATABLE :: dum(:)
    REAL(8) :: qcal,small,small0=1.d-6
    INTEGER :: firsttime=0

    n=Grid%n
    small=small0

    ALLOCATE(dum(n),STAT=k)
    IF (k /= 0) THEN
       WRITE(STD_OUT,*) 'Error in Get_FCKinCoul allocation; n,k=',n,k
       STOP
    ENDIF


     write(std_out,*) 'In Get_FCKinCoul ', firsttime
    !update density  and calculated  kinetic energy
    Orbit%den(1:n)=0.d0
    Orbit%tau(1:n)=0.d0
    FC%valeden(1:n)=0.d0     ;   tv=0 ; tc=0; SCF%eone=0
    DO io=1,Orbit%norbit
       IF (Orbit%occ(io).GT.small) THEN
          If(Orbit%iscore(io)) THEN
            If (firsttime==0) then
              If(present(noalt).or.needvtau) then
                  !CALL kinetic(Grid,Orbit%wfn(:,io),Orbit%l(io),x)
                   x=integrator(Grid,Orbit%otau(:,io))
              else
                  CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io),Pot%rv,x)
              endif
              tc=tc+Orbit%occ(io)*x
            Endif
          Else
              FC%valeden(:)=FC%valeden(:)+ Orbit%occ(io)*(Orbit%wfn(:,io)**2)
             IF (diracrelativistic) then
              FC%valeden(:)=FC%valeden(:)+ &
&              Orbit%occ(io)*((Orbit%lwfn(:,io))**2)
             ENDIF         
              If(present(noalt).or.needvtau) then
                  !CALL kinetic(Grid,Orbit%wfn(:,io),Orbit%l(io),x)
                  x=integrator(Grid,Orbit%otau(:,io))
              else
                  CALL altkinetic(Grid,Orbit%wfn(:,io),Orbit%eig(io),Pot%rv,x)
              endif
              tv=tv+Orbit%occ(io)*x
          ENDIF
          Orbit%den(:)=Orbit%den(:)+ Orbit%occ(io)*(Orbit%wfn(:,io)**2)
          IF (diracrelativistic) then
            Orbit%den(:)=Orbit%den(:)+ &
&           Orbit%occ(io)*((Orbit%lwfn(:,io))**2)
          ENDIF         
          CALL taufromwfn(Grid,Orbit%wfn(:,io),Orbit%l(io),Orbit%otau(:,io))
          Orbit%tau(:)=Orbit%tau(:)+ Orbit%occ(io)*Orbit%otau(:,io)
          SCF%eone=SCF%eone+Orbit%occ(io)*Orbit%eig(io)
       ENDIF
    ENDDO

    write(std_out,*) 'electron density integral',integrator(Grid,Orbit%den) 

    write(std_out,*) 'eone = ', SCF%eone

    if (firsttime==0) SCF%corekin=tc
    firsttime=firsttime+1
    SCF%valekin=tv
    SCF%ekin=SCF%corekin+tv
    qcal=integrator(Grid,FC%valeden)
    WRITE(STD_OUT,*) 'qcalval = ', qcal
    Write(STD_OUT,*) 'Core kin ', SCF%corekin,'  Vale kin  ', SCF%valekin

    electrons=FC%zvale
    rescale=electrons/qcal
    FC%valeden(1:n)=FC%valeden(1:n)*rescale
    WRITE(STD_OUT,*) 'rescaled qcalval = ', integrator(Grid,FC%valeden), FC%zvale

    x=FC%zvale
    CALL poisson(Grid,x,FC%valeden,dum,vcoul,v00=y)  !valence-valence

    dum(2:n)=(dum(2:n)*FC%coreden(2:n)+Pot%rvn(2:n)*FC%valeden(2:n))/Grid%r(2:n)
    dum(1)=0.d0
    SCF%valecoul=vcoul+integrator(Grid,dum)

    !Update complete Hartree potential
    qcal=integrator(Grid,Orbit%den)
    WRITE(STD_OUT,*) 'qcal = ', qcal
    rescale=Pot%q/qcal; Orbit%den=Orbit%den*rescale
    write(std_out,*) 'rescaled qcal ', integrator(Grid,Orbit%den),Pot%q
    Call poisson(Grid,Pot%q,Orbit%den,Pot%rvh,x,v00=Pot%v0)
        dum=0
    dum(2:n)=Pot%rvn(2:n)*Orbit%den(2:n)/Grid%r(2:n)
    SCF%estatic=integrator(Grid,dum)+x
    SCF%ecoul=x

    DEALLOCATE(dum)
  END SUBROUTINE Get_FCKinCoul


  SUBROUTINE Get_Nuclearpotential(Grid,Pot)
    TYPE(GridInfo), INTENT(INOUT) :: Grid
    TYPE(PotentialInfo), INTENT(INOUT) :: Pot
!  Various finite nuclear models follow the manuscript of Andrae
!   Physics Reports 336 (2000) 413-525
!    finitenucleusmodel 2,3,4,5 correspond to the options
!     described in that paper while finitenucleusmodel 0 corresponds to
!     Gaussian model originally programmed
!     Note that logarithmic grid is reset to be compatible with
!       nuclear model with approximately NN integration points within
!       finite nucleus
    INTEGER :: n
    REAL(8) :: h,q,v0,v0p
    REAL(8) :: r,RR,r0,a
    INTEGER :: i,j,k
    INTEGER, PARAMETER :: NN=651    ! number of grid points within RR
    REAL(8), PARAMETER :: gridrange=100.d0
    REAL(8), PARAMETER :: bohr=0.529177249d0  !Ang/Bohr from Andrae

    IF (.NOT.finitenucleus) THEN
      !  grid already set      
       DO i=1,Grid%n
          Pot%rvn(i)=-2*Pot%nz
       ENDDO
    ELSE
       write(std_out,*) 'Finite nucleus model  -- readjusting integration grid'
       a=bohr*1.d-5*(0.57d0+0.836*   &
        &     (-1.168d0+Pot%nz*(2.163d0+Pot%nz*0.004467d0)))
           !  From Eqs. A.3 and 51 in Andrae paper
       write(std_out,*) 'a parameter calculated to be', a
       call destroygrid(Grid)
       SELECT CASE(Pot%finitenucleusmodel)
          CASE DEFAULT
            write(std_out,*) 'Error in finitenucleusmodel',Pot%finitenucleusmodel
            write(std_out,*) ' Exiting '
            Stop
          CASE(0)      
            write(std_out,*) 'Original Gaussian model'      
            RR=Pot%nz
            RR=2.9d-5*(RR**0.3333333333333333333333333d0)
            h=log(FLOAT(NN))/(NN-1)
            r0=RR/(NN-1)
            write(std_out,*) 'calling InitGrid with h, r0 =',h,r0
            Call InitGrid(Grid,h,gridrange,r0=r0)
            write(std_out,*) 'New Grid ', Grid%n
            Call DestroyPot(Pot)
            Call InitPot(Pot,Grid%n)
            DO i=1,Grid%n
              Pot%rvn(i)=-2*Pot%nz*derf(Grid%r(i)/RR)
            ENDDO
              Pot%Nv0=-2*Pot%nz*sqrt(4.d0/pi)
              Pot%Nv0p=0.d0
          CASE(2)      
            write(std_out,*) 'Model 2 -- Breit'      
            RR=sqrt(2.d0)*a
            h=log(FLOAT(NN))/(NN-1)
            r0=RR/(NN-1)
            write(std_out,*) 'calling InitGrid with h, r0 =',h,r0
            Call InitGrid(Grid,h,gridrange,r0=r0)
            write(std_out,*) 'New Grid ', Grid%n
            Call DestroyPot(Pot)
            Call InitPot(Pot,Grid%n)
            DO i=1,Grid%n
               if (Grid%r(i)<RR) then
                 Pot%rvn(i)=-2*Pot%nz*Grid%r(i)*(2.d0-Grid%r(i)/RR)/RR
               else  
                 Pot%rvn(i)=-2*Pot%nz
               endif  
            ENDDO
              Pot%Nv0=-2*Pot%nz*2.0d0/RR
              Pot%Nv0p=2*Pot%nz/(RR**2)
          CASE(3)      
            write(std_out,*) 'Model 3 -- uniform'      
            RR=sqrt(5.d0/3.d0)*a
            h=log(FLOAT(NN))/(NN-1)
            r0=RR/(NN-1)
            write(std_out,*) 'calling InitGrid with h, r0 =',h,r0
            Call InitGrid(Grid,h,gridrange,r0=r0)
            write(std_out,*) 'New Grid ', Grid%n
            Call DestroyPot(Pot)
            Call InitPot(Pot,Grid%n)
            DO i=1,Grid%n
               if (Grid%r(i)<RR) then
                 Pot%rvn(i)=-3*Pot%nz*Grid%r(i)*&
                    &     (1.d0-(Grid%r(i)/RR)**2/3)/RR
               else  
                 Pot%rvn(i)=-2*Pot%nz
               endif  
            ENDDO
              Pot%Nv0=-3*Pot%nz/RR
              Pot%Nv0p=0.d0
          CASE(4)      
            write(std_out,*) 'Model 4 -- exponential'      
            RR=sqrt(1.d0/12.d0)*a
            h=log(FLOAT(NN))/(NN-1)
            r0=RR/(NN-1)
            write(std_out,*) 'calling InitGrid with h, r0 =',h,r0
            Call InitGrid(Grid,h,gridrange,r0=r0)
            write(std_out,*) 'New Grid ', Grid%n
            Call DestroyPot(Pot)
            Call InitPot(Pot,Grid%n)
            DO i=1,Grid%n
             Pot%rvn(i)=-2*Pot%nz*   &
               &  (1.d0-exp(-grid%r(i)/RR)*(1.d0+0.5d0*Grid%r(i)/RR))
            ENDDO
              Pot%Nv0=-Pot%nz/RR
              Pot%Nv0p=0.d0
          CASE(5)      
            write(std_out,*) 'Model 5 -- Gaussian'      
            RR=sqrt(2.d0/3.d0)*a
            h=log(FLOAT(NN))/(NN-1)
            r0=RR/(NN-1)
            write(std_out,*) 'calling InitGrid with h, r0 =',h,r0
            Call InitGrid(Grid,h,gridrange,r0=r0)
            write(std_out,*) 'New Grid ', Grid%n
            Call DestroyPot(Pot)
            Call InitPot(Pot,Grid%n)
            DO i=1,Grid%n
             Pot%rvn(i)=-2*Pot%nz*erf(Grid%r(i)/RR)
            ENDDO
              Pot%Nv0=-2*Pot%nz/RR*(sqrt(4.d0/pi))
              Pot%Nv0p=0.d0
    END SELECT     
    open(7,file='nuclearpot.dat',form='formatted')
       write(7,*) '#  ', Pot%Nv0, Pot%Nv0p
    do i=1,Grid%n
       write(7,'(1P2E16.7)') Grid%r(i),Pot%rvn(i)
    enddo
    close(7)   
    ENDIF
  END SUBROUTINE Get_Nuclearpotential

  SUBROUTINE ORTHONORMALIZE(Grid,Orbit)
    TYPE(GridInfo) ,INTENT(IN):: Grid
    TYPE(OrbitInfo) ,INTENT(INOUT):: Orbit

    Integer :: l,lmax,many,io,n
    REAL(8), ALLOCATABLE :: wfn(:,:)

    lmax=MAXVAL(Orbit%l)
    n=Grid%n
    allocate (wfn(n,Orbit%norbit))

    lmax=MAXVAL(Orbit%l)
    DO l=0,lmax
       wfn=0;many=0
       DO io=1,Orbit%norbit
          IF(Orbit%l(io)==l) THEN
            if (frozencorecalculation.AND.Orbit%iscore(io)) then
               many=many+1
               wfn(:,many)=Orbit%wfn(:,io)
            else
               CALL gramschmidt(Grid,many,wfn,Orbit%wfn(:,io))
               many=many+1
               Orbit%wfn(:,io)=wfn(:,many)
               CALL ADJUSTSIGN(Orbit%wfn(:,io),3)
            endif
          ENDIF
       ENDDO
    ENDDO

    deallocate (wfn)

  END SUBROUTINE ORTHONORMALIZE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Subroutine to input f(:) and return f(:) if f(index)>=0 or
!                                return -f(:) if f(index)<0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ADJUSTSIGN(f,index)
     REAL(8),INTENT(INOUT) :: f(:)
     INTEGER,INTENT(IN) :: index

     IF (f(index)<0.d0) then
        f=-f
     ENDIF
  END SUBROUTINE ADJUSTSIGN

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Adjustnodes(Grid,Orbit,whichone)
!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE Adjustnodes(Grid,Orbit,whichone)
!    TYPE(GridInfo),INTENT(IN)::Grid
!    TYPE(OrbitInfo),INTENT(INOUT)::Orbit
!    INTEGER, INTENT(IN), OPTIONAL :: whichone
!
!    INTEGER :: io, i,j,n,num,nodes
!
!    n=Grid%n
!
!    if (PRESENT(whichone)) then
!       io=whichone
!       num=0;j=n+1
!       nodes=Orbit%np(io)-Orbit%l(io)-1
!       do i=3,n
!          if(Orbit%wfn(i-1,io)*Orbit%wfn(i,io)<0.d0) then
!            num=num+1
!            if (num>nodes) then
!               j=i
!               exit
!            endif
!          endif
!       enddo
!       If (j<n+1) then
!            Orbit%wfn(j:n,io)=0.d0
!            write(std_out,*) 'Adjusted orbit ', io,'  at r = ', &
!&                  Grid%r(j), Orbit%wfn(j-1,io)
!       endif
!     else
!       do io=1,Orbit%norbit
!          num=0;j=n+1
!          nodes=Orbit%np(io)-Orbit%l(io)-1
!          do i=3,n
!             if(Orbit%wfn(i-1,io)*Orbit%wfn(i,io)<0.d0) then
!               num=num+1
!               if (num>nodes) then
!                  j=i
!                  exit
!               endif
!             endif
!          enddo
!          If (j<n+1) Orbit%wfn(j:n,io)=0.d0
!       enddo
!     endif
!  END SUBROUTINE Adjustnodes
!

END  MODULE general_mod

