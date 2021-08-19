!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module constains the following active subroutines:
!   InitOrbit, DestroyOrbit, CopyOrbit, InitFC, DestroyFC,
!      InitPot, DestroyPot, CopyPot, InitSCF, CopySCF
!
!  Note that all energies including tau are in Rydberg units
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE atomdata

  IMPLICIT NONE

  TYPE OrbitInfo
     CHARACTER(132) :: exctype
     INTEGER :: nps, npp, npd ,npf, npg, norbit
     INTEGER, POINTER :: np(:) => null()
     INTEGER, POINTER :: l(:) => null()
     INTEGER, POINTER :: kappa(:) => null()
     REAL(8), POINTER :: eig(:) => null()
     REAL(8), POINTER :: occ(:) => null()
     REAL(8), POINTER :: wfn(:,:) => null()
     REAL(8), POINTER :: lwfn(:,:) => null()
     REAL(8), POINTER :: otau(:,:) => null() ! kinetic energy density for orbital
     REAL(8), POINTER :: lqp(:,:) => null()  ! only used for HF
     REAL(8), POINTER :: X(:,:) => null()    ! identical to HF%SumY(:,:)
     LOGICAL, POINTER :: iscore(:) => null()
     REAL(8),POINTER :: den(:) => null() ! accumulated over states
     REAL(8),POINTER :: tau(:) => null() ! accumulated over states
     REAL(8),POINTER :: deltatau(:) => null() !tau-tauW   (tauW==Weizsacker)
                                      ! note this is evaluated with 1/4pir^2
  END TYPE OrbitInfo

  TYPE FCinfo
     REAL(8), POINTER :: coreden(:) => null()
     REAL(8), POINTER :: valeden(:) => null()
     REAL(8), POINTER :: coretau(:) => null()
     REAL(8), POINTER :: valetau(:) => null()
     REAL(8) :: zvale,zcore
  END TYPE FCinfo

  TYPE PotentialInfo
     CHARACTER(2) :: sym
     INTEGER :: nz     !  nz is nuclear charge     
     REAL(8) :: zz        !  zz=nz is nuclear charge
     REAL(8) :: q,v0,v0p  !  q is total electron charge
     !  v0,v0p are potential value and deriv at r=0
     REAL(8) :: Nv0,Nv0p    !  finite nucleus value and deriv at 0
     REAL(8), POINTER :: rv(:)  => null()
     REAL(8), POINTER :: rvn(:) => null()
     REAL(8), POINTER :: rvh(:) => null()
     REAL(8), POINTER :: rvx(:) => null()
     !  rv(n) is  veff * r
     !  rvh is hartree potential for den
     !  rvn is nuclear potential
     !  rvx is exchange-correlation potential
     REAL(8), POINTER :: vtau(:) => null() !for meta-gga
     INTEGER :: finitenucleusmodel
     ! Based on models 2, 3, 4, 5 discussed by Dirk Anrae ,
     !   Physics Reports 336 (2000) 413-525
     !    default is 0 for previous Gaussian model
     !    for finitenucleusmodel<0, finite nucleus is false
  END TYPE PotentialInfo

  TYPE SCFInfo
     INTEGER :: iter
     REAL(8) :: delta,eone,ekin,estatic,ecoul,eexc,oepcs,etot
     REAL(8) :: valekin,valecoul,valeexc,corekin,evale ! used in frozencore only
  END TYPE SCFInfo

  LOGICAL :: frozencorecalculation
  LOGICAL :: setupfrozencore
  LOGICAL :: scalarrelativistic
  LOGICAL :: diracrelativistic
  LOGICAL :: shapetcore
  LOGICAL :: needvtau
  LOGICAL :: usespline
  LOGICAL :: BDsolve
  LOGICAL :: finitenucleus
  LOGICAL :: gaussianshapefunction,besselshapefunction
  LOGICAL :: ColleSalvetti
  LOGICAL :: HFpostprocess
  LOGICAL :: localizedcoreexchange
  LOGICAL :: exploremode
  INTEGER :: nlogderiv
  REAL(8) :: minlogderiv,maxlogderiv


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Subroutine InitOrbit  -- used in CopyOrbit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE InitOrbit(Orbit,norbit,n,exctype)
    TYPE (OrbitInfo), INTENT(INOUT) :: Orbit
    INTEGER, INTENT(IN) :: n,norbit
    CHARACTER(*),INTENT(IN) :: exctype
    INTEGER :: ok
    CALL DestroyOrbit(Orbit)
    Orbit%norbit=norbit;Orbit%exctype=TRIM(exctype)
    Orbit%nps=0;Orbit%npp=0;Orbit%npd=0;Orbit%npf=0;Orbit%npg=0
    ALLOCATE(Orbit%np(norbit),Orbit%l(norbit),Orbit%eig(norbit),&
&            Orbit%occ(norbit),Orbit%iscore(norbit),&
&            stat=ok)
    IF (ok/=0) STOP 'Error in allocation of nl, l, occ...'
    Orbit%iscore=.false.
    Orbit%np=0;Orbit%l=0
    Orbit%eig=0.d0;Orbit%occ=0.d0
    ALLOCATE(Orbit%wfn(n,norbit),Orbit%otau(n,norbit), &
&             Orbit%den(n),Orbit%tau(n),Orbit%deltatau(n),stat=ok)
    IF (ok/=0) STOP 'Error in allocation of wfn, den...'
    Orbit%wfn=0.d0;Orbit%den=0.d0;Orbit%tau=0.d0;Orbit%otau=0.d0
    Orbit%deltatau=0.d0
    If (diracrelativistic) then
       ALLOCATE(Orbit%lwfn(n,norbit),Orbit%kappa(norbit),stat=ok)
       IF (ok/=0) STOP 'Error in allocation of lwfn,kappa'
       Orbit%lwfn=0.d0
       Orbit%kappa=0.d0
    else
       NULLIFY(Orbit%lwfn,Orbit%kappa)
    Endif
    If (exctype == "HF".or.exctype == "EXXKLI") then
       ALLOCATE(Orbit%lqp(norbit,norbit),Orbit%X(n,norbit),stat=ok)
       IF (ok/=0) STOP 'Error in allocation of lqp, X...'
    ELSE
       NULLIFY(Orbit%lqp,Orbit%X)
    ENDIF
  END SUBROUTINE InitOrbit

  SUBROUTINE DestroyOrbit(Orbit)
    TYPE (OrbitInfo), INTENT(INOUT) :: Orbit
    IF (ASSOCIATED(Orbit%np)) DEALLOCATE(Orbit%np)
    IF (ASSOCIATED(Orbit%l)) DEALLOCATE(Orbit%l)
    IF (ASSOCIATED(Orbit%kappa)) DEALLOCATE(Orbit%kappa)
    IF (ASSOCIATED(Orbit%iscore)) DEALLOCATE(Orbit%iscore)
    IF (ASSOCIATED(Orbit%eig)) DEALLOCATE(Orbit%eig)
    IF (ASSOCIATED(Orbit%occ)) DEALLOCATE(Orbit%occ)
    IF (ASSOCIATED(Orbit%wfn)) DEALLOCATE(Orbit%wfn)
    IF (ASSOCIATED(Orbit%otau)) DEALLOCATE(Orbit%otau)
    IF (ASSOCIATED(Orbit%lwfn)) DEALLOCATE(Orbit%lwfn)
    IF (ASSOCIATED(Orbit%den)) DEALLOCATE(Orbit%den)
    IF (ASSOCIATED(Orbit%tau)) DEALLOCATE(Orbit%tau)
    IF (ASSOCIATED(Orbit%deltatau)) DEALLOCATE(Orbit%deltatau)
    IF (ASSOCIATED(Orbit%lqp)) DEALLOCATE(Orbit%lqp)
    IF (ASSOCIATED(Orbit%X)) DEALLOCATE(Orbit%X)
  END SUBROUTINE DestroyOrbit

!!!!!!!!!!!!!!!!!!!!!!!!!
!  CopyOrbit(source,copy)
!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CopyOrbit(SOrbit,COrbit)
    TYPE(OrbitInfo),INTENT(INOUT)::SOrbit
    TYPE(OrbitInfo),INTENT(INOUT)::COrbit
    INTEGER::n
    n=SIZE(SOrbit%den,1)
    CALL InitOrbit(COrbit,SOrbit%norbit,n,SOrbit%exctype)
    COrbit%nps=SOrbit%nps
    COrbit%npp=SOrbit%npp
    COrbit%npd=SOrbit%npd
    COrbit%npf=SOrbit%npf
    COrbit%npg=SOrbit%npg
    COrbit%np(1:SOrbit%norbit)=SOrbit%np(1:SOrbit%norbit)
    COrbit%l(1:SOrbit%norbit)=SOrbit%l(1:SOrbit%norbit)
    COrbit%eig(1:SOrbit%norbit)=SOrbit%eig(1:SOrbit%norbit)
    COrbit%occ(1:SOrbit%norbit)=SOrbit%occ(1:SOrbit%norbit)
    COrbit%wfn(:,1:SOrbit%norbit)=SOrbit%wfn(:,1:SOrbit%norbit)
    COrbit%otau(:,1:SOrbit%norbit)=SOrbit%otau(:,1:SOrbit%norbit)
    COrbit%iscore(1:SOrbit%norbit)=SOrbit%iscore(1:SOrbit%norbit)
    COrbit%den=SOrbit%den
    COrbit%tau=SOrbit%tau
    COrbit%deltatau=SOrbit%deltatau
    If (diracrelativistic) then
       COrbit%lwfn(:,1:SOrbit%norbit)=SOrbit%lwfn(:,1:SOrbit%norbit)
       COrbit%kappa(1:SOrbit%norbit)=SOrbit%kappa(1:SOrbit%norbit)
    Endif

    IF (SOrbit%exctype == "HF".or.SOrbit%exctype == "EXXKLI") then
        COrbit%X(:,1:SOrbit%norbit)=SOrbit%X(:,1:SOrbit%norbit)
        COrbit%lqp(1:SOrbit%norbit,1:SOrbit%norbit)= &
&                SOrbit%lqp(1:SOrbit%norbit,1:SOrbit%norbit)
    ENDIF
  END SUBROUTINE CopyOrbit

  SUBROUTINE InitFC(FC,n)
    INTEGER, INTENT(IN) :: n
    TYPE (FCInfo), INTENT(INOUT) :: FC
    INTEGER :: ok
    CALL DestroyFC(FC)
    FC%zvale=0.d0;FC%zcore=0.d0
    ALLOCATE(FC%coreden(n),FC%valeden(n),stat=ok)
    IF (ok/=0) STOP 'Error in allocation of coreden, valeden,...'
    FC%coreden=0.d0;FC%valeden=0.d0
    ALLOCATE(FC%coretau(n),FC%valetau(n),stat=ok)
    IF (ok/=0) STOP 'Error in allocation of coretau, valetau,...'
    FC%coretau=0.d0;FC%valetau=0.d0
  END SUBROUTINE InitFC

  SUBROUTINE DestroyFC(FC)
    TYPE (FCInfo), INTENT(INOUT) :: FC
    IF (ASSOCIATED(FC%valeden)) DEALLOCATE(FC%valeden)
    IF (ASSOCIATED(FC%coreden)) DEALLOCATE(FC%coreden)
    IF (ASSOCIATED(FC%valetau)) DEALLOCATE(FC%valetau)
    IF (ASSOCIATED(FC%coretau)) DEALLOCATE(FC%coretau)
  END SUBROUTINE DestroyFC

!!!!!!!!!!!!!!!!!!!!!!!!!
!  CopyFC(source,copy)
!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CopyFC(SFC,CFC)
    TYPE(FCInfo),INTENT(IN) :: SFC
    TYPE(FCInfo),INTENT(INOUT) :: CFC
    INTEGER :: n
    n=SIZE(SFC%coreden,1)
    CALL InitFC(CFC,n)
    CFC%zvale=SFC%zvale
    CFC%zcore=SFC%zcore
    CFC%coreden(1:n)=SFC%coreden(1:n)
    CFC%valeden(1:n)=SFC%valeden(1:n)
    CFC%coretau(1:n)=SFC%coretau(1:n)
    CFC%valetau(1:n)=SFC%valetau(1:n)
  END SUBROUTINE CopyFC

  SUBROUTINE InitPot(Pot,n)
    INTEGER, INTENT(IN) :: n
    TYPE (PotentialInfo), INTENT(INOUT) :: Pot
    INTEGER :: ok
    CALL DestroyPot(Pot)
!   Pot%sym="";Pot%nz=0;Pot%zz=0.d0;Pot%q=0.d0;Pot%v0=0.d0;Pot%v0p=0.d0
    ALLOCATE(Pot%rv(n),Pot%rvn(n),Pot%rvh(n),Pot%rvx(n),Pot%vtau(n),stat=ok)
    IF (ok/=0) STOP 'Error in allocation of Pot%rv, Pot%rvh...'
    Pot%rv=0.d0;Pot%rvn=0.d0;Pot%rvh=0.d0;Pot%rvx=0.d0;Pot%vtau=0.d0
  END SUBROUTINE InitPot

  SUBROUTINE DestroyPot(Pot)
    TYPE (PotentialInfo), INTENT(INOUT) :: Pot
    IF (ASSOCIATED(Pot%rv)) DEALLOCATE(Pot%rv)
    IF (ASSOCIATED(Pot%rvn)) DEALLOCATE(Pot%rvn)
    IF (ASSOCIATED(Pot%rvh)) DEALLOCATE(Pot%rvh)
    IF (ASSOCIATED(Pot%rvx)) DEALLOCATE(Pot%rvx)
    IF (ASSOCIATED(Pot%vtau)) DEALLOCATE(Pot%vtau)
  END SUBROUTINE DestroyPot

!!!!!!!!!!!!!!!!!!!!!!!!!
!  CopyPot(source,copy)
!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CopyPot(SPot,CPot)
    TYPE(PotentialInfo),INTENT(IN) :: SPot
    TYPE(PotentialInfo),INTENT(INOUT) :: CPot
    INTEGER :: n
    n=SIZE(SPot%rv,1)
    CALL InitPot(CPot,n)
    CPot%nz=SPot%nz
    CPot%zz=SPot%zz
    CPot%sym=SPot%sym
    CPot%q=SPot%q
    CPot%v0=SPot%v0
    CPot%v0p=SPot%v0p
    CPot%finitenucleusmodel=SPot%finitenucleusmodel
    CPot%Nv0=SPot%Nv0
    CPot%Nv0p=SPot%Nv0p
    CPot%rv(1:n)=SPot%rv(1:n)
    CPot%rvn(1:n)=SPot%rvn(1:n)
    CPot%rvh(1:n)=SPot%rvh(1:n)
    CPot%rvx(1:n)=SPot%rvx(1:n)
    CPot%vtau(1:n)=SPot%vtau(1:n)
  END SUBROUTINE CopyPot

  SUBROUTINE InitSCF(SCF)
    TYPE(SCFInfo),INTENT(INOUT)::SCF
    SCF%iter=0
    SCF%delta=0.d0;SCF%eone=0.d0;SCF%ekin=0.d0;SCF%estatic=0.d0
    SCF%ecoul=0.d0;SCF%eexc=0.d0;SCF%oepcs=0.d0;SCF%etot=0.d0
    SCF%valekin=0.d0;SCF%valecoul=0.d0;SCF%valeexc=0.d0
    SCF%corekin=0.d0;SCF%evale=0.d0
  END SUBROUTINE InitSCF

!!!!!!!!!!!!!!!!!!!!!!!!!
!  CopySCF(source,copy)
!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CopySCF(SSCF,CSCF)
    TYPE(SCFInfo),INTENT(IN)::SSCF
    TYPE(SCFInfo),INTENT(INOUT)::CSCF
    CALL InitSCF(CSCF)
    CSCF%iter=SSCF%iter
    CSCF%delta=SSCF%delta
    CSCF%eone=SSCF%eone
    CSCF%ekin=SSCF%ekin
    CSCF%ecoul=SSCF%ecoul
    CSCF%estatic=SSCF%estatic
    CSCF%eexc=SSCF%eexc
    CSCF%etot=SSCF%etot
    CSCF%valekin=SSCF%valekin
    CSCF%valecoul=SSCF%valecoul
    CSCF%valeexc=SSCF%valeexc
    CSCF%corekin=SSCF%corekin
    CSCF%evale=SSCF%evale
  END SUBROUTINE CopySCF

END MODULE atomdata
