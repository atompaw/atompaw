MODULE atomdata

  IMPLICIT NONE

  TYPE OrbitInfo
     CHARACTER(132) :: exctype
     INTEGER :: nps, npp, npd ,npf, npg, norbit
     INTEGER, POINTER :: np(:),l(:)
     REAL(8), POINTER :: eig(:),occ(:),wfn(:,:)
     REAL(8), POINTER :: lqp(:,:)     ! only used for HF
     REAL(8), POINTER :: X(:,:)       ! identical to HF%SumY(:,:)
     LOGICAL , POINTER :: iscore(:)
     REAL(8),POINTER :: den(:)
  END TYPE OrbitInfo

  TYPE FCinfo
     REAL(8), POINTER :: coreden(:),valeden(:)
     REAL(8) :: zvale,zcore
  END TYPE FCinfo

  TYPE PotentialInfo
     CHARACTER(2) :: sym
     INTEGER :: nz     !  nz is nuclear charge     
     REAL(8) :: zz        !  zz=nz is nuclear charge
     REAL(8) :: q,v0,v0p  !  q is total electron charge
     !  v0,v0p are potential value and deriv at r=0
     REAL(8) , POINTER :: rv(:),rvn(:),rvh(:),rvx(:)
     !  rv(n) is  veff * r
     !  rvh is hartree potential for den
     !  rvn is nuclear potential
     !  rvx is exchange-correlation potential
  END TYPE PotentialInfo

  TYPE SCFInfo
     INTEGER :: iter
     REAL(8) :: delta,eone,ekin,estatic,ecoul,eexc,oepcs,etot
     REAL(8) :: valekin,valecoul,valeexc,corekin,evale ! used in frozencore only
  END TYPE SCFInfo

  LOGICAL :: frozencorecalculation
  LOGICAL :: setupfrozencore
  LOGICAL :: scalarrelativistic
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
    Orbit%norbit=norbit;Orbit%exctype=TRIM(exctype)
    Orbit%nps=0;Orbit%npp=0;Orbit%npd=0;Orbit%npf=0;Orbit%npg=0
    ALLOCATE(Orbit%np(norbit),Orbit%l(norbit),Orbit%eig(norbit),&
&            Orbit%occ(norbit),Orbit%iscore(norbit),&
&            stat=ok)
    IF (ok/=0) STOP 'Error in allocation of nl, l, occ...'
    Orbit%iscore=.false.
    Orbit%np=0;Orbit%l=0
    Orbit%eig=0.d0;Orbit%occ=0.d0
    ALLOCATE(Orbit%wfn(n,norbit),Orbit%den(n),stat=ok)
    IF (ok/=0) STOP 'Error in allocation of wfn, den...'
    Orbit%wfn=0.d0;Orbit%den=0.d0
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
    IF (ASSOCIATED(Orbit%iscore)) DEALLOCATE(Orbit%iscore)
    IF (ASSOCIATED(Orbit%eig)) DEALLOCATE(Orbit%eig)
    IF (ASSOCIATED(Orbit%occ)) DEALLOCATE(Orbit%occ)
    IF (ASSOCIATED(Orbit%wfn)) DEALLOCATE(Orbit%wfn)
    IF (ASSOCIATED(Orbit%den)) DEALLOCATE(Orbit%den)
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
    COrbit%iscore(1:SOrbit%norbit)=SOrbit%iscore(1:SOrbit%norbit)
    COrbit%den=SOrbit%den
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
    FC%zvale=0.d0;FC%zcore=0.d0
    ALLOCATE(FC%coreden(n),FC%valeden(n),stat=ok)
    IF (ok/=0) STOP 'Error in allocation of coreden, valeden,...'
    FC%coreden=0.d0;FC%valeden=0.d0
  END SUBROUTINE InitFC

  SUBROUTINE DestroyFC(FC)
    TYPE (FCInfo), INTENT(INOUT) :: FC
    IF (ASSOCIATED(FC%valeden)) DEALLOCATE(FC%valeden)
    IF (ASSOCIATED(FC%coreden)) DEALLOCATE(FC%coreden)
  END SUBROUTINE DestroyFC

  SUBROUTINE InitPot(Pot,n)
    INTEGER, INTENT(IN) :: n
    TYPE (PotentialInfo), INTENT(INOUT) :: Pot
    INTEGER :: ok
!   Pot%sym="";Pot%nz=0;Pot%zz=0.d0;Pot%q=0.d0;Pot%v0=0.d0;Pot%v0p=0.d0
    ALLOCATE(Pot%rv(n),Pot%rvn(n),Pot%rvh(n),Pot%rvx(n),stat=ok)
    IF (ok/=0) STOP 'Error in allocation of Pot%rv, Pot%rvh...'
    Pot%rv=0.d0;Pot%rvn=0.d0;Pot%rvh=0.d0;Pot%rvx=0.d0
  END SUBROUTINE InitPot

  SUBROUTINE DestroyPot(Pot)
    TYPE (PotentialInfo), INTENT(INOUT) :: Pot
    IF (ASSOCIATED(Pot%rv)) DEALLOCATE(Pot%rv)
    IF (ASSOCIATED(Pot%rvn)) DEALLOCATE(Pot%rvn)
    IF (ASSOCIATED(Pot%rvh)) DEALLOCATE(Pot%rvh)
    IF (ASSOCIATED(Pot%rvx)) DEALLOCATE(Pot%rvx)
  END SUBROUTINE DestroyPot

!!!!!!!!!!!!!!!!!!!!!!!!!
!  CopyPot(source,copy)
!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CopyPot(SPot,CPot)
    TYPE(PotentialInfo),INTENT(IN) :: SPot
    TYPE(PotentialInfo),INTENT(INOUT) :: CPot
    INTEGER :: n
    CPot%nz=SPot%nz
    CPot%zz=SPot%zz
    CPot%sym=SPot%sym
    CPot%q=SPot%q
    CPot%v0=SPot%v0
    CPot%v0p=SPot%v0p
    n=SIZE(SPot%rv,1)
    ALLOCATE(CPot%rv(n),CPot%rvn(n),CPot%rvh(n),CPot%rvx(n))
    CPot%rv(1:n)=SPot%rv(1:n)
    CPot%rvn(1:n)=SPot%rvn(1:n)
    CPot%rvh(1:n)=SPot%rvh(1:n)
    CPot%rvx(1:n)=SPot%rvx(1:n)
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
