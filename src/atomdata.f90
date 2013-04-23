MODULE atomdata

  TYPE OrbitInfo
     CHARACTER(132) :: exctype
     INTEGER :: nps, npp, npd ,npf, npg, norbit
     INTEGER, POINTER :: np(:),l(:)
     REAL(8), POINTER :: eig(:),occ(:),wfn(:,:)
     REAL(8), POINTER :: lqp(:,:)     ! only used for HF
     REAL(8), POINTER :: X(:,:)   ! identical to HF%SumY(:,:)
     LOGICAL , POINTER :: iscore(:)
     REAL(8),POINTER::den(:)
  END TYPE OrbitInfo

  TYPE SCFInfo
     INTEGER :: iter
     REAL(8) :: delta,eone,ekin,estatic,ecoul,eexc,oepcs,etot
     REAL(8) :: valekin,valecoul,valeexc,corekin,evale ! used in frozencore only
  END TYPE SCFInfo

  TYPE FCinfo
     REAL(8), POINTER :: coreden(:),valeden(:)
     REAL(8) :: zvale,zcore
  END TYPE FCinfo

  TYPE PotentialInfo
     INTEGER :: nz       !  nz is nuclear charge
     CHARACTER(2) :: sym
     REAL(8) :: q,v0,v0p  !  q is total electron charge
     !  v0,v0p are potential value and deriv at r=0
     REAL(8) , POINTER :: rv(:),rvn(:),rvh(:),rvx(:)
     !  rv(n) is  veff * r
     !  rvh is hartree potential for den
     !  rvn is nuclear potential
     !  rvx is exchange-correlation potential
  END TYPE PotentialInfo


  LOGICAL :: frozencorecalculation
  LOGICAL :: setupfrozencore
  LOGICAL :: scalarrelativistic
  LOGICAL :: finitenucleus
  LOGICAL :: gaussianshapefunction,besselshapefunction
  LOGICAL :: ColleSalvetti
  LOGICAL :: HFpostprocess
  LOGICAL :: localizedcoreexchange
  LOGICAL :: exploremode


CONTAINS
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

  SUBROUTINE DestroyFC(FC)
    TYPE (FCInfo), INTENT(INOUT) :: FC

       IF (ASSOCIATED(FC%valeden)) DEALLOCATE(FC%valeden)
       IF (ASSOCIATED(FC%coreden)) DEALLOCATE(FC%coreden)

  END SUBROUTINE DestroyFC

  SUBROUTINE DestroyPot(Pot)
    TYPE (PotentialInfo), INTENT(INOUT) :: Pot

       IF (ASSOCIATED(Pot%rv)) DEALLOCATE(Pot%rv)
       IF (ASSOCIATED(Pot%rvn)) DEALLOCATE(Pot%rvn)
       IF (ASSOCIATED(Pot%rvh)) DEALLOCATE(Pot%rvh)
       IF (ASSOCIATED(Pot%rvx)) DEALLOCATE(Pot%rvx)

  END SUBROUTINE DestroyPot



END MODULE atomdata
