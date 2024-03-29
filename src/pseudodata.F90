!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!        InitPAW, DestroyPAW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE pseudodata

  Use io_tools
  Use gridmod
  Use atomdata

  IMPLICIT NONE

  TYPE  Pseudoinfo
     CHARACTER(132) :: exctype
     INTEGER  :: lmax,irc,irc_shap,irc_vloc,irc_core,coretailpoints,mesh_size
     INTEGER  :: ivale,itau,ivion
     CHARACTER(15) :: orthogonalization_scheme
     CHARACTER(132) :: Vloc_description
     CHARACTER(132) :: Proj_description
     CHARACTER(132) :: Comp_description
     LOGICAL :: multi_rc,poscorenhat
     REAL(8) :: rc,rc_shap,rc_vloc,rc_core,energyoflmax,gausslength
     REAL(8), POINTER :: rcio(:)
     REAL(8), POINTER :: vloc(:),abinitvloc(:),abinitnohat(:)
     REAL(8), POINTER :: rveff(:),AErefrv(:),rvx(:),trvx(:)
     REAL(8), POINTER :: Ktvtau(:),Krveff(:),Kunscreen(:) ! Kresse form
     REAL(8), POINTER :: projshape(:),hatshape(:),hatden(:),hatpot(:)
     REAL(8), POINTER :: den(:),tden(:),core(:),tcore(:),nhatv(:)
     REAL(8), POINTER :: coretau(:),tcoretau(:)
     REAL(8), POINTER :: valetau(:),tvaletau(:)
     REAL(8), POINTER :: vtau(:),tvtau(:)
     INTEGER :: nbase,ncoreshell
     INTEGER, POINTER :: np(:),l(:),nodes(:),kappa(:)
     INTEGER, POINTER :: rng(:)       ! rng particularly of continuum states
     CHARACTER(8), POINTER :: label(:)
     REAL(8), POINTER :: phi(:,:),tphi(:,:),tp(:,:) ! before orthog
     REAL(8), POINTER :: ophi(:,:),otphi(:,:),otp(:,:) ! after orthog
     REAL(8), POINTER :: Kop(:,:)    ! for storing K|phi>
     REAL(8), POINTER :: eig(:),occ(:),ck(:),vrc(:)
     REAL(8), POINTER :: oij(:,:),dij(:,:),wij(:,:)
     !********** modified parameters for use with KS and HF
     REAL(8), POINTER :: rVf(:),rtVf(:),g(:,:)
     REAL(8), POINTER :: Kij(:,:),Vfij(:,:),mLij(:,:,:),DR(:,:,:,:,:)
     REAL(8), POINTER :: DRVC(:,:,:),TXVC(:,:)   ! now output for DFT also
     REAL(8) :: lambshielding    
     REAL(8) :: XCORECORE     ! output for DFT
     INTEGER, POINTER :: valencemap(:)   ! valencemap({occ. states})={basis}
     Type(OrbitInfo), POINTER :: OCCwfn
     Type(OrbitInfo), POINTER :: TOCCwfn
     REAL(8) :: tkin,tion,tvale,txc,Ea,Etotal,Eaion,Eaionhat,Eaxc
     REAL(8) :: VlocCoef,VlocRad
     !***********for HF only
     REAL(8), POINTER :: lmbd(:,:) !(Eq. 72) lmbd({occ. states},{basis states})
     REAL(8), POINTER :: DRC(:,:,:,:),mLic(:,:,:)
     REAL(8), POINTER :: DRCC(:,:,:,:),DRCjkl(:,:,:,:,:),mLcc(:,:,:),Dcj(:,:)
     REAL(8) :: coretol
  END  TYPE Pseudoinfo

  CONTAINS

    SUBROUTINE InitPAW(PAW,Grid,Orbit)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(OrbitInfo), INTENT(IN) :: Orbit
      Type(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER :: io,l,n,mxbase,nbase,ok
!     Initialize logical variables
      PAW%multi_rc=.false.
      PAW%poscorenhat=.true.
      CALL DestroyPAW(PAW)
!     Compute initial size of basis
      n=Grid%n
      nbase=0
      DO l=0,PAW%lmax
         DO io=1,Orbit%norbit    ! cycle through all configurations
            IF (Orbit%l(io).EQ.l.AND.(.NOT.Orbit%iscore(io))) THEN
               nbase=nbase+1
            ENDIF
         ENDDO
      ENDDO
      mxbase=nbase+5*max(1,PAW%lmax) !Estimate excess
      PAW%nbase=nbase
      WRITE(STD_OUT,*) 'Found ', nbase,' valence basis functions '
      WRITE(STD_OUT,*) 'Allocating for ', mxbase, ' total basis functions'
      ALLOCATE(PAW%projshape(n),PAW%hatden(n),PAW%hatpot(n),&
&        PAW%hatshape(n),PAW%vloc(n),PAW%rveff(n),PAW%abinitvloc(n),&
&        PAW%abinitnohat(n),PAW%AErefrv(n),PAW%rvx(n),PAW%trvx(n),&
&        PAW%den(n),PAW%tden(n),PAW%core(n),PAW%tcore(n),&
&        PAW%coretau(n),PAW%tcoretau(n),&
&        PAW%valetau(n),PAW%tvaletau(n),&
&        PAW%vtau(n),PAW%tvtau(n),&
&        PAW%nhatv(n),PAW%Ktvtau(n),PAW%Krveff(n),PAW%Kunscreen(n),stat=ok)
      IF (ok/=0) STOP 'Allocation error 1 in InitPAW'
      PAW%projshape=0.d0;PAW%hatden=0.d0;PAW%hatpot=0.d0
      PAW%hatshape=0.d0;PAW%vloc=0.d0;PAW%rveff=0.d0
      PAW%abinitvloc=0.d0;PAW%abinitnohat=0.d0
      PAW%AErefrv=0.d0;PAW%rvx=0.d0;PAW%trvx=0.d0
      PAW%den=0.d0;PAW%tden=0.d0;PAW%core=0.d0;PAW%tcore=0.d0
      PAW%XCORECORE=0.d0;PAW%nhatv=0.d0
      PAW%lambshielding=0.d0
      PAW%coretau=0.d0;PAW%tcoretau=0.d0
      PAW%valetau=0.d0;PAW%tvaletau=0.d0
      PAW%vtau=0.d0;PAW%tvtau=0.d0
      PAW%Ktvtau=0.d0;PAW%Krveff=0.d0;PAW%Kunscreen=0.d0
      ALLOCATE(PAW%phi(n,mxbase),PAW%tphi(n,mxbase),PAW%tp(n,mxbase),&
&        PAW%ophi(n,mxbase),PAW%otphi(n,mxbase),PAW%otp(n,mxbase),&
&        PAW%np(mxbase),PAW%l(mxbase),PAW%eig(mxbase),PAW%occ(mxbase),&
&        PAW%ck(mxbase),PAW%vrc(mxbase),PAW%Kop(n,mxbase),PAW%rng(mxbase),&
&        PAW%rcio(mxbase),PAW%nodes(mxbase),stat=ok)
      IF (ok/=0) STOP 'Allocation error 2 in InitPAW'
      PAW%phi=0.d0;PAW%tphi=0.d0;PAW%tp=0.d0
      PAW%ophi=0.d0;PAW%otphi=0.d0;PAW%otp=0.d0
      PAW%eig=0.d0;PAW%occ=0.d0;PAW%vrc=0.d0;PAW%ck=0.d0;PAW%Kop=0.d0
      PAW%rcio=0.d0;PAW%np=0;PAW%l=0
      if(diracrelativistic) then
        ALLOCATE(PAW%kappa(mxbase))
        PAW%kappa=0
      endif        
      PAW%rng=Grid%n
      ALLOCATE(PAW%oij(mxbase,mxbase),PAW%dij(mxbase,mxbase),&
&      PAW%wij(mxbase,mxbase), stat=ok)
      IF (ok/=0) STOP 'Allocation error 3 in InitPAW'
      PAW%oij=0.d0;PAW%dij=0.d0;PAW%wij=0.d0
      ALLOCATE(PAW%rVf(n),PAW%rtVf(n),PAW%Kij(mxbase,mxbase),&
&      PAW%Vfij(mxbase,mxbase),stat=ok)
      IF (ok/=0) STOP 'Allocation error 4 in InitPAW'
      PAW%rVf=0.d0;PAW%rtVf=0.d0;PAW%Kij=0.d0;PAW%Vfij=0.d0
      IF (Orbit%exctype=='HF') THEN
         ALLOCATE(PAW%lmbd(Orbit%norbit,mxbase),stat=ok)
         IF (ok/=0) STOP 'Allocation error 5 in InitPAW'
         PAW%lmbd=0.d0
      ELSE
         nullify(PAW%lmbd)
      ENDIF
      ALLOCATE(PAW%valencemap(Orbit%norbit),stat=ok)
      IF (ok/=0) STOP 'Allocation error 6 in InitPAW'
      ALLOCATE(PAW%OCCwfn,PAW%TOCCwfn,stat=ok)
      IF (ok/=0) STOP 'Allocation error 7 in InitPAW'
    END SUBROUTINE InitPAW

  Subroutine DestroyPAW(PAW)
    Type(PseudoInfo), INTENT(INOUT) :: PAW
    IF (ASSOCIATED(PAW%rcio)) DEALLOCATE(PAW%rcio)
    If (ASSOCIATED(PAW%vloc)) DEALLOCATE(PAW%vloc)
    If (ASSOCIATED(PAW%abinitvloc)) DEALLOCATE(PAW%abinitvloc)
    If (ASSOCIATED(PAW%abinitnohat)) DEALLOCATE(PAW%abinitnohat)
    If (ASSOCIATED(PAW%rveff)) DEALLOCATE(PAW%rveff)
    If (ASSOCIATED(PAW%AErefrv)) DEALLOCATE(PAW%AErefrv)
    If (ASSOCIATED(PAW%rvx)) DEALLOCATE(PAW%rvx)
    If (ASSOCIATED(PAW%trvx)) DEALLOCATE(PAW%trvx)
    If (ASSOCIATED(PAW%projshape)) DEALLOCATE(PAW%projshape)
    If (ASSOCIATED(PAW%hatshape)) DEALLOCATE(PAW%hatshape)
    If (ASSOCIATED(PAW%hatden)) DEALLOCATE(PAW%hatden)
    If (ASSOCIATED(PAW%hatpot)) DEALLOCATE(PAW%hatpot)
    If (ASSOCIATED(PAW%den)) DEALLOCATE(PAW%den)
    If (ASSOCIATED(PAW%tden)) DEALLOCATE(PAW%tden)
    If (ASSOCIATED(PAW%core)) DEALLOCATE(PAW%core)
    If (ASSOCIATED(PAW%tcore)) DEALLOCATE(PAW%tcore)
    If (ASSOCIATED(PAW%coretau)) DEALLOCATE(PAW%coretau)
    If (ASSOCIATED(PAW%tcoretau)) DEALLOCATE(PAW%tcoretau)
    If (ASSOCIATED(PAW%valetau)) DEALLOCATE(PAW%valetau)
    If (ASSOCIATED(PAW%tvaletau)) DEALLOCATE(PAW%tvaletau)
    If (ASSOCIATED(PAW%vtau)) DEALLOCATE(PAW%vtau)
    If (ASSOCIATED(PAW%tvtau)) DEALLOCATE(PAW%tvtau)
    If (ASSOCIATED(PAW%nhatv)) DEALLOCATE(PAW%nhatv)
    If (ASSOCIATED(PAW%np)) DEALLOCATE(PAW%np)
    If (ASSOCIATED(PAW%l)) DEALLOCATE(PAW%l)
    If (ASSOCIATED(PAW%nodes)) DEALLOCATE(PAW%nodes)
    If (ASSOCIATED(PAW%kappa)) DEALLOCATE(PAW%kappa)
    If (ASSOCIATED(PAW%rng)) DEALLOCATE(PAW%rng)
    If (ASSOCIATED(PAW%label)) DEALLOCATE(PAW%label)
    If (ASSOCIATED(PAW%phi)) DEALLOCATE(PAW%phi)
    If (ASSOCIATED(PAW%tphi)) DEALLOCATE(PAW%tphi)
    If (ASSOCIATED(PAW%tp)) DEALLOCATE(PAW%tp)
    If (ASSOCIATED(PAW%ophi)) DEALLOCATE(PAW%ophi)
    If (ASSOCIATED(PAW%otphi)) DEALLOCATE(PAW%otphi)
    If (ASSOCIATED(PAW%otp)) DEALLOCATE(PAW%otp)
    If (ASSOCIATED(PAW%Kop)) DEALLOCATE(PAW%Kop)
    If (ASSOCIATED(PAW%eig)) DEALLOCATE(PAW%eig)
    If (ASSOCIATED(PAW%occ)) DEALLOCATE(PAW%occ)
    If (ASSOCIATED(PAW%ck)) DEALLOCATE(PAW%ck)
    If (ASSOCIATED(PAW%vrc)) DEALLOCATE(PAW%vrc)
    If (ASSOCIATED(PAW%oij)) DEALLOCATE(PAW%oij)
    If (ASSOCIATED(PAW%dij)) DEALLOCATE(PAW%dij)
    If (ASSOCIATED(PAW%wij)) DEALLOCATE(PAW%wij)
    If (ASSOCIATED(PAW%rVf)) DEALLOCATE(PAW%rVf)
    If (ASSOCIATED(PAW%rtVf)) DEALLOCATE(PAW%rtVf)
    If (ASSOCIATED(PAW%g)) DEALLOCATE(PAW%g)
    If (ASSOCIATED(PAW%Kij)) DEALLOCATE(PAW%Kij)
    If (ASSOCIATED(PAW%Vfij)) DEALLOCATE(PAW%Vfij)
    If (ASSOCIATED(PAW%mLij)) DEALLOCATE(PAW%mLij)
    If (ASSOCIATED(PAW%DR)) DEALLOCATE(PAW%DR)
    If (ASSOCIATED(PAW%DRVC)) DEALLOCATE(PAW%DRVC)
    If (ASSOCIATED(PAW%TXVC)) DEALLOCATE(PAW%TXVC)
    If (ASSOCIATED(PAW%valencemap)) DEALLOCATE(PAW%valencemap)
    If (ASSOCIATED(PAW%lmbd)) DEALLOCATE(PAW%lmbd)
    If (ASSOCIATED(PAW%DRC)) DEALLOCATE(PAW%DRC)
    If (ASSOCIATED(PAW%DRCC)) DEALLOCATE(PAW%DRCC)
    If (ASSOCIATED(PAW%DRCjkl)) DEALLOCATE(PAW%DRCjkl)
    If (ASSOCIATED(PAW%mLic)) DEALLOCATE(PAW%mLic)
    If (ASSOCIATED(PAW%mLcc)) DEALLOCATE(PAW%mLcc)
    If (ASSOCIATED(PAW%Dcj)) DEALLOCATE(PAW%Dcj)
    If (ASSOCIATED(PAW%OCCwfn)) then
      call DestroyOrbit(PAW%OCCwfn)
      DEALLOCATE(PAW%OCCwfn)
    end if
    If (ASSOCIATED(PAW%TOCCwfn)) then
      call DestroyOrbit(PAW%TOCCwfn)
      DEALLOCATE(PAW%TOCCwfn)
    end if
  End Subroutine DestroyPAW

End module pseudodata
