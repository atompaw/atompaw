MODULE pseudodata
  Use atomdata

  IMPLICIT NONE

  TYPE  Pseudoinfo
     CHARACTER(132) :: exctype
     INTEGER  :: lmax,irc,irc_shap,irc_vloc,irc_core,coretailpoints,mesh_size
     INTEGER  :: ivale, ivion
     CHARACTER(132) :: Vloc_description
     CHARACTER(132) :: Proj_description
     CHARACTER(132) :: Comp_description
     LOGICAL :: multi_rc
     REAL(8) :: rc,rc_shap,rc_vloc,rc_core,energyoflmax,gausslength
     REAL(8), POINTER :: vloc(:),abinitvloc(:),abinitnohat(:)
     REAL(8), POINTER :: rveff(:),AErefrv(:),rvx(:),trvx(:)
     REAL(8), POINTER :: projshape(:),hatshape(:),hatden(:),hatpot(:)
     REAL(8), POINTER :: den(:),tden(:),core(:),tcore(:)
     INTEGER :: nbase,ncoreshell
     INTEGER, POINTER :: np(:),l(:)
     INTEGER, POINTER :: rng(:)       ! rng particularly of continuum states
     CHARACTER(8), POINTER :: label(:)
     REAL(8), POINTER :: phi(:,:),tphi(:,:),tp(:,:) ! before orthog
     REAL(8), POINTER :: ophi(:,:),otphi(:,:),otp(:,:) ! after orthog
     REAL(8), POINTER :: Kop(:,:)    ! for storing K|phi>
     REAL(8), POINTER :: eig(:),occ(:),ck(:),vrc(:)
     REAL(8), POINTER :: oij(:,:),dij(:,:),wij(:,:)
     !** L=0 matrix elements for atomic SC calculations
     REAL(8), POINTER :: tvij(:,:),vhatij(:,:),kin(:,:)
     REAL(8), POINTER :: v0ij(:,:),vhijkl(:,:,:,:)
     !********** modified parameters for use with KS and HF
     REAL(8), POINTER :: rVf(:),rtVf(:),g(:,:)
     REAL(8), POINTER :: Kij(:,:),Vfij(:,:),mLij(:,:,:),DR(:,:,:,:,:)
     REAL(8), POINTER :: DRVC(:,:,:)
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

  Subroutine DestroyPAW(PAW)
    Type(PseudoInfo), INTENT(INOUT) :: PAW

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
    If (ASSOCIATED(PAW%np)) DEALLOCATE(PAW%np)
    If (ASSOCIATED(PAW%l)) DEALLOCATE(PAW%l)
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
    If (ASSOCIATED(PAW%tvij)) DEALLOCATE(PAW%tvij)
    If (ASSOCIATED(PAW%vhatij)) DEALLOCATE(PAW%vhatij)
    If (ASSOCIATED(PAW%kin)) DEALLOCATE(PAW%kin)
    If (ASSOCIATED(PAW%v0ij)) DEALLOCATE(PAW%v0ij)
    If (ASSOCIATED(PAW%vhijkl)) DEALLOCATE(PAW%vhijkl)
    If (ASSOCIATED(PAW%rVf)) DEALLOCATE(PAW%rVf)
    If (ASSOCIATED(PAW%rtVf)) DEALLOCATE(PAW%rtVf)
    If (ASSOCIATED(PAW%g)) DEALLOCATE(PAW%g)
    If (ASSOCIATED(PAW%Kij)) DEALLOCATE(PAW%Kij)
    If (ASSOCIATED(PAW%Vfij)) DEALLOCATE(PAW%Vfij)
    If (ASSOCIATED(PAW%mLij)) DEALLOCATE(PAW%mLij)
    If (ASSOCIATED(PAW%DR)) DEALLOCATE(PAW%DR)
    If (ASSOCIATED(PAW%DRVC)) DEALLOCATE(PAW%DRVC)
    If (ASSOCIATED(PAW%valencemap)) DEALLOCATE(PAW%valencemap)
    If (ASSOCIATED(PAW%lmbd)) DEALLOCATE(PAW%lmbd)
    If (ASSOCIATED(PAW%DRC)) DEALLOCATE(PAW%DRC)
    If (ASSOCIATED(PAW%DRCC)) DEALLOCATE(PAW%DRCC)
    If (ASSOCIATED(PAW%DRCjkl)) DEALLOCATE(PAW%DRCjkl)
    If (ASSOCIATED(PAW%mLic)) DEALLOCATE(PAW%mLic)
    If (ASSOCIATED(PAW%mLcc)) DEALLOCATE(PAW%mLcc)
    If (ASSOCIATED(PAW%Dcj)) DEALLOCATE(PAW%Dcj)
    If (ASSOCIATED(PAW%OCCwfn)) Call DestroyOrbit(PAW%OCCwfn)
    If (ASSOCIATED(PAW%TOCCwfn)) Call DestroyOrbit(PAW%TOCCwfn)
  
  End Subroutine DestroyPAW

End module pseudodata
