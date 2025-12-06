!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!      SetPAWOptions1, SetPAWOptions2, StoreTOCCWFN, Troullier, kerker,
!        nonncps, checkghosts, sethat, coretailselfenergy, setcoretail,
!        fixtcorewfn, selfhatpot, setbasis, makebasis_bloechl,
!        makebasis_custom, makebasis_modrrkj, readmatchradius,
!        makebasis_V_setvloc, formprojectors, bsolv, ftprod, fthatpot,
!        ftkin, ftvloc, unboundsep, boundsep, PStoAE, Set_PAW_MatrixElements,
!        logderiv, FindVlocfromVeff, SCFPAW, PAWIter_LDA, exploreparms,
!        EXPLORElogderiv, Report_PseudobasisRP, phase_unwrap
!        Check_overlap_of_projectors
!        smoothcore, smoothtau, smoothexpcore, setttau , calculate_tvtau
!        resettcore, VPSmatch, Smoothtvtau
!
! 5/2018 phase_unwrap contributed by Casey Brock from Vanderbilt U. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE pseudo

  USE io_tools
  USE GlobalMath
  USE atomdata
  USE aeatom
  USE excor
  USE exx_pseudo
  USE gridmod 
  USE hf_pseudo
  USE numerov_mod
  USE paw_sub
  USE pseudodata
  USE pseudo_sub
  USE radialDirac
  USE radialked
  USE radialsr
  USE input_dataset_mod

  IMPLICIT NONE

  Type(Pseudoinfo), TARGET :: PAW

 !  Parameters controlling PAW options
  INTEGER,PRIVATE,PARAMETER :: BLOECHL=1, VANDERBILT=2, CUSTOM=3, MODRRKJ=7, HFPROJ=8
  INTEGER,PRIVATE,PARAMETER :: BLOECHLPS=0, POLYNOM=1, POLYNOM2=2, RRKJ=3
  INTEGER,PRIVATE,PARAMETER :: VANDERBILTORTHO=0, GRAMSCHMIDTORTHO=1
  INTEGER,PRIVATE,PARAMETER :: SVDORTHO=2, HFORTHO=-13
  INTEGER,PRIVATE,PARAMETER :: MTROULLIER=1, ULTRASOFT=2, BESSEL=3, KERKER_E=4, KERKER_P=5
  INTEGER,PRIVATE,PARAMETER :: HARTREE_FOCK=6, SETVLOC=7, VPSMATCHNC=8, VPSMATCHNNC=9

  REAL(8),PRIVATE, PARAMETER :: coretailtol=1.d-12
  INTEGER, PRIVATE :: coretailpoints=-1,besselopt=-1
  INTEGER, PRIVATE :: Projectorindex=-1,PSindex=-1,Orthoindex=-1,Vlocalindex=-1
  REAL(8), PRIVATE :: gaussparam=-1

CONTAINS

 SUBROUTINE SetPAWOptions1(ifen,Grid)
   INTEGER, INTENT(IN) :: ifen
   Type(Gridinfo), INTENT(IN) :: Grid

  INTEGER :: n,i,j,l,irc
  REAL(8):: h,rc
  LOGICAL :: multi_rc

  n=Grid%n

! Maximum l for basis functions (from input dataset)
  PAW%lmax=input_dataset%lmax

! Cut-off radii (from input dataset)
  multi_rc=(input_dataset%rc/=input_dataset%rc_shap.or.&
&           input_dataset%rc/=input_dataset%rc_vloc.or.&
&           input_dataset%rc/=input_dataset%rc_core.or.&
&           input_dataset%rc_shap/=input_dataset%rc_vloc.or.&
&           input_dataset%rc_shap/=input_dataset%rc_core.or.&
&           input_dataset%rc_vloc/=input_dataset%rc_core)
  PAW%irc     =FindGridIndex(Grid,input_dataset%rc)
  PAW%irc_shap=FindGridIndex(Grid,input_dataset%rc_shap)
  PAW%irc_vloc=FindGridIndex(Grid,input_dataset%rc_vloc)
  PAW%irc_core=FindGridIndex(Grid,input_dataset%rc_core)

  PAW%rc     =Grid%r(PAW%irc)
  PAW%rc_shap=Grid%r(PAW%irc_shap)
  PAW%rc_vloc=Grid%r(PAW%irc_vloc)
  PAW%rc_core=Grid%r(PAW%irc_core)
  irc=PAW%irc;rc=PAW%rc

  if (irc>n-Grid%ishift) stop 'error -- rc is too big !'
  WRITE(STD_OUT,*) ' adjusted rc ',rc, Grid%r(irc)
  WRITE(STD_OUT,*) ' irc,rc = ',irc,rc
  if (multi_rc) then
   WRITE(STD_OUT,*) ' adjusted rc_shape ',PAW%rc_shap
   WRITE(STD_OUT,*) ' adjusted rc_vloc  ',PAW%rc_vloc
   WRITE(STD_OUT,*) ' adjusted rc_core  ',PAW%rc_core
  endif
  
  if(input_dataset%coreshapemod) then
     PAW%tcoreshapeexp=.true.
     WRITE(STD_OUT,*) ' Using modified core shape functions'
  endif   
             
  WRITE(ifen,*) ' paw parameters: '
  WRITE(ifen,*) '      lmax = ',PAW%lmax
  WRITE(ifen,*) '        rc = ',PAW%rc
  WRITE(ifen,*) '       irc = ',PAW%irc
  if (multi_rc) then
   WRITE(ifen,*) '  rc_shape = ',PAW%rc_shap
   WRITE(ifen,*) '   rc_vloc = ',PAW%rc_vloc
   WRITE(ifen,*) '   rc_core = ',PAW%rc_core
  endif
 End Subroutine SetPAWOptions1

 SUBROUTINE SetPAWOptions2(ifen,Grid,Orbit,Pot,success)
   INTEGER, INTENT(IN) :: ifen
   Type(Gridinfo), INTENT(IN) :: Grid
   Type(OrbitInfo), INTENT(IN) :: Orbit
   Type(PotentialInfo), INTENT(IN) :: Pot
   LOGICAL , INTENT(OUT) :: success

  INTEGER :: i,pdeg,l
  REAL(8) :: qcut,x,y,e

  success=.true.

  !Pseudization and orthogonalization parameters (from input dataset)
   pdeg=input_dataset%pseudo_polynom2_pdeg
   qcut=input_dataset%pseudo_polynom2_qcut
   IF (input_dataset%projector_type==PROJECTOR_TYPE_BLOECHL) Projectorindex=BLOECHL  
   IF (input_dataset%projector_type==PROJECTOR_TYPE_VANDERBILT) Projectorindex=VANDERBILT  
   IF (input_dataset%projector_type==PROJECTOR_TYPE_MODRRKJ) Projectorindex=MODRRKJ  
   IF (input_dataset%projector_type==PROJECTOR_TYPE_CUSTOM) Projectorindex=CUSTOM  
   IF (input_dataset%projector_type==PROJECTOR_TYPE_HF) Projectorindex=HFPROJ
   IF (input_dataset%pseudo_type==PSEUDO_TYPE_BLOECHL) PSindex=BLOECHLPS
   IF (input_dataset%pseudo_type==PSEUDO_TYPE_POLYNOM) PSindex=POLYNOM
   IF (input_dataset%pseudo_type==PSEUDO_TYPE_POLYNOM2) PSindex=POLYNOM2
   IF (input_dataset%pseudo_type==PSEUDO_TYPE_RRKJ) PSindex=RRKJ
   IF (input_dataset%pseudo_type==PSEUDO_TYPE_HF) PSindex=HARTREE_FOCK
   IF (input_dataset%ortho_type==ORTHO_TYPE_GRAMSCHMIDT) Orthoindex=GRAMSCHMIDTORTHO
   IF (input_dataset%ortho_type==ORTHO_TYPE_VANDERBILT) Orthoindex=VANDERBILTORTHO
   IF (input_dataset%ortho_type==ORTHO_TYPE_SVD) Orthoindex=SVDORTHO
   IF (input_dataset%ortho_type==ORTHO_TYPE_HF) Orthoindex=HFORTHO

   write(PAW%Proj_description,'("Projector type:")')

   if (PSindex==BLOECHLPS) then
    write(PAW%Proj_description,'(a," Bloechl")') trim(PAW%Proj_description)
   else if (Projectorindex==MODRRKJ) then
    write(PAW%Proj_description,'(a," modified RKKJ projectors")') &
&            trim(PAW%Proj_description)
   else if (PSindex==POLYNOM) then
    write(PAW%Proj_description,'(a," polynomial pseudization")') &
&            trim(PAW%Proj_description)
   else if (PSindex==POLYNOM2) then
    write(PAW%Proj_description,'(a," improved polynomial pseudization")') &
&            trim(PAW%Proj_description)
   else if (PSindex==RRKJ) then
    write(PAW%Proj_description,'(a," RRKJ pseudization")') &
&         trim(PAW%Proj_description)
   else if (PSindex==HARTREE_FOCK) then
    write(PAW%Proj_description,'(" HF projectors using Vanderbilt-like scheme")')
   endif

   if (Orthoindex==VANDERBILTORTHO) then
    PAW%orthogonalization_scheme='vanderbilt'
    write(PAW%Proj_description,'(a," + Vanderbilt ortho.")') &
&         trim(PAW%Proj_description)
   else if (Orthoindex==GRAMSCHMIDTORTHO) then
    PAW%orthogonalization_scheme='gramschmidt'
    write(PAW%Proj_description,'(a," + Gram-Schmidt ortho.")') &
&         trim(PAW%Proj_description)
   else if (Orthoindex==SVDORTHO) then
    PAW%orthogonalization_scheme='svd'
    write(PAW%Proj_description,'(a," + SVD ortho.")') &
&         trim(PAW%Proj_description)
   else if (Orthoindex==HFORTHO) then
    PAW%orthogonalization_scheme='hartree-fock'
    write(PAW%Proj_description,'(a," + HF ortho.")') &
&         trim(PAW%Proj_description)
   end if

   write(std_out,*) PAW%Proj_description

  !Shape function parameters (from input dataset)
   gaussianshapefunction=.false.;besselshapefunction=.false.
   if (input_dataset%shapefunc_type==SHAPEFUNC_TYPE_GAUSSIAN) then
     gaussianshapefunction=.true.
     gaussparam=input_dataset%shapefunc_gaussian_param
     CALL sethat(Grid,PAW,gaussparam=gaussparam)    ! Gaussian shape function
     write(PAW%Comp_description,&
&      '("Gaussian compensation charge shape with gausstol = ",1p,1e12.4)') gaussparam
   else if (input_dataset%shapefunc_type==SHAPEFUNC_TYPE_BESSEL) then
     besselshapefunction=.true.
     CALL sethat(Grid,PAW,besselopt=i)               ! Bessel shape function
     if (PAW%irc_shap/=PAW%irc) then
      write(PAW%Comp_description,&
&      '("Bessel compensation charge shape zeroed at ",1p,1e12.4)') PAW%rc_shap
     else
      write(PAW%Comp_description,&
&      '("Bessel compensation charge shape zeroed at rc")')
     endif
   else
     CALL sethat(Grid,PAW)                          ! sinc^2 shape function
     if (PAW%irc_shap/=PAW%irc) then
      write(PAW%Comp_description,&
&      '("Sinc^2 compensation charge shape zeroed at ",1p,1e12.4)') PAW%rc_shap
     else
      write(PAW%Comp_description,&
&      '("Sinc^2 compensation charge shape zeroed at rc")')
     endif
   endif

   !Smooth core shape
   shapetcore=input_dataset%shapetcore
   if (shapetcore) then
     write(std_out,*) 'Resetting tcore with shapetcore'
     call resettcore(Grid,PAW)
   endif

   !Core tolerance for HF (from input dataset)
   PAW%coretol=MAX(input_dataset%hf_coretol,0.d0)
   if (Projectorindex==HFPROJ) write(std_out,*) 'Resetting coretol to ', PAW%coretol

   !Vlocal parameters (from input dataset)
   pdeg=input_dataset%pseudo_polynom2_pdeg
   qcut=input_dataset%pseudo_polynom2_qcut
   IF (input_dataset%vloc_type==VLOC_TYPE_MTROULLIER)  Vlocalindex=MTROULLIER
   IF (input_dataset%vloc_type==VLOC_TYPE_ULTRASOFT)   Vlocalindex=ULTRASOFT
   IF (input_dataset%vloc_type==VLOC_TYPE_BESSEL)      Vlocalindex=BESSEL
   IF (input_dataset%vloc_type==VLOC_TYPE_VPSMATCHNC)  Vlocalindex=VPSMATCHNC
   IF (input_dataset%vloc_type==VLOC_TYPE_VPSMATCHNNC) Vlocalindex=VPSMATCHNNC
   IF (input_dataset%vloc_type==VLOC_TYPE_SETVLOC)     Vlocalindex=SETVLOC
   IF (input_dataset%vloc_type==VLOC_TYPE_KERKER_EXPF) Vlocalindex=KERKER_E
   IF (input_dataset%vloc_type==VLOC_TYPE_KERKER_POLY) Vlocalindex=KERKER_P

   !Store the description of Vloc scheme in a string
   if (Vlocalindex==MTROULLIER) then
     l=input_dataset%vloc_l ; e=input_dataset%vloc_ene
     WRITE(PAW%Vloc_description,&
&      '("Vloc: Norm-conserving Troullier-Martins with l= ",i1,";e= ",1p,1e12.4)')l,e
   endif
   if (Vlocalindex==ULTRASOFT) then
     l=input_dataset%vloc_l ; e=input_dataset%vloc_ene
     WRITE(PAW%Vloc_description,&
&      '("Vloc: Non norm-conserving form with l= ",i1,";e= ",1p,1e12.4)')l,e
   endif
   if (Vlocalindex==VPSMATCHNC) then
     l=input_dataset%vloc_l ; e=input_dataset%vloc_ene
     WRITE(PAW%Vloc_description,&
&      '("Vloc: VPS match (norm-conservation) with l= ",i1,";e= ",1p,1e12.4)')l,e
   endif
   if (Vlocalindex==VPSMATCHNNC) then
     l=input_dataset%vloc_l ; e=input_dataset%vloc_ene
     WRITE(PAW%Vloc_description,&
&      '("Vloc: VPS match (without norm-conservation) with l= ",i1,";e= ",1p,1e12.4)')l,e
   endif
   if (Vlocalindex==BESSEL) then
     WRITE(PAW%Vloc_description,&
&        '("Vloc: truncated form - Vps(r)=A.sin(qr)/r for r<rc")')
   endif
   if (Vlocalindex==SETVLOC) then
     PAW%VlocCoef=input_dataset%vloc_setvloc_coef
     i=FindGridIndex(Grid,input_dataset%vloc_setvloc_rad)
     PAW%VlocRad=Grid%r(i)
     WRITE(PAW%Vloc_description, &
&      '("Vloc == VlocCoef*shapfunc , VlocCoef,Rad  = ",1p,2e15.7)') PAW%VlocCoef,PAW%VlocRad
     PAW%vloc= 0.d0;y=PAW%VlocRad
     PAW%vloc(1)=PAW%VlocCoef
     do i=2,Grid%n
       if (Grid%r(i)<PAW%VlocRad) then
         PAW%vloc(i)=PAW%VlocCoef*(SIN(pi*Grid%r(i)/y)/(pi*Grid%r(i)/y))**2
       endif
     enddo
   endif
   if (Vlocalindex==KERKER_E.or.Vlocalindex==KERKER_P) then
     l=input_dataset%vloc_l ; e=input_dataset%vloc_ene
     if (Vlocalindex==KERKER_E) PAW%Vloc_description="Norm-conserving Exp Vloc"
     if (Vlocalindex==KERKER_P) PAW%Vloc_description="Norm-conserving Poly Vloc"
     WRITE(PAW%Vloc_description,'(a,"; l = ",i1,"; powers = ",4i3,"; e = ",1p,e12.3)')&
&       trim(PAW%Vloc_description),l,input_dataset%vloc_kerker_power(1:4),e
   ENDIF

   WRITE(STD_OUT,*) PAW%Vloc_description

   !Call the routine computing Vloc - Not mGGA
   IF (.NOT.needvtau) THEN
	 IF (Vlocalindex==MTROULLIER.and.Projectorindex/=HFPROJ) THEN
       CALL troullier(Grid,Pot,PAW,l,e)
         ENDIF
	 IF (Vlocalindex==MTROULLIER.and.Projectorindex==HFPROJ) THEN
	   CALL make_hf_basis_only(Grid,Pot,PAW)
	   CALL troullier_HF(Grid,Pot,PAW,l,e)
	 ENDIF
	 IF (Vlocalindex==ULTRASOFT) CALL nonncps(Grid,Pot,PAW,l,e)
	 IF (Vlocalindex==BESSEL) CALL besselps(Grid,Pot,PAW)
	 IF (Vlocalindex==VPSMATCHNC) CALL VPSmatch(Grid,Pot,PAW,l,e,.true.)
	 IF (Vlocalindex==VPSMATCHNNC) CALL VPSmatch(Grid,Pot,PAW,l,e,.false.)
	 IF (Vlocalindex==KERKER_E.or.Vlocalindex==KERKER_P) CALL kerker(Grid,Pot,PAW)
	 IF (Projectorindex==BLOECHL) THEN
	   CALL makebasis_bloechl(Grid,Pot,0)
	 ELSE IF (Projectorindex==CUSTOM.AND.PSindex==BLOECHLPS) THEN
	   CALL makebasis_bloechl(Grid,Pot,1)
	 ELSE IF (Projectorindex==VANDERBILT.OR.Projectorindex==CUSTOM) THEN
	   !! Note: if SETVLOC only HF or VANDERBILT schemes work
	   if (Vlocalindex==SETVLOC) then
		  Call makebasis_V_setvloc(Grid,Pot,PAW)
	   else
		  CALL makebasis_custom(Grid,Pot,PSindex,Orthoindex,pdeg,qcut)
	   endif
	 ELSE IF (Projectorindex==HARTREE_FOCK) THEN
	  CALL make_hf_tp_only(Grid,Pot,PAW)
	 ELSE IF (Projectorindex==MODRRKJ) THEN
	  CALL makebasis_modrrkj(Grid,Pot,Orthoindex,success)
	 ENDIF
   ENDIF

   !Call the routine computing Vloc - mGGA case
   IF (needvtau) THEN
     !All compatibility checks in input_dataset_read routine
     WRITE(STD_OUT,*) 'Sequence of dataset construction modified for MGGA'      
     WRITE(STD_OUT,*) ' Not all possibilites tested carefully yet.... '
     WRITE(STD_OUT,*) ' Some possibilites not yet programmed.... '
     !Build basis
     IF (Projectorindex==VANDERBILT.OR.Projectorindex==CUSTOM) THEN
       CALL makebasis_custom_vtau1(Grid,Pot,PSindex,Orthoindex,pdeg,qcut)
     ELSE IF (Projectorindex==MODRRKJ) THEN
      CALL makebasis_modrrkj_vtau1(Grid,Pot,Orthoindex,success)
     ENDIF
     !Calculate PAW%vtau and PAW%tvtau     
     CALL calculate_tvtau(Grid)
     CALL Smoothtvtau(Grid,PAW,0.001d0)
     !Set pseudoptentials     
     IF (Vlocalindex==MTROULLIER.and.(TRIM(Orbit%exctype)/='HF')) then
       !!!WRITE(STD_OUT,*) 'TROULLIER PS not available for MGGA '
       !!!WRITE(STD_OUT,*) ' calling VPSmatch with norm conservation instead '     
       !!!CALL VPSmatch(Grid,Pot,PAW,l,e,.true.)
       CALL troullier(Grid,Pot,PAW,l,e)
     ENDIF         
     IF (Vlocalindex==VPSMATCHNNC) CALL VPSmatch(Grid,Pot,PAW,l,e,.false.)
     IF (Vlocalindex==VPSMATCHNC) CALL VPSmatch(Grid,Pot,PAW,l,e,.true.)
     IF (Vlocalindex==ULTRASOFT) CALL nonncps(Grid,Pot,PAW,l,e)
     IF (Vlocalindex==BESSEL) CALL besselps(Grid,Pot,PAW)
     !Calculate projectors
     CALL makeprojectors_vtau(Grid)
   ENDIF

   !Output in summary file
   IF (needvtau) THEN
     WRITE(ifen,*) 'Sequence of dataset construction steps modified for mGGA'      
     WRITE(ifen,*) 'Only projectors from Vanderbilt scheme available'
   ENDIF
   WRITE(ifen,*) TRIM(PAW%Vloc_description)
   WRITE(ifen,*) TRIM(PAW%Proj_description)
   WRITE(ifen,*) TRIM(PAW%Comp_description)
   if (shapetcore) WRITE(ifen,*) 'tcore reset with hatshape form'

   CALL StoreTOCCWFN(PAW)

 END SUBROUTINE SetPAWOptions2

     SUBROUTINE StoreTOCCWFN(PAW)
       TYPE(PseudoInfo), INTENT(INOUT) :: PAW

       INTEGER :: io,ib

       do io=1,PAW%TOCCWFN%norbit
          if (PAW%valencemap(io)>0) then
              ib=PAW%valencemap(io)
              PAW%TOCCWFN%wfn(:,io)=PAW%tphi(:,ib)
          endif
       enddo
     END SUBROUTINE StoreTOCCWFN
!!!!##  !***************************************************************
!!!!##  ! SUBROUTINE troullier(lmax,Grid,Pot)
!!!!##  !  Creates  screened norm-conserving pseudopotential following
!!!!##  !    approach of N. Troullier and J. L. Martins, PRB 43, 1993 (1991)
!!!!##  !    Uses p(r)=a0+f(r); f(r)=SUMm(Coef(m)*r^(2*m), where
!!!!##  !          m=1,2..6
!!!!##  !    Psi(r) = r^(l+1)*exp(p(r))
!!!!##  !***************************************************************
!!!!##  SUBROUTINE Troullier(Grid,Pot,PAW,l,e)
!!!!##    TYPE(Gridinfo), INTENT(IN) :: Grid
!!!!##    TYPE(Potentialinfo), INTENT(IN) :: Pot
!!!!##    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
!!!!##    INTEGER,INTENT(IN) :: l
!!!!##    REAL(8),INTENT(IN) :: e
!!!!##
!!!!##    REAL(8), ALLOCATABLE :: VNC(:)
!!!!##    REAL(8) :: A0,A,B,B0,C,C0,D,F,S
!!!!##    REAL(8) :: Coef(6),Coef0,Coef0old
!!!!##    REAL(8) :: h,rc,delta,x,pp,dpp,ddpp,dddpp,ddddpp
!!!!##    REAL(8) :: gam,bet
!!!!##    INTEGER :: i,j,k,n,iter,nr,nodes,irc,ok,m,wavetype
!!!!##    INTEGER, PARAMETER :: niter=5000
!!!!##    REAL(8), PARAMETER :: small=1.0d-9
!!!!##    REAL(8), ALLOCATABLE ::  wfn(:),p(:),dum(:),aux(:),extra1(:),extra2(:)
!!!!##    REAL(8), POINTER :: r(:),rv(:)
!!!!##    CHARACTER(132) :: line
!!!!##
!!!!##    n=Grid%n
!!!!##    h=Grid%h
!!!!##    r=>Grid%r
!!!!##    rv=>Pot%rv
!!!!##    nr=min(PAW%irc_vloc+10,n)
!!!!##    irc=PAW%irc_vloc
!!!!##    rc=PAW%rc_vloc
!!!!##
!!!!##    ALLOCATE(VNC(n),wfn(nr),p(nr),dum(nr),aux(nr),extra1(nr),extra2(nr),stat=ok)
!!!!##    IF (ok /=0) THEN
!!!!##       WRITE(STD_OUT,*) 'Error in troullier  -- in allocating wfn,p', nr,ok
!!!!##       STOP
!!!!##    ENDIF
!!!!##
!!!!##    !write(std_out,*) ' Troullier ', n,nr,irc
!!!!##    !call flush_unit(std_out)
!!!!##    if (scalarrelativistic) then
!!!!##       CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
!!!!##    else if (needvtau) then
!!!!##       CALL unboundked(Grid,Pot,nr,l,e,wfn,nodes)
!!!!##    else
!!!!##       CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,e,wfn,nodes)
!!!!##    endif
!!!!##
!!!!##    IF (wfn(irc)<0) wfn=-wfn
!!!!##    dum(1:irc)=(wfn(1:irc)**2)
!!!!##    S=integrator(Grid,dum(1:irc),1,irc)
!!!!##    A0=LOG(wfn(irc)/(rc**(l+1)))
!!!!##    B0=(rc*Gfirstderiv(Grid,irc,wfn)/wfn(irc)-(l+1))
!!!!##    C0=rc*(rv(irc)-rc*e)-B0*(B0+2*l+2)
!!!!##    D=-rc*(rv(irc)-rc*Gfirstderiv(Grid,irc,rv))-2*B0*C0-2*(l+1)*(C0-B0)
!!!!##    F=rc*(2*rv(irc)-rc*(2*Gfirstderiv(Grid,irc,rv) &
!!!!##&        -rc*Gsecondderiv(Grid,irc,rv)))+&
!!!!##&        4*(l+1)*(C0-B0)-2*(l+1)*D-2*C0**2-2*B0*D
!!!!##
!!!!##    WRITE(STD_OUT,*) 'In troullier -- matching parameters',S,A0,B0,C0,D,F
!!!!##
!!!!##    delta=1.d10
!!!!##    iter=0
!!!!##    Coef0=0
!!!!##
!!!!##    DO WHILE(delta>small.AND.iter<=niter)
!!!!##       iter=iter+1
!!!!##       A=A0-Coef0
!!!!##       B=B0
!!!!##       C=C0
!!!!##       CALL EvaluateTp(l,A,B,C,D,F,coef)
!!!!##
!!!!##       dum=0
!!!!##       DO  i=1,irc
!!!!##          x=(r(i)/rc)**2
!!!!##          p(i)=x*(Coef(1)+x*(Coef(2)+x*(Coef(3)+&
!!!!##&              x*(Coef(4)+x*(Coef(5)+x*Coef(6))))))
!!!!##          dum(i)=((r(i)**(l+1))*EXP(p(i)))**2
!!!!##       ENDDO
!!!!##       Coef0old=Coef0
!!!!##
!!!!##       x=integrator(Grid,dum(1:irc),1,irc)
!!!!##       Coef0=(LOG(S/x))/2
!!!!##
!!!!##       delta=ABS(Coef0-Coef0old)
!!!!##       !WRITE(STD_OUT,'(" VNC: iter Coef0 delta",i5,1p,2e15.7)') iter,Coef0,delta
!!!!##    ENDDO
!!!!##
!!!!##    WRITE(STD_OUT,*) '  VNC converged in ', iter,'  iterations'
!!!!##    WRITE(STD_OUT,*) '  Coefficients  -- ', Coef0,Coef(1:6)
!!!!##    !
!!!!##    ! Now  calculate VNC
!!!!##    extra1=0.d0;extra2=0.d0
!!!!##    if (needvtau) then
!!!!##       aux=0.d0
!!!!##       call derivative(Grid,PAW%tvtau,aux,1,nr)
!!!!##       WRITE(STD_OUT,*) 'In subroutine Troullier -- aux(1) = ',aux(1)
!!!!##       WRITE(STD_OUT,*) ' Resetting aux(1) to 0 '
!!!!##       aux(1)=0.d0 
!!!!##    endif     
!!!!##    OPEN(88,file='NC',form='formatted')
!!!!##    write(88,*) '# rc = ',r(irc)
!!!!##    !
!!!!##    VNC=0.d0;extra1=0.d0;extra2=0.d0
!!!!##    DO  i=2,nr
!!!!##       x=(r(i)/rc)**2
!!!!##       p(i)=Coef0+x*(Coef(1)+x*(Coef(2)+&
!!!!##&           x*(Coef(3)+x*(Coef(4)+x*(Coef(5)+x*Coef(6))))))
!!!!##       dpp=2*r(i)/(rc**2)*(Coef(1)+x*(2*Coef(2)+x*(3*Coef(3)+&
!!!!##&           x*(4*Coef(4)+x*(5*Coef(5)+x*6*Coef(6))))))
!!!!##       ddpp=(1/(rc**2))*(2*Coef(1)+x*(12*Coef(2)+x*(30*Coef(3)+&
!!!!##&           x*(56*Coef(4)+x*(90*Coef(5)+x*132*Coef(6))))))
!!!!##       dddpp=(r(i)/rc**4)*(24*Coef(2)+x*(120*Coef(3)+x*(336*Coef(4)+&
!!!!##&           x*(720*Coef(5)+x*1320*Coef(6)))))
!!!!##       ddddpp=(1/(rc**4)*(24*Coef(2)+x*(360*Coef(3)+x*(1680*Coef(4)+&
!!!!##&           x*(5040*Coef(5)+x*11880*Coef(6))))))
!!!!##       IF (i==irc) THEN
!!!!##          WRITE(STD_OUT,*) 'check  dp ', dpp,  B0/rc
!!!!##          WRITE(STD_OUT,*) 'check ddp ', ddpp, C0/rc**2
!!!!##          WRITE(STD_OUT,*) 'check dddp', dddpp, D/rc**3
!!!!##          WRITE(STD_OUT,*) 'check ddddp', ddddpp, F/rc**4
!!!!##       ENDIF
!!!!##       extra1(i)=e+ddpp+dpp*(dpp+2*(l+1)/r(i))
!!!!##       if (needvtau) then
!!!!##               extra2(i)=PAW%tvtau(i)*extra1(i)+aux(i)*(dpp+l/r(i))
!!!!##!               VNC(i)=e+(1.d0+PAW%tvtau(i))*(ddpp+dpp*(dpp+2*(l+1)/r(i))) &
!!!!##!&                +aux(i)*(dpp+l/r(i))                       
!!!!##                VNC(i)=extra1(i)+extra2(i)
!!!!##       else        
!!!!##!               VNC(i)=e+ddpp+dpp*(dpp+2*(l+1)/r(i))
!!!!##                VNC(i)=extra1(i)
!!!!##       endif
!!!!##               dum(i)=(r(i)**(l+1))*EXP(p(i))
!!!!##       WRITE(88,'(1p,7e15.7)') r(i),wfn(i),dum(i),VNC(i)*r(i),rv(i),extra1(i)*r(i),extra2(i)*r(i)
!!!!##    ENDDO
!!!!##    CLOSE(88)
!!!!##    x=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
!!!!##    WRITE(STD_OUT,*) 'check norm ',x,S
!!!!##
!!!!##    VNC(irc:n)=rv(irc:n)/r(irc:n)
!!!!##    PAW%rveff(1:n)=VNC(1:n)*r(1:n)
!!!!##
!!!!##    DEALLOCATE(VNC,wfn,p,dum,aux)
!!!!##  END SUBROUTINE troullier

  !***************************************************************
  ! SUBROUTINE troullier(lmax,Grid,Pot)
  !  Creates  screened norm-conserving pseudopotential following
  !    approach of N. Troullier and J. L. Martins, PRB 43, 1993 (1991)
  !    Uses p(r)=a0+f(r); f(r)=SUMm(Coef(m)*r^(2*m), where
  !          m=1,2..6
  !    Psi(r) = r^(l+1)*exp(p(r))
  !    modified 10/12/2023 by NAWH to accomodate possible metaGGA
  !***************************************************************
  SUBROUTINE Troullier(Grid,Pot,PAW,l,e)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
    INTEGER,INTENT(IN) :: l
    REAL(8),INTENT(IN) :: e

    REAL(8), ALLOCATABLE :: VNC(:)
    REAL(8) :: A0,A,B,B0,C,C0,D,F,S
    REAL(8) :: Coef(6),Coef0,Coef0old
    REAL(8) :: h,rc,delta,x,pp,dpp,ddpp,dddpp,ddddpp
    REAL(8) :: gam,bet
    INTEGER :: i,j,k,n,iter,nr,nodes,irc,ok,m,wavetype
    INTEGER, PARAMETER :: niter=5000
    REAL(8), PARAMETER :: small=1.0d-9
    REAL(8), ALLOCATABLE ::  wfn(:),p(:),dum(:),aux(:),Kaux(:),extra1(:),extra2(:),extra3(:)
    REAL(8), ALLOCATABLE :: pote(:),dpote(:),ddpote(:),logvtau(:),dlogvtau(:),ddlogvtau(:),dddlogvtau(:)
    REAL(8), POINTER :: r(:),rv(:)
    CHARACTER(132) :: line

    n=Grid%n
    h=Grid%h
    r=>Grid%r
    rv=>Pot%rv
    nr=min(PAW%irc_vloc+10,n)
    irc=PAW%irc_vloc
    rc=PAW%rc_vloc

    ALLOCATE(VNC(n),wfn(nr),p(nr),dum(nr),aux(nr),Kaux(nr),extra1(nr),extra2(nr),extra3(nr),stat=ok)
    ALLOCATE(pote(nr),dpote(nr),ddpote(nr),logvtau(nr),dlogvtau(nr),ddlogvtau(nr),dddlogvtau(nr),stat=ok)
    IF (ok /=0) THEN
       WRITE(STD_OUT,*) 'Error in troullier  -- in allocating wfn,p', nr,ok
       STOP
    ENDIF

    !write(std_out,*) ' Troullier ', n,nr,irc
    !call flush_unit(std_out)
    if (scalarrelativistic) then
       CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
    else if (needvtau) then
       CALL unboundked(Grid,Pot,nr,l,e,wfn,nodes)
    else
       CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,e,wfn,nodes)
    endif

    IF (wfn(irc)<0) wfn=-wfn
    dum(1:irc)=(wfn(1:irc)**2)
    S=integrator(Grid,dum(1:irc),1,irc)
    A0=LOG(wfn(irc)/(rc**(l+1)))
    B0=(rc*Gfirstderiv(Grid,irc,wfn)/wfn(irc)-(l+1))
    pote=0.d0;pote(2:nr)=rv(2:nr)/r(2:nr)-e
    if (needvtau) pote(1:nr)=pote(1:nr)/(1.d0+Pot%vtau(1:nr))
    logvtau=0.d0;dlogvtau=0.d0;ddlogvtau=0.d0;dddlogvtau=0.d0
    if (needvtau) then  !   find logvtau terms and their derivatives
       logvtau(1:nr)=LOG(1.d0+Pot%vtau(1:nr))     
       call derivative(Grid,logvtau,dlogvtau,1,nr)
       call derivative(Grid,dlogvtau,ddlogvtau,1,nr)
       call derivative(Grid,ddlogvtau,dddlogvtau,1,nr)
    endif
    C0=rc*rc*pote(irc)-B0*(B0+2*l+2)-rc*dlogvtau(irc)*(B0+l)
    D=rc*rc*rc*Gfirstderiv(Grid,irc,pote)-2*B0*C0-2*(l+1)*(C0-B0)
    if (needvtau) then
        D=D-rc*rc*ddlogvtau(irc)*(B0+l)-rc*dlogvtau(irc)*(C0-l)    
    endif
    F=rc*rc*rc*rc*Gsecondderiv(Grid,irc,pote)+&
&        4*(l+1)*(C0-B0)-2*(l+1)*D-2*C0**2-2*B0*D
    if (needvtau) then
        F=F-rc*rc*rc*dddlogvtau(irc)*(B0+l)-2*rc*rc*ddlogvtau(irc)*(C0-l) &
&           -rc*dlogvtau(irc)*(D+2*l)                
    endif

    WRITE(STD_OUT,*) 'In troullier -- matching parameters',S,A0,B0,C0,D,F

    delta=1.d10
    iter=0
    Coef0=0

    DO WHILE(delta>small.AND.iter<=niter)
       iter=iter+1
       A=A0-Coef0
       B=B0
       C=C0
       CALL EvaluateTp(l,A,B,C,D,F,coef)

       dum=0
       DO  i=1,irc
          x=(r(i)/rc)**2
          p(i)=x*(Coef(1)+x*(Coef(2)+x*(Coef(3)+&
&              x*(Coef(4)+x*(Coef(5)+x*Coef(6))))))
          dum(i)=((r(i)**(l+1))*EXP(p(i)))**2
       ENDDO
       Coef0old=Coef0

       x=integrator(Grid,dum(1:irc),1,irc)
       Coef0=(LOG(S/x))/2

       delta=ABS(Coef0-Coef0old)
       !WRITE(STD_OUT,'(" VNC: iter Coef0 delta",i5,1p,2e15.7)') iter,Coef0,delta
    ENDDO

    WRITE(STD_OUT,*) '  VNC converged in ', iter,'  iterations'
    WRITE(STD_OUT,*) '  Coefficients  -- ', Coef0,Coef(1:6)
    !
    ! Now  calculate VNC
    extra1=0.d0;extra2=0.d0;extra3=0.d0
    !##    if (needvtau) then
    !##       aux=0.d0
    !##       call derivative(Grid,PAW%tvtau,aux,1,nr)
    !##       WRITE(STD_OUT,*) 'In subroutine Troullier -- aux(1:10) = ',aux(1:10)
    !##       WRITE(STD_OUT,*) ' Resetting aux(1:10) to 0 '
    !##       aux(1:10)=0.d0 
    !##       Kaux=0.d0
    !##       call derivative(Grid,PAW%Ktvtau,Kaux,1,nr)
    !##       WRITE(STD_OUT,*) 'In subroutine Troullier -- Kaux(1:10) = ',Kaux(1:10)
    !##       WRITE(STD_OUT,*) ' Resetting Kaux(1:10) to 0 '
    !##       Kaux(1:10)=0.d0 
    !##    endif     
    OPEN(88,file='NC',form='formatted')
    write(88,*) '# rc = ',r(irc)
    write(88,*) 'r    wfn    dum    VNC*r    rv    extra1*r    extra2*r    extra3*r   '
    !
    VNC=0.d0;extra1=0.d0;extra2=0.d0;extra3=0.d0
    DO  i=2,nr
       x=(r(i)/rc)**2
       p(i)=Coef0+x*(Coef(1)+x*(Coef(2)+&
&           x*(Coef(3)+x*(Coef(4)+x*(Coef(5)+x*Coef(6))))))
       dpp=2*r(i)/(rc**2)*(Coef(1)+x*(2*Coef(2)+x*(3*Coef(3)+&
&           x*(4*Coef(4)+x*(5*Coef(5)+x*6*Coef(6))))))
       ddpp=(1/(rc**2))*(2*Coef(1)+x*(12*Coef(2)+x*(30*Coef(3)+&
&           x*(56*Coef(4)+x*(90*Coef(5)+x*132*Coef(6))))))
       dddpp=(r(i)/rc**4)*(24*Coef(2)+x*(120*Coef(3)+x*(336*Coef(4)+&
&           x*(720*Coef(5)+x*1320*Coef(6)))))
       ddddpp=(1/(rc**4)*(24*Coef(2)+x*(360*Coef(3)+x*(1680*Coef(4)+&
&           x*(5040*Coef(5)+x*11880*Coef(6))))))
       IF (i==irc) THEN
          WRITE(STD_OUT,*) 'check  dp ', dpp,  B0/rc  ; call flush_unit(std_out)
          WRITE(STD_OUT,*) 'check ddp ', ddpp, C0/rc**2; call flush_unit(std_out)
          WRITE(STD_OUT,*) 'check dddp', dddpp, D/rc**3; call flush_unit(std_out)
          WRITE(STD_OUT,*) 'check ddddp', ddddpp, F/rc**4; call flush_unit(std_out)
       ENDIF
       extra1(i)=ddpp+dpp*(dpp+2*(l+1)/r(i))
       if (needvtau) then
               extra2(i)=PAW%tvtau(i)*extra1(i)+PAW%dtvtaudr(i)*(dpp+l/r(i))
               extra3(i)=(1.d0+PAW%Ktvtau(i))*extra1(i)+PAW%Kdtvtaudr(i)*(dpp+l/r(i))  !  Kresse form
                VNC(i)=e+extra1(i)+extra2(i)   ! Bloechl form
       else        
                VNC(i)=e+extra1(i)
       endif
               dum(i)=(r(i)**(l+1))*EXP(p(i))
       WRITE(88,'(1p,9e15.7)') r(i),wfn(i),dum(i),VNC(i)*r(i),rv(i),extra1(i)*r(i),extra2(i)*r(i),extra3(i)*r(i)
    ENDDO
    CLOSE(88)
    x=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
    WRITE(STD_OUT,*) 'check norm ',x,S; call flush_unit(std_out)

    PAW%rveff=rv
    PAW%Krveff=rv
    PAW%rveff(1:irc-1)=VNC(1:irc-1)*r(1:irc-1)
    PAW%Krveff(1:irc-1)=extra3(1:irc-1)*r(1:irc-1)

    DEALLOCATE(VNC,wfn,p,dum,aux,Kaux,extra1,extra2,extra3)
    DEALLOCATE(pote,dpote,ddpote,logvtau,dlogvtau,ddlogvtau,dddlogvtau)
  END SUBROUTINE troullier

!++  !***************************************************************
!++  ! SUBROUTINE kerker(lmax,Grid,Pot)
!++  !  Creates  screened norm-conserving pseudopotential following
!++  !    approach of G. P. Kerker, J. Phys. C. 13,L189-L194 (1980)
!++  !    Uses p(r)=a0+f(r); f(r)=SUMi(Coef(i)*r^m(i)), where m(i)
!++  !          are input powers
!++  !    Psi(r) = r^(l+1)*exp(p(r)) if PStype = EXPF
!++  !    Psi(r) = r^(l+1)*(p(r))    if PStype = POLY
!++  !    Note that this subroutine may not be called in this
!++  !      version of ATOMPAW ...
!++  !***************************************************************
!++  SUBROUTINE kerker(Grid,Pot,PAW)
!++    TYPE(Gridinfo), INTENT(IN) :: Grid
!++    TYPE(Potentialinfo), INTENT(IN) :: Pot
!++    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
!++
!++    REAL(8), ALLOCATABLE :: VNC(:)
!++    REAL(8) :: A0,A,B,C,D,S,Coef(4),Coef0,Coef0old
!++    REAL(8) :: h,e,rc,delta,x,pp,dpp,ddpp,dddpp
!++    REAL(8) :: gam,bet
!++    INTEGER :: i,j,k,n,iter,nr,nodes,irc,l,ok,m(4),wavetype
!++    INTEGER, PARAMETER :: EXPF=1, POLY=2
!++    INTEGER, PARAMETER :: niter=5000
!++    REAL(8), PARAMETER :: small=1.0d-12
!++    CHARACTER(10) :: vtype
!++    REAL(8), ALLOCATABLE ::  wfn(:),p(:),dum(:),aux(:)
!++    REAL(8), POINTER :: r(:),rv(:)
!++
!++    !Read data from input dataset
!++    IF (input_dataset%vloc_type==VLOC_TYPE_KERKER_EXPF) wavetype=EXPF
!++    IF (input_dataset%vloc_type==VLOC_TYPE_KERKER_POLY) wavetype=POLY
!++    l=input_dataset%vloc_l ; e=input_dataset%vloc_ene
!++    m(1:4)=input_dataset%vloc_kerker_power(1:4)
!++
!++    n=Grid%n
!++    n=Grid%n
!++    h=Grid%h
!++    r=>Grid%r
!++    rv=>Pot%rv
!++    nr=min(PAW%irc_vloc+10,n)
!++    irc=PAW%irc_vloc
!++    rc=PAW%rc_vloc
!++    ALLOCATE(VNC(n),wfn(nr),p(nr),dum(nr),aux(nr),stat=ok)
!++    IF (ok /=0) THEN
!++       WRITE(STD_OUT,*) 'Error in kerker  -- in allocating wfn,p', nr,ok
!++       STOP
!++    ENDIF
!++
!++    if (scalarrelativistic) then
!++       CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
!++    else if (needvtau) then
!++       CALL unboundked(Grid,Pot,nr,l,e,wfn,nodes)
!++    else
!++       CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,e,wfn,nodes)
!++    endif
!++
!++
!++
!++    IF (wfn(irc)<0) wfn=-wfn
!++    dum(1:irc)=(wfn(1:irc)**2)
!++    S=integrator(Grid,dum(1:irc),1,irc)
!++    IF (wavetype==EXPF) THEN
!++       A0=LOG(wfn(irc)/(rc**(l+1)))
!++       B=(rc*Gfirstderiv(Grid,irc,wfn)/wfn(irc)-(l+1))
!++       C=rc*(rv(irc)-rc*e)-B*(B+2*l+2)
!++       D=-rc*(rv(irc)-rc*Gfirstderiv(Grid,irc,rv))-2*B*C-2*(l+1)*(C-B)
!++    ENDIF
!++
!++    IF (wavetype==POLY) THEN
!++       A0=(wfn(irc)/(rc**(l+1)))
!++       B=(rc*Gfirstderiv(Grid,irc,wfn))/(rc**(l+1))-(l+1)*A0
!++       C=rc*(rv(irc)-rc*e)*A0-2*(l+1)*B
!++       D=-rc*(rv(irc)-rc*Gfirstderiv(Grid,irc,rv))*A0+2*(l+1)*(B-C)+&
!++&           rc*(rv(irc)-rc*e)*B
!++    ENDIF
!++
!++
!++    WRITE(STD_OUT,*) 'In kerker -- matching parameters',S,A0,B,C,D
!++
!++    delta=1.d10
!++    iter=0
!++    Coef0=0
!++
!++    DO WHILE(delta>small.AND.iter<=niter)
!++       iter=iter+1
!++       A=A0-Coef0
!++       CALL EvaluateP(m,A,B,C,D,Coef)
!++
!++       dum=0
!++       DO  i=1,irc
!++          x=(r(i)/rc)
!++          p(i)=(x**m(1))*Coef(1)+(x**m(2))*Coef(2)+(x**m(3))*Coef(3)+(x**m(4))*Coef(4)
!++          IF (wavetype==EXPF)dum(i)=((r(i)**(l+1))*EXP(p(i)))**2
!++          IF (wavetype==POLY)dum(i)=(wfn(i))**2-((r(i)**(l+1))*(p(i)))**2
!++       ENDDO
!++       Coef0old=Coef0
!++       IF (wavetype==EXPF) THEN
!++          x=integrator(Grid,dum(1:irc),1,irc)
!++          Coef0=(LOG(S/x))/2
!++       ENDIF
!++       IF (wavetype==POLY) THEN
!++          gam=(2*l+3)*integrator(Grid,dum(1:irc),1,irc)/(rc**(2*l+3))
!++          bet=(2*l+3)*(Coef(1)/(2*l+3+m(1))+Coef(2)/(2*l+3+m(2))+&
!++&              Coef(3)/(2*l+3+m(3))+Coef(4)/(2*l+3+m(4)))
!++          !WRITE(STD_OUT,'("VNC: iter -- bet,gam = ",i5,1p,4e15.7)') iter,bet,gam
!++          x=bet**2+gam
!++          Coef0old=Coef0
!++          IF (x<0.d0) THEN
!++             WRITE(STD_OUT,*) 'Warning in Kerker subroutine x = ',x
!++               Coef0=Coef0+0.1*A0
!++            ELSE
!++               Coef0=SQRT(x)-bet
!++            ENDIF
!++         ENDIF
!++         delta=ABS(Coef0-Coef0old)
!++         !WRITE(STD_OUT,*) '  VNC: iter  Coef0  delta', iter,Coef0,delta
!++      ENDDO
!++
!++      WRITE(STD_OUT,*) '  VNC converged in ', iter,'  iterations'
!++      WRITE(STD_OUT,*) '  Coefficients  -- ', Coef0,Coef(1:4)
!++      !
!++      ! Now  calculate VNC
!++    if (needvtau) then
!++       aux=0.d0
!++       call derivative(Grid,PAW%tvtau,aux,1,nr)
!++    endif     
!++      OPEN(88,file='NC',form='formatted')
!++      !
!++      VNC=0
!++      DO  i=1,nr
!++         x=(r(i)/rc)
!++         p(i)=Coef0+(x**m(1))*Coef(1)+(x**m(2))*Coef(2)+&
!++&             (x**m(3))*Coef(3)+(x**m(4))*Coef(4)
!++         dpp=(m(1)*(x**(m(1)-1))*Coef(1)+m(2)*(x**(m(2)-1))*Coef(2)+&
!++&             m(3)*(x**(m(3)-1))*Coef(3)+m(4)*(x**(m(4)-1))*Coef(4))/rc
!++         ddpp=(m(1)*(m(1)-1)*(x**(m(1)-2))*Coef(1)+&
!++&             m(2)*(m(2)-1)*(x**(m(2)-2))*Coef(2)+&
!++&             m(3)*(m(3)-1)*(x**(m(3)-2))*Coef(3)+&
!++&             m(4)*(m(4)-1)*(x**(m(4)-2))*Coef(4))/(rc**2)
!++         dddpp=(m(1)*(m(1)-1)*(m(1)-2)*(x**(m(1)-3))*Coef(1)+&
!++&             m(2)*(m(2)-1)*(m(2)-2)*(x**(m(2)-3))*Coef(2)+&
!++&             m(3)*(m(3)-1)*(m(3)-2)*(x**(m(3)-3))*Coef(3)+&
!++&             m(4)*(m(4)-1)*(m(4)-2)*(x**(m(4)-3))*Coef(4))/(rc**3)
!++         IF (i==irc) THEN
!++            WRITE(STD_OUT,*) 'check  dp ', dpp,  B/rc
!++            WRITE(STD_OUT,*) 'check ddp ', ddpp, C/rc**2
!++            WRITE(STD_OUT,*) 'check dddp', dddpp,  D/rc**3
!++         ENDIF
!++         IF (wavetype==EXPF) THEN
!++           if (needvtau) then
!++               VNC(i)=e+(1.d0+PAW%tvtau(i))*(ddpp+dpp*(dpp+2*(l+1)/r(i))) &
!++&                +aux(i)*(dpp+l/r(i))                       
!++           else        
!++               VNC(i)=e+ddpp+dpp*(dpp+2*(l+1)/r(i))
!++           endif
!++            dum(i)=(r(i)**(l+1))*EXP(p(i))
!++         ENDIF
!++         IF (wavetype==POLY) THEN
!++           if (needvtau) then
!++               VNC(i)=e+(1.d0+PAW%tvtau(i))*(ddpp+2*(l+1)*dpp/r(i))/p(i) &
!++&                +aux(i)*(dpp/p(i)+l/r(i))                       
!++           else        
!++               VNC(i)=e+(ddpp+2*(l+1)*dpp/r(i))/p(i)
!++           endif
!++            dum(i)=(r(i)**(l+1))*(p(i))
!++         ENDIF
!++         WRITE(88,'(1p,5e15.7)') r(i),wfn(i),dum(i),VNC(i)*r(i),rv(i)
!++      ENDDO
!++      CLOSE(88)
!++      x=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
!++      WRITE(STD_OUT,*) 'check norm ',x,S
!++
!++      VNC(irc:n)=rv(irc:n)/r(irc:n)
!++      PAW%rveff(1:n)=VNC(1:n)*r(1:n)
!++
!++      DEALLOCATE(VNC,wfn,p,dum,aux)
!++    END SUBROUTINE kerker
!++ 

  !***************************************************************
  ! SUBROUTINE kerker(lmax,Grid,Pot)
  !  Creates  screened norm-conserving pseudopotential following
  !    approach of G. P. Kerker, J. Phys. C. 13,L189-L194 (1980)
  !    Uses p(r)=a0+f(r); f(r)=SUMi(Coef(i)*r^m(i)), where m(i)
  !          are input powers
  !    Psi(r) = r^(l+1)*exp(p(r)) if PStype = EXPF
  !    Psi(r) = r^(l+1)*(p(r))    if PStype = POLY
  !    Note that this subroutine may not be called in this
  !      version of ATOMPAW ...
  !  Modified 10/25/2023 by NAWH hopefully for compatibility
  !   with meta-GGA
  !***************************************************************
  SUBROUTINE kerker(Grid,Pot,PAW)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW

    REAL(8), ALLOCATABLE :: VNC(:)
    REAL(8) :: A0,A,B,C,D,S,Coef(4),Coef0,Coef0old
    REAL(8) :: h,e,rc,delta,x,pp,dpp,ddpp,dddpp
    REAL(8) :: gam,bet
    INTEGER :: i,j,k,n,iter,nr,nodes,irc,l,ok,m(4),wavetype
    INTEGER, PARAMETER :: EXPF=1, POLY=2
    INTEGER, PARAMETER :: niter=5000
    REAL(8), PARAMETER :: small=1.0d-12
    CHARACTER(10) :: vtype
    REAL(8), ALLOCATABLE ::  wfn(:),p(:),dum(:),aux(:)
    REAL(8), POINTER :: r(:),rv(:)
    REAL(8), ALLOCATABLE :: pote(:),dpote(:),ddpote(:),logvtau(:),dlogvtau(:),ddlogvtau(:),dddlogvtau(:)

    !Read data from input dataset
    IF (input_dataset%vloc_type==VLOC_TYPE_KERKER_EXPF) wavetype=EXPF
    IF (input_dataset%vloc_type==VLOC_TYPE_KERKER_POLY) wavetype=POLY
    l=input_dataset%vloc_l ; e=input_dataset%vloc_ene
    m(1:4)=input_dataset%vloc_kerker_power(1:4)

    n=Grid%n
    n=Grid%n
    h=Grid%h
    r=>Grid%r
    rv=>Pot%rv
    nr=min(PAW%irc_vloc+10,n)
    irc=PAW%irc_vloc
    rc=PAW%rc_vloc
    ALLOCATE(VNC(n),wfn(nr),p(nr),dum(nr),aux(nr),stat=ok)
    ALLOCATE(pote(nr),dpote(nr),ddpote(nr),logvtau(nr),dlogvtau(nr),ddlogvtau(nr),dddlogvtau(nr),stat=ok)
    IF (ok /=0) THEN
       WRITE(STD_OUT,*) 'Error in kerker  -- in allocating wfn,p', nr,ok
       STOP
    ENDIF

    if (scalarrelativistic) then
       CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
    else if (needvtau) then
       CALL unboundked(Grid,Pot,nr,l,e,wfn,nodes)
    else
       CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,e,wfn,nodes)
    endif



    IF (wfn(irc)<0) wfn=-wfn
    dum(1:irc)=(wfn(1:irc)**2)
    S=integrator(Grid,dum(1:irc),1,irc)
       pote=0.d0;pote(2:nr)=rv(2:nr)/r(2:nr)-e
       if (needvtau) pote(1:nr)=pote(1:nr)/(1.d0+Pot%vtau(1:nr))
       logvtau=0.d0;dlogvtau=0.d0;ddlogvtau=0.d0;dddlogvtau=0.d0
       if (needvtau) then  !   find logvtau terms and their derivatives
          logvtau(1:nr)=LOG(1.d0+Pot%vtau(1:nr))
          call derivative(Grid,logvtau,dlogvtau,1,nr)
          call derivative(Grid,dlogvtau,ddlogvtau,1,nr)
          call derivative(Grid,ddlogvtau,dddlogvtau,1,nr)
       endif
    IF (wavetype==EXPF) THEN
       A0=LOG(wfn(irc)/(rc**(l+1)))
       B=(rc*Gfirstderiv(Grid,irc,wfn)/wfn(irc)-(l+1))
       C=rc*rc*pote(irc)-B*(B+2*l+2)-rc*dlogvtau(irc)*(B+l)
       D=rc*rc*rc*Gfirstderiv(Grid,irc,pote)-2*B*C-2*(l+1)*(C-B)
       if (needvtau) then
           D=D-rc*rc*ddlogvtau(irc)*(B+l)-rc*dlogvtau(irc)*(C-l)
       endif
    ENDIF

    IF (wavetype==POLY) THEN
       A0=(wfn(irc)/(rc**(l+1)))
       B=(rc*Gfirstderiv(Grid,irc,wfn))/(rc**(l+1))-(l+1)*A0
       C=rc*rc*pote(irc)*A0-B*(2*l+2)-rc*dlogvtau(irc)*(B+l*A0)
       D=rc*rc*rc*Gfirstderiv(Grid,irc,pote)*A0+rc*rc*pote(irc)*B-2*(l+1)*(C-B)
       if (needvtau) then
           D=D-rc*rc*ddlogvtau(irc)*(B+l*A0)-rc*dlogvtau(irc)*(C+l*B-l*A0)
       endif
    ENDIF


    WRITE(STD_OUT,*) 'In kerker -- matching parameters',S,A0,B,C,D

    delta=1.d10
    iter=0
    Coef0=0

    DO WHILE(delta>small.AND.iter<=niter)
       iter=iter+1
       A=A0-Coef0
       CALL EvaluateP(m,A,B,C,D,Coef)

       dum=0
       DO  i=1,irc
          x=(r(i)/rc)
          p(i)=(x**m(1))*Coef(1)+(x**m(2))*Coef(2)+(x**m(3))*Coef(3)+(x**m(4))*Coef(4)
          IF (wavetype==EXPF)dum(i)=((r(i)**(l+1))*EXP(p(i)))**2
          IF (wavetype==POLY)dum(i)=(wfn(i))**2-((r(i)**(l+1))*(p(i)))**2
       ENDDO
       Coef0old=Coef0
    ENDDO


    WRITE(STD_OUT,*) 'In kerker -- matching parameters',S,A0,B,C,D

    delta=1.d10
    iter=0
    Coef0=0

    DO WHILE(delta>small.AND.iter<=niter)
       iter=iter+1
       A=A0-Coef0
       CALL EvaluateP(m,A,B,C,D,Coef)

       dum=0
       DO  i=1,irc
          x=(r(i)/rc)
          p(i)=(x**m(1))*Coef(1)+(x**m(2))*Coef(2)+(x**m(3))*Coef(3)+(x**m(4))*Coef(4)
          IF (wavetype==EXPF)dum(i)=((r(i)**(l+1))*EXP(p(i)))**2
          IF (wavetype==POLY)dum(i)=(wfn(i))**2-((r(i)**(l+1))*(p(i)))**2
       ENDDO
       Coef0old=Coef0
       IF (wavetype==EXPF) THEN
          x=integrator(Grid,dum(1:irc),1,irc)
          Coef0=(LOG(S/x))/2
       ENDIF
       IF (wavetype==POLY) THEN
          gam=(2*l+3)*integrator(Grid,dum(1:irc),1,irc)/(rc**(2*l+3))
          bet=(2*l+3)*(Coef(1)/(2*l+3+m(1))+Coef(2)/(2*l+3+m(2))+&
&              Coef(3)/(2*l+3+m(3))+Coef(4)/(2*l+3+m(4)))
          !WRITE(STD_OUT,'("VNC: iter -- bet,gam = ",i5,1p,4e15.7)') iter,bet,gam
          x=bet**2+gam
          Coef0old=Coef0
          IF (x<0.d0) THEN
             WRITE(STD_OUT,*) 'Warning in Kerker subroutine x = ',x
               Coef0=Coef0+0.1*A0
            ELSE
               Coef0=SQRT(x)-bet
          ENDIF
       ENDIF
         delta=ABS(Coef0-Coef0old)
         !WRITE(STD_OUT,*) '  VNC: iter  Coef0  delta', iter,Coef0,delta
      ENDDO

      WRITE(STD_OUT,*) '  VNC converged in ', iter,'  iterations'
      WRITE(STD_OUT,*) '  Coefficients  -- ', Coef0,Coef(1:4)
      !
      ! Now  calculate VNC
      !##    if (needvtau) then
      !##aux=0.d0
      !##call derivative(Grid,PAW%tvtau,aux,1,nr)
      !##write(std_out,*) ' In Kerker -- aux ', aux(1:10)
      !##aux(1:10)=0.d0
      !##endif     
      OPEN(88,file='NC',form='formatted')
      !
      VNC=0
      DO  i=1,nr
         x=(r(i)/rc)
         p(i)=Coef0+(x**m(1))*Coef(1)+(x**m(2))*Coef(2)+&
&             (x**m(3))*Coef(3)+(x**m(4))*Coef(4)
         dpp=(m(1)*(x**(m(1)-1))*Coef(1)+m(2)*(x**(m(2)-1))*Coef(2)+&
&             m(3)*(x**(m(3)-1))*Coef(3)+m(4)*(x**(m(4)-1))*Coef(4))/rc
         ddpp=(m(1)*(m(1)-1)*(x**(m(1)-2))*Coef(1)+&
&             m(2)*(m(2)-1)*(x**(m(2)-2))*Coef(2)+&
&             m(3)*(m(3)-1)*(x**(m(3)-2))*Coef(3)+&
&             m(4)*(m(4)-1)*(x**(m(4)-2))*Coef(4))/(rc**2)
         dddpp=(m(1)*(m(1)-1)*(m(1)-2)*(x**(m(1)-3))*Coef(1)+&
&             m(2)*(m(2)-1)*(m(2)-2)*(x**(m(2)-3))*Coef(2)+&
&             m(3)*(m(3)-1)*(m(3)-2)*(x**(m(3)-3))*Coef(3)+&
&             m(4)*(m(4)-1)*(m(4)-2)*(x**(m(4)-3))*Coef(4))/(rc**3)
         IF (i==irc) THEN
            WRITE(STD_OUT,*) 'check  dp ', dpp,  B/rc
            WRITE(STD_OUT,*) 'check ddp ', ddpp, C/rc**2
            WRITE(STD_OUT,*) 'check dddp', dddpp,  D/rc**3
         ENDIF
         IF (wavetype==EXPF) THEN
           if (needvtau) then
               VNC(i)=e+(1.d0+PAW%tvtau(i))*(ddpp+dpp*(dpp+2*(l+1)/r(i))) &
&                +PAW%dtvtaudr(i)*(dpp+l/r(i))                       
           else        
               VNC(i)=e+ddpp+dpp*(dpp+2*(l+1)/r(i))
           endif
            dum(i)=(r(i)**(l+1))*EXP(p(i))
         ENDIF
         IF (wavetype==POLY) THEN
           if (needvtau) then
               VNC(i)=e+(1.d0+PAW%tvtau(i))*(ddpp+2*(l+1)*dpp/r(i))/p(i) &
&                ++PAW%dtvtaudr(i)*(dpp/p(i)+l/r(i))                       
           else        
               VNC(i)=e+(ddpp+2*(l+1)*dpp/r(i))/p(i)
           endif
            dum(i)=(r(i)**(l+1))*(p(i))
         ENDIF
         WRITE(88,'(1p,5e15.7)') r(i),wfn(i),dum(i),VNC(i)*r(i),rv(i)
      ENDDO
      CLOSE(88)
      x=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
      WRITE(STD_OUT,*) 'check norm ',x,S

      VNC(irc:n)=rv(irc:n)/r(irc:n)
      PAW%rveff(1:n)=VNC(1:n)*r(1:n)

      DEALLOCATE(VNC,wfn,p,dum,aux)
    END SUBROUTINE kerker

  !***************************************************************
  ! SUBROUTINE VPSmatch(lmax,Grid,Pot,NC)
  !  Creates  screened norm-conserving pseudopotential similar to the
  !    approach of N. Troullier and J. L. Martins, PRB 43, 1993 (1991)
  !    Uses p(r)=a0+f(r); f(r)=SUMm(Coef(m)*r^(2*m), where
  !          m=1,2..6
  !    Psi(r) = r^(l+1)*exp(p(r))
  !    Modified for MGGA case and norm conserving condition is optional
  !    norm conservation controlled with optional variable NC
  !      defaults to no norm conservation
  !  Note this program assumes that wfn keeps the same sign for
  !    all matching points r(irc-match+1)....r(irc)
  !***************************************************************
  SUBROUTINE VPSmatch(Grid,Pot,PAW,l,e,NC)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
    INTEGER,INTENT(IN) :: l
    REAL(8),INTENT(IN) :: e
    LOGICAL,INTENT(IN), OPTIONAL :: NC

    REAL(8), ALLOCATABLE :: VNC(:)
    REAL(8) :: A0,A,B,B0,C,C0,D,F,S
    REAL(8) :: Coef(6),Coef0,Coef0old,corr
    REAL(8) :: h,rc,delta,x
    REAL(8) :: gam,bet
    INTEGER :: i,j,k,n,iter,nr,nodes,irc,ok,m,wavetype
    INTEGER, parameter :: match=6
    REAL(8) :: AAA(match,match),BBB(match)
    REAL(8) :: AAAA(match-1,match-1),BBBB(match-1)
    INTEGER, PARAMETER :: niter=5000
    REAL(8), PARAMETER :: small=1.0d-9
    REAL(8), ALLOCATABLE :: wfn(:),p(:),dum(:),aux(:),Kaux(:),v(:),dp(:),ddp(:)
    REAL(8), POINTER :: r(:),rv(:)
    CHARACTER(132) :: line
    LOGICAL :: normcons

    normcons=.false.
    if(PRESENT(NC)) normcons=NC

    write(std_out,*) 'Entering VPSmatch with normcons ', normcons

    n=Grid%n
    h=Grid%h
    r=>Grid%r
    rv=>Pot%rv
    nr=min(PAW%irc_vloc+10,n)
    irc=PAW%irc_vloc
    rc=PAW%rc_vloc

    ALLOCATE(VNC(n),wfn(n),p(n),dum(n),aux(n),dp(n),ddp(n),Kaux(n),v(n),stat=ok)
    IF (ok /=0) THEN
       WRITE(std_out,*) 'Error in VPSmatch  -- in allocating wfn,p', n,ok
       STOP
    ENDIF

    wfn=0
    !write(std_out,*) ' VPSmatch ', n,nr,irc,match
    !call flush_unit(std_out)
    if (scalarrelativistic) then
       CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
    else if (needvtau) then
       CALL unboundked(Grid,Pot,nr,l,e,wfn,nodes)
    else
       CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,e,wfn,nodes)
    endif

    IF (wfn(irc)<0) wfn=-wfn

!  first solve non-norm conserving results

    AAA=0;BBB=0
    Do i=1,match
       x=r(irc-match+i)
       BBB(i)=log(wfn(irc-match+i)/(x**(l+1)))
       do j=1,match
          AAA(i,j)=x**(2*(j-1))
       enddo
    Enddo

    CALL  SolveAXeqB(match,AAA,BBB,1.d20)

    write(std_out,*) 'Returned from SolveAXeqB in VPSmatch with Coefficients '
    write(std_out,'(1p,50e16.7)') (BBB(i),i=1,match)
    ! Now  calculate VNC
    !##    if (needvtau) then
    !##aux=0.d0
    !##call derivative(Grid,PAW%tvtau,aux,1,nr)
    !##Kaux=0.d0
    !##call derivative(Grid,PAW%Ktvtau,Kaux,1,nr)   ! Kresse form
    !##
    !##write(std_out,*) 'tvtau slope near origin -- Bloechl  ',aux(1:10) 
    !##write(std_out,*) 'tvtau slope near origin -- Kresse  ',Kaux(1:10) 
    !##write(std_out,*) 'reset to zero'
    !##aux(1:10)=0.d0; Kaux(1:10)=0.d0
    !##
    !##endif     
    OPEN(88,file='NNC',form='formatted')
    !
    VNC=0;p=0;dp=0;ddp=0;dum=0
!specific for match=6
    DO  i=2,nr
       x=(r(i))**2
       p(i)=BBB(1)+x*(BBB(2)+x*(BBB(3)+x*(BBB(4)+x*(BBB(5)+x*BBB(6)))))
       dp(i)=2*r(i)*(BBB(2)+x*(2*BBB(3)+x*(3*BBB(4)+x*(4*BBB(5)+5*x*BBB(6)))))
       ddp(i)=2*(BBB(2)+x*(6*BBB(3)+x*(15*BBB(4)+x*(28*BBB(5)+45*x*BBB(6)))))
       if (needvtau) then
               VNC(i)=e+(1.d0+PAW%tvtau(i))*(ddp(i)+ &
&                dp(i)*(dp(i)+2*(l+1)/r(i))) &
&                +PAW%dtvtaudr(i)*(dp(i)+l/r(i))                       
               v(i)=e+(1.d0+PAW%Ktvtau(i))*(ddp(i)+ &
&                dp(i)*(dp(i)+2*(l+1)/r(i))) &
&                +PAW%Kdtvtaudr(i)*(dp(i)+l/r(i))                       
       else        
               VNC(i)=e+ddp(i)+dp(i)*(dp(i)+2*(l+1)/r(i))
               v(i)=e+ddp(i)+dp(i)*(dp(i)+2*(l+1)/r(i))
       endif
               dum(i)=(r(i)**(l+1))*EXP(p(i))
    ENDDO
 
    S=overlap(Grid,wfn(1:irc),wfn(1:irc),1,irc)
    C0=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
    WRITE(STD_OUT,*) 'check norm ',C0,S

    VNC(irc:n)=rv(irc:n)/r(irc:n)
    v(irc:n)=rv(irc:n)/r(irc:n)
    PAW%rveff(1:n)=VNC(1:n)*r(1:n)
    PAW%Krveff(1:n)=v(1:n)*r(1:n)

    write(88,*) '#r   wfn and stuff   VPS(B)    r*VPS(B)    r*VPS(K)   r*v' 
    DO  i=1,n
       WRITE(88,'(1p,50e15.7)') r(i),dum(i),wfn(i),VNC(i),PAW%rveff(i),PAW%Krveff(i),rv(i)
    ENDDO
    CLOSE(88)
    
    write(std_out,*) 'Completed non-norm-conserving PS '

   if(.not.normcons) return

!  Iterate to find norm conserving results

    C0=C0/(EXP(2*BBB(1)))
    delta=1.d10
    iter=0
    Coef0=0.5d0*log(S/C0)

    DO WHILE(delta>small.AND.iter<=niter)
       Coef0old=Coef0
       iter=iter+1
       AAAA=0;BBBB=0
       Do i=1,match-1
          x=r(irc-match+1+i)
          BBBB(i)=log(wfn(irc-match+1+i)/(x**(l+1)))-Coef0old
          do j=1,match-1
             AAAA(i,j)=x**(2*(j))
          enddo
       Enddo

       CALL  SolveAXeqB(match-1,AAAA,BBBB,1.d20)

       write(std_out,*) 'Returned from SolveAXeqB in VPSmatch with Coefficients '
       write(std_out,'(1p,50e16.7)') (BBBB(i),i=1,match-1)

       !specific for match-1=5
           p=0;dum=0
           DO  i=2,nr
              x=(r(i))**2
              p(i)=x*(BBBB(1)+x*(BBBB(2)+x*(BBBB(3)+x*(BBBB(4)+x*BBBB(5)))))
              dum(i)=(r(i)**(l+1))*EXP(Coef0+p(i))
           ENDDO
           C0=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
           WRITE(std_out,*) 'check norm ',C0,S
           corr=0.5d0*log(S/C0);Coef0=Coef0+corr
           delta=ABS(Coef0-Coef0old)
           WRITE(std_out,'(" VNC: iter Coef0 delta",i5,1p,3e15.7)') iter,Coef0,delta,corr
    ENDDO
    WRITE(std_out,*) '  VNC converged in ', iter,'  iterations'
    WRITE(std_out,*) '  Coefficients  -- ', Coef0,BBBB(1:match-1)
    !
    ! Now  calculate VNC
    !##    if (needvtau) then
    !##aux=0.d0
    !##call derivative(Grid,PAW%tvtau,aux,1,nr)
    !##Kaux=0.d0
    !##call derivative(Grid,PAW%Ktvtau,Kaux,1,nr)
    !##
    !##write(std_out,*) 'tvtau slope near origin -- Bloechl  ',aux(1:10) 
    !##write(std_out,*) 'tvtau slope near origin -- Kresse  ',Kaux(1:10) 
    !##write(std_out,*) 'reset to zero'
    !##aux(1:10)=0.d0; Kaux(1:10)=0.d0
    !##
    !##endif     
    OPEN(88,file='NC',form='formatted')
    !
    VNC=0;p=0;dp=0;ddp=0;v=0
    
    DO  i=2,nr
       x=(r(i))**2
       p(i)=x*(BBBB(1)+x*(BBBB(2)+x*(BBBB(3)+x*(BBBB(4)+x*BBBB(5)))))
       dp(i)=2*r(i)*(BBBB(1)+x*(2*BBBB(2)+x*(3*BBBB(3)+x*(4*BBBB(4)+5*x*BBBB(5)))))
       ddp(i)=2*(BBBB(1)+x*(6*BBBB(2)+x*(15*BBBB(3)+x*(28*BBBB(4)+45*x*BBBB(5)))))
       if (needvtau) then
               VNC(i)=e+(1.d0+PAW%tvtau(i))*(ddp(i)+ &
&                dp(i)*(dp(i)+2*(l+1)/r(i))) &
&                +PAW%dtvtaudr(i)*(dp(i)+l/r(i))                       
               v(i)=e+(1.d0+PAW%Ktvtau(i))*(ddp(i)+ &
&                dp(i)*(dp(i)+2*(l+1)/r(i))) &
&                +PAW%Kdtvtaudr(i)*(dp(i)+l/r(i))                       
       else        
               VNC(i)=e+ddp(i)+dp(i)*(dp(i)+2*(l+1)/r(i))
               v(i)=e+ddp(i)+dp(i)*(dp(i)+2*(l+1)/r(i))
       endif
               dum(i)=(r(i)**(l+1))*EXP(Coef0+p(i))
    ENDDO
 
    C0=overlap(Grid,dum(1:irc),dum(1:irc),1,irc)
    WRITE(std_out,*) 'check norm ',C0,S

    VNC(irc:n)=rv(irc:n)/r(irc:n)
    v(irc:n)=rv(irc:n)/r(irc:n)
    PAW%rveff(1:n)=VNC(1:n)*r(1:n)
    PAW%Krveff(1:n)=v(1:n)*r(1:n)

    write(88,*) '#r   wfn and stuff   VPS(B)    r*VPS(B)    r*VPS(K)   r*v' 
    DO  i=1,n
       WRITE(88,'(1p,50e15.7)') r(i),dum(i),wfn(i),VNC(i),PAW%rveff(i),PAW%Krveff(i),rv(i)
    ENDDO
    CLOSE(88)
    
    write(std_out,*) 'Completed norm-conserving PS '

    if (iter.ge.niter) then
         write(std_out,*) 'Failed to converged norm-conserving VPS '
         stop
    endif     
    DEALLOCATE(VNC,wfn,p,dp,ddp,dum,aux,Kaux,v)
  END SUBROUTINE VPSmatch


  !***************************************************************
  ! SUBROUTINE nonncps(lmax,Grid,Pot)
  !  Creates  screened pseudopotential by inverting Schroedinger
  !    equation from a pseudized radial wave function of the form:
  !        Psi(r) = r**(l+1) * exp (a + b*r**2 + c*r**4 + d*r**6)
  !  No norm-conserving condition is imposed on Psi
  !***************************************************************
  SUBROUTINE nonncps(Grid,Pot,PAW,l,e)
    TYPE(Gridinfo), INTENT(IN) :: Grid
    TYPE(Potentialinfo), INTENT(IN) :: Pot
    TYPE(Pseudoinfo), INTENT(INOUT) ::  PAW
    INTEGER,INTENT(IN) :: l
    REAL(8),INTENT(IN) :: e

    INTEGER :: i,irc,n,nr,ok,nodes,i1,i2,i3,i4
    REAL(8) :: rc,x,y1,y2,y3,p0,p1,p2,p3,sgn
    REAL(8) :: b(4),c(4),d(4),amat(4,4)
    REAL(8),ALLOCATABLE ::  VNC(:),wfn(:),aux(:),Kaux(:)
    REAL(8),POINTER :: r(:),rv(:)
    CHARACTER(132) :: line

   !Polynomial definitions
    p0(x,y1,y2,y3)=(x-y1)*(x-y2)*(x-y3)
    p1(x,y1,y2,y3)=(x-y2)*(x-y3)+(x-y1)*(x-y3)+(x-y1)*(x-y2)
    p2(x,y1,y2,y3)=2.0d0*((x-y1)+(x-y2)+(x-y3))
    p3(x,y1,y2,y3)=6.0d0

    n=Grid%n
    r=>Grid%r
    rv=>Pot%rv
    nr=min(PAW%irc_vloc+10,n)
    irc=PAW%irc_vloc
    rc=PAW%rc_vloc

    ALLOCATE(VNC(n),wfn(nr),aux(nr),Kaux(nr),stat=ok)
    IF (ok/=0) stop 'Error in uspseudo -- allocating arrays'

    if (scalarrelativistic) then
       CALL unboundsr(Grid,Pot,nr,l,e,wfn,nodes)
    else if (needvtau) then
       CALL unboundked(Grid,Pot,nr,l,e,wfn,nodes)
    else
       CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,e,wfn,nodes)
    endif
    IF (wfn(irc)<0) wfn=-wfn

    DO i=2,nr
      wfn(i)=wfn(i)/r(i)**dble(l+1)
    ENDDO

    i1=irc-1;i2=i1+1;i3=i2+1;i4=i3+1
    c(1)=wfn(i1)/p0(r(i1),r(i2),r(i3),r(i4))
    c(2)=wfn(i2)/p0(r(i2),r(i3),r(i4),r(i1))
    c(3)=wfn(i3)/p0(r(i3),r(i4),r(i1),r(i2))
    c(4)=wfn(i4)/p0(r(i4),r(i1),r(i2),r(i3))
    d(1)=c(1)*p0(rc,r(i2),r(i3),r(i4)) + c(2)*p0(rc,r(i3),r(i4),r(i1)) + &
&        c(3)*p0(rc,r(i4),r(i1),r(i2)) + c(4)*p0(rc,r(i1),r(i2),r(i3))
    d(2)=c(1)*p1(rc,r(i2),r(i3),r(i4)) + c(2)*p1(rc,r(i3),r(i4),r(i1)) + &
&        c(3)*p1(rc,r(i4),r(i1),r(i2)) + c(4)*p1(rc,r(i1),r(i2),r(i3))
    d(3)=c(1)*p2(rc,r(i2),r(i3),r(i4)) + c(2)*p2(rc,r(i3),r(i4),r(i1)) + &
&        c(3)*p2(rc,r(i4),r(i1),r(i2)) + c(4)*p2(rc,r(i1),r(i2),r(i3))
    d(4)=c(1)*p3(rc,r(i2),r(i3),r(i4)) + c(2)*p3(rc,r(i3),r(i4),r(i1)) + &
&        c(3)*p3(rc,r(i4),r(i1),r(i2)) + c(4)*p3(rc,r(i1),r(i2),r(i3))

    sgn=d(1)/abs(d(1));d(1:4)=d(1:4)*sgn
    b(1)=log(d(1));b(2:4)=d(2:4)
    amat(1,1)= 1.0d0
    amat(2:4,1)= 0.0d0
    amat(1,2)= rc**2
    amat(2,2)= 2.0d0*d(1)*rc
    amat(3,2)= 2.0d0*d(1)   +2.0d0*d(2)*rc
    amat(4,2)=               4.0d0*d(2)   +2.0d0*d(3)*rc
    amat(1,3)= rc**4
    amat(2,3)=  4.0d0*d(1)*rc**3
    amat(3,3)= 12.0d0*d(1)*rc**2+ 4.0d0*d(2)*rc**3
    amat(4,3)= 24.0d0*d(1)*rc   +24.0d0*d(2)*rc**2+4.0d0*d(3)*rc**3
    amat(1,4)= rc**6
    amat(2,4)=   6.0d0*d(1)*rc**5
    amat(3,4)=  30.0d0*d(1)*rc**4+ 6.0d0*d(2)*rc**5
    amat(4,4)= 120.0d0*d(1)*rc**3+60.0d0*d(2)*rc**4+6.0d0*d(3)*rc**5

    CALL linsol(amat,b,4,4,4,4)
    !write(std_out,*) 'Completed linsol with coefficients'
    !write(std_out,'(1p,10e15.7)') (b(i),i=1,4)

    !##    if (needvtau) then
    !##  aux=0.d0
    !##  call derivative(Grid,PAW%tvtau,aux,1,nr)
    !##  Kaux=0.d0
    !##  call derivative(Grid,PAW%Ktvtau,Kaux,1,nr)
    !##
    !##  write(std_out,*) 'tvtau slope near origin -- Bloechl  ',aux(1:10) 
    !##  write(std_out,*) 'tvtau slope near origin -- Kresse  ',Kaux(1:10) 
    !##  write(std_out,*) 'reset to zero'
    !##  aux(1:10)=0.d0; Kaux(1:10)=0.d0
    !##
    !##endif     

    PAW%rveff(1)=0.d0
    DO i=2,irc-1
     c(1)=2.0d0*b(2)*r(i)+ 4.0d0*b(3)*r(i)**3+ 6.0d0*b(4)*r(i)**5
     c(2)=2.0d0*b(2)     +12.0d0*b(3)*r(i)**2+30.0d0*b(4)*r(i)**4
     if (needvtau) then
        PAW%rveff(i)=r(i)*(e+(dble(2*l+2)*c(1)/r(i)+c(1)**2+&
&           c(2))*(1.d0+PAW%tvtau(i))++PAW%dtvtaudr(i)*(c(1)+dble(l)/r(i)))
        PAW%Krveff(i)=r(i)*(e+(dble(2*l+2)*c(1)/r(i)+c(1)**2+&
&           c(2))*(1.d0+PAW%Ktvtau(i))++PAW%Kdtvtaudr(i)*(c(1)+dble(l)/r(i)))
     else        
        PAW%rveff(i)=r(i)*(e+dble(2*l+2)*c(1)/r(i)+c(1)**2+c(2))
        PAW%Krveff(i)=r(i)*(e+dble(2*l+2)*c(1)/r(i)+c(1)**2+c(2))
     endif
    ENDDO
    PAW%rveff(irc:n)=rv(irc:n)
    PAW%Krveff(irc:n)=rv(irc:n)


    open(88,file='NNC',form='formatted')
    write(88,*) '   r   rv(AE)    rps(B)        rps(K)  '
    do i=1,n
       write(88,'(1p,50e16.7)') r(i),rv(i),PAW%rveff(i),PAW%Krveff(i)
    enddo
    close(88)   

    DEALLOCATE(VNC,wfn,aux,Kaux)

  END SUBROUTINE nonncps



    SUBROUTINE checkghosts(Grid,Orbit,FC,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(OrbitInfo), INTENT(IN) :: Orbit
      TYPE(FCInfo), INTENT(IN) :: FC
      TYPE(PseudoInfo), INTENT(in) :: PAW
      INTEGER :: l,nr,nodes,i,io
      REAL(8) :: energy,h
      REAL(8), POINTER :: r(:)
      REAL(8), ALLOCATABLE :: wfn(:),VNC(:)
      TYPE(PotentialInfo) :: Pot

      h=Grid%h
      r=>Grid%r
      nr=min(PAW%irc+5,Grid%n)
      call InitPot(Pot,nr)
      ALLOCATE(VNC(nr),wfn(nr),stat=i)
      IF (i/=0) STOP 'Error in checkghosts allocation'
      POT%rv(1:nr)=PAW%rveff(1:nr)
      call zeropot(Grid,POT%rv,POT%v0,POT%v0p)

      DO l=0,PAW%lmax
         DO io=1,Orbit%norbit
            IF((.NOT.Orbit%iscore(io)).AND.(Orbit%l(io)==l)) THEN
               energy=Orbit%eig(io)
               WRITE(STD_OUT,*) 'Check  for ghosts with  l', l,energy
               CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,energy,wfn,nodes)
               !DO i=1,nr
               !   WRITE(l+17,'(1p,2e15.7)') Grid%r(i),wfn(i)
               !ENDDO
               WRITE(STD_OUT,*) '    Found # nodes = ', nodes
               EXIT
            ENDIF
         ENDDO
      ENDDO

      CALL DestroyPot(Pot)
      DEALLOCATE(VNC,wfn)

    END SUBROUTINE checkghosts

    SUBROUTINE sethat(Grid,PAW,gaussparam,besselopt)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER,INTENT(IN), OPTIONAL :: besselopt
      REAL(8),INTENT(IN), OPTIONAL :: gaussparam

      INTEGER :: n,irc,irc_shap,i
      REAL(8), POINTER :: r(:)
      REAL(8) :: h,con,rc,rc_shap,selfen,d,dd,jbes1,jbes2,qr
      REAL(8) :: al(2),ql(2)

      n=Grid%n
      h=Grid%h
      irc=PAW%irc
      rc=PAW%rc
      irc_shap=PAW%irc_shap
      rc_shap=PAW%rc_shap
      r=>Grid%r

      PAW%hatden=0
      PAW%projshape=0
      PAW%hatshape=0
      PAW%projshape(1)=1
      PAW%hatshape(1)=1
      DO i=2,irc-1
       PAW%projshape(i)=(SIN(pi*r(i)/rc)/(pi*r(i)/rc))**2
      ENDDO
      if(present(gaussparam)) then
       d=rc_shap/SQRT(LOG(1.d0/gaussparam))
       PAW%gausslength=d
       DO i=2,irc
        PAW%hatshape(i)=EXP(-(r(i)/d)**2)
       ENDDO
       PAW%irc_shap=PAW%irc
       PAW%rc_shap=PAW%rc
      else if(present(besselopt)) then
       call shapebes(al,ql,0,rc_shap)
       DO i=1,irc_shap-1
        qr=ql(1)*r(i);CALL jbessel(jbes1,d,dd,0,0,qr)
        qr=ql(2)*r(i);CALL jbessel(jbes2,d,dd,0,0,qr)
        PAW%hatshape(i)=al(1)*jbes1+al(2)*jbes2
       ENDDO
      else
       DO i=2,irc_shap-1
        PAW%hatshape(i)=(SIN(pi*r(i)/rc_shap)/(pi*r(i)/rc_shap))**2
       ENDDO
      endif
      PAW%hatden(1:irc)=PAW%hatshape(1:irc)*(r(1:irc)**2)

      !  normalize
      if (.not.besselshapefunction) then
       con=integrator(Grid,PAW%hatden,1,PAW%irc_shap)
       WRITE(STD_OUT,*) ' check hatden normalization', con
       PAW%hatden=PAW%hatden/con
      endif

      CALL poisson(Grid,con,PAW%hatden,PAW%hatpot,selfen)
      WRITE(STD_OUT,*) 'Self energy for L=0 hat density  ', selfen
      WRITE(STD_OUT,*) 'hatden charge  ', con

      OPEN(1001,file='hatpot',form='formatted')
      write(1001,'(a)') '#     r    hatden       hatpot'
      do i=1,n
           write(1001,'(1p, 3e15.7)') r(i), PAW%hatden(i), PAW%hatpot(i)
      enddo   
      CLOSE(1001)


    END SUBROUTINE sethat

    SUBROUTINE coretailselfenergy(Grid,PAW,ctctse,cthatse)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      REAL(8), INTENT(OUT) :: ctctse,cthatse

      INTEGER :: i,irc,n
      REAL(8) :: rc,h,x,y,z
      REAL(8), allocatable :: d1(:),d2(:)

      n=Grid%n
      h=Grid%h
      irc=PAW%irc

      allocate(d1(n),d2(n),stat=i)
          if (i /= 0) then
             write(std_out,*) 'coretailselfenergy: allocation error -- ', n,i
             stop
          endif

      x=integrator(Grid,PAW%tcore)
      write(std_out,*) 'tcore charge ' , x
      CALL poisson(Grid,x,PAW%tcore,d1,ctctse)
      d2(2:n)=PAW%hatden(2:n)*d1(2:n)/Grid%r(2:n)
      d2(1)=0
      cthatse=integrator(Grid,d2(1:irc),1,PAW%irc_shap)
      write(std_out,*) 'ctctse,cthatse = ', ctctse,cthatse

      deallocate(d1,d2)

    END SUBROUTINE coretailselfenergy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SUBROUTINE smoothcore
!    program to take array orig (all electron den or tau
!        and return smooth polynomial function for 0 \le r \le rc_core
!        with first and second derivatives continuous
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE smoothcore(Grid,orig,smooth)
      TYPE(GridInfo), INTENT(IN) :: Grid
      REAL(8), INTENT(IN) :: orig(:)
      REAL(8), INTENT(INOUT) :: smooth(:)

      REAL(8) :: rc,h,x,y,z,u0,u2,u4
      REAL(8), allocatable :: d1(:),d2(:)
      INTEGER :: i,j,k,n,irc
      INTEGER :: counter=0

      n=Grid%n
      h=Grid%h
      irc=PAW%irc_core
      rc=PAW%rc_core
      write(std_out,*) 'In smoothcore ', irc,rc

      allocate(d1(n),d2(n),stat=i)
          if (i /= 0) then
             write(std_out,*) 'smoothcore: allocation error -- ', n,i
             stop
          endif
      CALL derivative(Grid,orig,d1)
      CALL derivative(Grid,d1,d2)

      x=orig(irc)
      y=d1(irc)*rc
      z=d2(irc)*(rc*rc)
      write(std_out,*) 'smoothcore: x,y,z = ', x,y,z

      u0=3*x - 9*y/8 + z/8
      u2=-3*x + 7*y/4 - z/4
      u4=x - 5*y/8 + z/8

      write(std_out,*) 'smoothcore: u0,u2,u4 = ', u0,u2,u4

      smooth=orig
      do i=1,irc
         x=(Grid%r(i)/rc)**2
         smooth(i)= x*(u0+x*(u2+x*u4))
      enddo

         open(1001,file='smoothcore',form='formatted')
         do i=1,Grid%n
          write(1001,'(1p,50e15.7)') Grid%r(i), orig(i),smooth(i)
         enddo  
         close(1001)  

      deallocate(d1,d2)
  END SUBROUTINE smoothcore

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SUBROUTINE smoothtau
!    program to take array orig (all electron) coretau
!        and return smooth polynomial function for 0 \le r \le rc_core
!        similar to smoothcore, but using r**4 as leading term
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE smoothtau(Grid,orig,smooth)
      TYPE(GridInfo), INTENT(IN) :: Grid
      REAL(8), INTENT(IN) :: orig(:)
      REAL(8), INTENT(INOUT) :: smooth(:)

      REAL(8) :: rc,h,x,y,z,u0,u2,u4
      REAL(8), allocatable :: d1(:),d2(:)
      INTEGER :: i,j,k,n,irc
      INTEGER :: counter=0

      n=Grid%n
      h=Grid%h
      irc=PAW%irc_core
      rc=PAW%rc_core
      write(std_out,*) 'In smoothtau ', irc,rc

      allocate(d1(n),d2(n),stat=i)
          if (i /= 0) then
             write(std_out,*) 'smoothtau: allocation error -- ', n,i
             stop
          endif
      CALL derivative(Grid,orig,d1)
      CALL derivative(Grid,d1,d2)

      x=orig(irc)
      y=d1(irc)*rc
      z=d2(irc)*(rc*rc)
      write(std_out,*) 'smoothtau: x,y,z = ', x,y,z

      u0=6*x - 13*y/8 + z/8
      u2=-8*x + 11*y/4 - z/4
      u4=3*x - 9*y/8 + z/8

      write(std_out,*) 'smoothtau: u0,u2,u4 = ', u0,u2,u4

      smooth=orig
      do i=1,irc
         x=(Grid%r(i)/rc)**2
         y=x*x
         smooth(i)= y*(u0+x*(u2+x*u4))
      enddo

         open(1001,file='smoothtau',form='formatted')
         do i=1,Grid%n
          write(1001,'(1p,50e15.7)') Grid%r(i), orig(i),smooth(i)
         enddo  
         close(1001)  

      deallocate(d1,d2)
  END SUBROUTINE smoothtau


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! SUBROUTINE smoothpower
!    program to take array orig (all electron tau or den)
!        and return smooth polynomial function for 0 \le r \le rc_core
!        matching 4 points less than and equal to rc_core
!          power could be 2 or 4, representing the leading power
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE smoothpower(Grid,power,orig,smooth)
      TYPE(GridInfo), INTENT(IN) :: Grid
      INTEGER, INTENT(IN) :: power
      REAL(8), INTENT(IN) :: orig(:)
      REAL(8), INTENT(INOUT) :: smooth(:)

      INTEGER, parameter :: terms=5
      REAL(8) :: rc,h,x,y,z,u0,u2,u4
      REAL(8) :: aa(terms,terms),Ci(terms)
      INTEGER :: i,j,k,n,irc,many

      n=Grid%n
      h=Grid%h
      irc=PAW%irc_core
      rc=Grid%r(irc)
      write(std_out,*) 'smoothpower -- ', power,irc,rc

      aa=0; Ci=0
      do i=1,terms
         x=Grid%r(irc-terms+i)
         Ci(i)=orig(irc-terms+i)/x**power
         do j=1,terms
            aa(i,j)=x**(2*(j-1))
         enddo   
      enddo
      call SolveAXeqBM(terms,aa,Ci,terms-1)
         write(std_out,*) 'Completed SolveAXeqB with coefficients'
         write(std_out,'(1p,10e15.7)') (Ci(i),i=1,terms)
      smooth=orig
      do i=1,irc-1
         smooth(i)=0
         x=Grid%r(i)
         do j=1,terms
            smooth(i)=smooth(i)+Ci(j)*(x**(power+2*(j-1)))
         enddo
      enddo

      if (power==2) open(1001,file='smoothcoreden',form='formatted')
      if (power==4) open(1001,file='smoothcoretau',form='formatted')
      do i=1,n
        write(1001,'(1p,50e16.7)') Grid%r(i),orig(i),smooth(i)
      enddo   
      close(1001)

  END SUBROUTINE smoothpower

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! SUBROUTINE smoothexpcore
!    program to take all electron coreden and coretau
!      and return simplified tcoreden and tcoretau in form r^2* C*exp(-alpha*r)
!                                                    for r>rc
!  Modified 10/20/2023 by NAWH to ensure Weizsacker form for tcoretau
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
    SUBROUTINE smoothexpcore(Grid,coredenin,coretauin,tcoredenout,tcoretauout)
      TYPE(GridInfo), INTENT(IN) :: Grid
      REAL(8), INTENT(IN) :: coredenin(:),coretauin(:)
      REAL(8), INTENT(INOUT) :: tcoredenout(:),tcoretauout(:)

      INTEGER :: i,j,k,n,irc
      REAL(8) :: x,rc,rc1,alpha,term1,term2,const
      REAL(8), allocatable :: tmpd(:),grad(:),dum(:)

      tcoredenout=coredenin
      tcoretauout=coretauin
      n=Grid%n
      allocate(tmpd(n),grad(n),dum(n))
      irc=PAW%irc_core
      rc=Grid%r(irc)
      rc1=Grid%r(irc-1)
      term1=tcoredenout(irc)/((rc**2))
      term2=tcoredenout(irc-1)/((rc1**2))
      x=(Grid%r(irc)-Grid%r(irc-1))
      alpha=-log(term1/term2)/x
      const=0.5d0*(term1*exp(alpha*rc)+term2*exp(alpha*rc1))
      write(std_out,*) 'smoothexpcore -- den ', irc,rc,const,alpha
      do i=1,irc-1
       tcoredenout(i)=(Grid%r(i)**2)*const*exp(-alpha*Grid%r(i))
      enddo 
!!! Insert exactly Weizacker form into augmentation range  using
!!!    same code as used in excor subroutines
        tmpd(2:n)=tcoredenout(2:n)/(4*pi*(Grid%r(2:n)**2))
        call extrapolate(Grid,tmpd)
        do i=1,n
           dum(i)=ddlog(tmpd(i))
        enddo
        call derivative(Grid,dum,grad,1,n)
        grad(1:n)=grad(1:n)*tmpd(1:n)     ! perhaps more accurate???
      do i=1,irc-1
         tcoretauout(i)=0.25d0*(4*pi*Grid%r(i)**2)*(grad(i)**2)/tmpd(i)        !!!Weizsacker form
      enddo 
         write(std_out,*) ' in smoothexpcore -- checking tcoretau mismatch '
           write(std_out,*) rc, tcoretauout(irc-1),coretauin(irc-1)
        ! take average at last point before rc:
          tcoretauout(irc-1)=0.5d0*(tcoretauout(irc-1)+coretauin(irc-1))

      open(1001,file='smoothexpcore.dat',form='formatted')
      write(1001,*) 'r             coreden           tcoreden          coretau       tcoretau'       
      do i=1,n
        write(1001,'(1p,50e16.7)') Grid%r(i),coredenin(i),tcoredenout(i),coretauin(i),tcoretauout(i)
      enddo   
      close(1001)
      deallocate(tmpd,grad,dum)


    END SUBROUTINE smoothexpcore         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE setcoretail(Grid,coreden)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE setcoretail(Grid,coreden)
      TYPE(GridInfo), INTENT(IN) :: Grid
      REAL(8), INTENT(IN) :: coreden(:)

      REAL(8) :: rc,h,x,y,z,u0,u2,u4
      INTEGER :: i,j,k,n,irc

      n=Grid%n
      h=Grid%h
      irc=PAW%irc_core
      rc=PAW%rc_core

      ! changed 6/6/2023 by NAWH to original smoothcore if .not.needvtau
      !!  changed again for possible alternate coreshape
      If(.not.needvtau) then
              CALL smoothcore(Grid,coreden,PAW%tcore) 
      else if (.not.PAW%tcoreshapeexp) then       
              CALL smoothpower(Grid,2,coreden,PAW%tcore) 
      endif
      PAW%core=coreden

  ! Find coretailpoints
     z = integrator(Grid,coreden)

     coretailpoints=PAW%irc+Grid%ishift    !! coretailpoints should be>=PAW%irc
        do i=PAW%irc+Grid%ishift,n
           if(ABS(z-integrator(Grid,coreden,1,i))<coretailtol) then
             coretailpoints=i
             exit
           endif
        enddo
     write(std_out,*) 'coretailpoints = ',coretailpoints
     PAW%coretailpoints=coretailpoints

    END SUBROUTINE setcoretail
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  resettcore(Grid,PAW)
!     This routine attempts to construct tcore so that it is likely
!      to cancel a negative compensation charge for the purpose of
!      dealing with Kresse exchange-correlation form
!     Performs the matching between irc_core and irc_shape, checking
!      to make sure that tcore is positive.   Stops if no positive
!      tcore is found
!  Note: assumes that PAW_irc_core<=irc_shape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE resettcore(Grid,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      
      REAL(8), allocatable :: d(:),g(:)
      REAL(8) :: denrc,ddenrc,shaperc,dshaperc
      INTEGER :: i,j
      LOGICAL :: ok

      allocate(d(Grid%n),g(Grid%n))
      call derivative(Grid,PAW%core,d)      
      call derivative(Grid,PAW%hatden,g)      

      Do j=PAW%irc_core,PAW%irc_shap
         denrc=PAW%core(j)
         ddenrc=d(j)
         shaperc=PAW%hatden(j)
         dshaperc=g(j)
         PAW%tcore=PAW%core
         ok=.true.
         do i=1,j-1
           PAW%tcore(i)=denrc+(ddenrc/dshaperc)*(PAW%hatden(i)-shaperc)
           if (PAW%tcore(i)<0.d0) ok=.false.
           if (abs(PAW%tcore(i))<machine_zero) PAW%tcore(i)=0.d0
         enddo      
         write(std_out,*) 'resettcore at j =', j,' r = ', Grid%r(j), ' good?  ',ok
         if (ok) exit
      enddo     

      open(2005,file='tcoreden.dat',form='formatted')
      do i=1,Grid%n
         write(2005,'(1p,30e16.7)')Grid%r(i),PAW%tcore(i),PAW%core(i), &
&                  PAW%hatden(i)    
      enddo
      close(2005)

      if (.not.ok) then
          write(std_out,*) 'resettcore terminating due to negative values '
          stop
      else  
          write(std_out,*)' successfully completed resettcore with irc_core = ',j    
      endif    
! presumably it is not necessary to recalculate coretailpoints
     deallocate(d,g)
End subroutine resettcore    

    SUBROUTINE setttau(Grid,coreden,coretau)
      TYPE(GridInfo), INTENT(IN) :: Grid
      REAL(8), INTENT(IN) :: coreden(:),coretau(:)
      REAL(8) :: sqr4pi

      write(std_out,*) 'in setttau '
      !CALL smoothtau(Grid,coretau,PAW%tcoretau) 

      if (PAW%tcoreshapeexp) then
        CALL smoothexpcore(Grid,coreden,coretau,PAW%tcore,PAW%tcoretau)     
      else        
        CALL smoothpower(Grid,4,coretau,PAW%tcoretau) 
      endif

      PAW%coretau=coretau
      write(std_out,*) 'completed setttau'

      sqr4pi=sqrt(4*pi)*1.d-10
      PAW%itau=max(PAW%coretailpoints,1)
      do while (PAW%itau<Grid%n.and. &
&       abs(PAW%tcoretau(PAW%itau))>sqr4pi*Grid%r(PAW%itau)**2)
       PAW%itau=PAW%itau+1
      end do

   END SUBROUTINE setttau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set smooth core functions for EXXKLI case  using procedure
!      identical to HF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE fixtcorewfn(Grid,PAW)
     TYPE(GridInfo), INTENT(IN) :: Grid
     TYPE(PseudoInfo), INTENT(INOUT) :: PAW

     INTEGER :: i,j,k,l,n,io,jo,ib,ic,irc
     LOGICAL :: last
     INTEGER , PARAMETER :: np=5

     n=Grid%n
     irc=PAW%irc

     !Do i=1,n
     !   write(444,'(1p,10e15.7)') Grid%r(i),PAW%tcore(i)
     !enddo
     Do io=1,PAW%TOCCWFN%norbit
        write(std_out,*) 'fixtcore', io,PAW%valencemap(io),PAW%TOCCWFN%iscore(io)
        IF(PAW%valencemap(io)<0) THEN     ! core state
          PAW%TOCCWFN%wfn(:,io)=0
          If (PAW%TOCCWFN%iscore(io)) then
             last=.true. ; l=PAW%TOCCWFN%l(io)
           ! check to find if this is the outermost core state for this l
             if (io<PAW%TOCCWFN%norbit) then
                do jo=io+1,PAW%TOCCWFN%norbit
                   if (PAW%TOCCWFN%iscore(jo).and.PAW%TOCCWFN%l(jo)==l) &
&                        last=.false.
                enddo
             endif
             if (last) then
                IF (MAXVAL(ABS(PAW%OCCWFN%wfn(irc-np/2:irc+np/2,io))) &
&                      >PAW%coretol) THEN
                   CALL Smoothfunc(Grid%r,l,np,irc,n,PAW%OCCWFN%wfn(:,io), &
&                          PAW%TOCCWFN%wfn(:,io))
                write(std_out,*) 'last',io
                ENDIF
             endif
          else
              write(std_out,*) 'Something wrong -- should be core state ',&
&                io,PAW%TOCCWFN%l(io),PAW%TOCCWFN%eig(io),PAW%TOCCWFN%occ(io)
              stop
          endif
       ENDIF
    Enddo

   ! reset PAW%tcore to be consistent
   PAW%tcore=0
   do io=1,PAW%TOCCWFN%norbit
      if (PAW%TOCCWFN%iscore(io)) then
         PAW%tcore=PAW%tcore+PAW%TOCCWFN%occ(io)*(PAW%TOCCWFN%wfn(:,io)**2)
      endif
   enddo
     !Do i=1,n
     !   write(445,'(1p,10e15.7)') Grid%r(i),PAW%tcore(i)
     !enddo
 END  SUBROUTINE fixtcorewfn

    SUBROUTINE selfhatpot(Grid,PAW,l,eself)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER, INTENT(IN) :: l
      REAL(8), INTENT(OUT) :: eself


      INTEGER :: n,irc,i
      REAL(8), POINTER :: r(:)
      REAL(8), ALLOCATABLE :: den(:),a(:)
      REAL(8) :: h,con
      REAL(8) :: qr,jbes1,jbes2,dum1,dum2,al(2),ql(2)

      n=Grid%n
      h=Grid%h
      r=>Grid%r
      irc=PAW%irc

      ALLOCATE(den(n),a(n),stat=i)
      IF (i/=0) THEN
         WRITE(STD_OUT,*) 'Error in selfhatpot allocation',irc,i
         STOP
      ENDIF

      if (besselshapefunction) then
       call shapebes(al,ql,l,PAW%rc_shap)
       DO i=1,PAW%irc_shap
        qr=ql(1)*r(i);CALL jbessel(jbes1,dum1,dum2,l,0,qr)
        qr=ql(2)*r(i);CALL jbessel(jbes2,dum1,dum2,l,0,qr)
        den(i)=(al(1)*jbes1+al(2)*jbes2)*r(i)**2
       ENDDO
       if (n>PAW%irc_shap) den(PAW%irc_shap+1:n)=0.d0
      else
       DO i=1,n
        den(i)=(r(i)**l)*PAW%hatden(i)
       ENDDO
       a(1:n)=den(1:n)*(r(1:n)**l)
       con=integrator(Grid,a,1,PAW%irc_shap)
       den=den/con
      endif

      a=0

      CALL apoisson(Grid,l,n,den,a)

      ! apoisson returns a*r
      a(2:n)=a(2:n)/r(2:n)
      a(1)=0

      eself=0.5d0*overlap(Grid,a,den)

      DEALLOCATE(den,a)

    END SUBROUTINE selfhatpot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For diracrelativistic case, this subroutine is not yet ready
!     Suggested future changes are added as comments, such as
!     Upper component of wave function is loaded with renormalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setbasis(Grid,Pot,Orbit)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    TYPE(OrbitInfo), INTENT(INOUT) :: Orbit


    INTEGER :: n,irc,nbase,l,lmax,lng,currentnode,mxbase
    INTEGER :: i,j,k,io,ok,nbl,nr,nodes,ib,loop,niter,iter,ibasis_add
    REAL(8) :: h,rc,q00,energy,rat,delta,thisconv,qeff,tq,range
    REAL(8) :: ecoul,etxc,eexc
    REAL(8),allocatable :: checkden(:),valden(:)
    REAL(8), POINTER  :: r(:)

    n=Grid%n
    h=Grid%h
    r=>Grid%r
    irc=PAW%irc
    nr=MIN(irc+100,n-100)
    rc=PAW%rc
    lmax=PAW%lmax

    allocate(checkden(n),valden(n));checkden=0.d0;valden=0.d0

!Check beginning valence density
    DO io=1,Orbit%norbit 
       if (.not.Orbit%iscore(io)) then
         valden=valden+Orbit%occ(io)*(Orbit%wfn(:,io)**2)        
       endif
    ENDDO   

  ! set AErefrv
    PAW%AErefrv=Pot%rv
    PAW%rvx=Pot%rvx

    PAW%exctype=Orbit%exctype

    nbase=PAW%nbase
    mxbase=PAW%mxbase

  ! "filter" occupied states for long-range noise
    DO io=1,Orbit%norbit
       Call Filter(n,Orbit%wfn(:,io),machine_zero)
    ENDDO

    call CopyOrbit(Orbit,PAW%OCCwfn)
    call CopyOrbit(Orbit,PAW%TOCCwfn)
    PAW%valencemap=-13
    IF (Orbit%exctype=='HF') THEN
       range=r(n)
       rat=-1.d30; k=1
       do io=1,Orbit%norbit
          if ((Orbit%occ(io)>1.d-5).and.Orbit%eig(io)>rat) then
              rat=Orbit%eig(io)
              k=io
          endif
       enddo
       do i=n,irc+1,-1
          if (ABS(Orbit%wfn(i,k))>1.d-4) then
             range=r(i)
             exit
          endif
       enddo
       write(std_out,*) 'range ', k,range; call flush_unit(std_out)
    ENDIF

    WRITE(STD_OUT,*) '  basis functions:'
    WRITE(STD_OUT,*)' No.   n     l         energy         occ   '

    nbase=0 ; ibasis_add=1
    DO l=0,lmax
       currentnode=-1
       nbl=0
       DO io=1,Orbit%norbit    ! cycle through all configuration
          IF (Orbit%l(io).EQ.l) THEN
             currentnode=Orbit%np(io)-l-1     
            IF (.NOT.Orbit%iscore(io)) THEN
             nbl=nbl+1
             nbase=nbase+1
             PAW%np(nbase)=Orbit%np(io)
             PAW%l(nbase)=l
             PAW%nodes(nbase)=PAW%np(nbase)-l-1
             write(std_out,*) 'l,nbase,node',l,nbase,currentnode
             PAW%eig(nbase)=Orbit%eig(io)
             PAW%occ(nbase)=Orbit%occ(io)
             PAW%phi(:,nbase)=Orbit%wfn(:,io)

             if(diracrelativistic) then 
               STOP 'Error -- setbasis subroutine not ready for diracrelativistic!'
               !rat=1.d0/sqrt(overlap(Grid,PAW%phi(:,nbase),PAW%phi(:,nbase)))
               !PAW%phi(:,nbase)=rat*PAW%phi(:,nbase)
               !PAW%kappa(nbase)=kappa
             endif  

             PAW%valencemap(io)=nbase
             IF(Orbit%exctype=='HF') then
               PAW%eig(nbase)=HF%lmbd(io,io)
               PAW%lmbd(:,nbase)=HF%lmbd(:,io)
               write(std_out,*) 'lmbd for nbase ', nbase
               write(std_out,'(1p,20e15.7)') PAW%lmbd(1:Orbit%norbit,nbase)
               PAW%lmbd(io,nbase)=0.d0         !  to avoid double counting
             endif
             WRITE(STD_OUT,'(3i6,1p,2e15.6)') nbase,PAW%np(nbase),l,             &
&                 PAW%eig(nbase),PAW%occ(nbase); call flush_unit(std_out)
            ENDIF
          ENDIF
       ENDDO

       generalizedloop: DO
         IF (ibasis_add>input_dataset%nbasis_add) EXIT generalizedloop
         IF (input_dataset%basis_add_l(ibasis_add)/=l) EXIT generalizedloop
         energy=input_dataset%basis_add_energy(ibasis_add)
         IF (energy<0.d0) &
&            WRITE(STD_OUT,*) 'energy is negative',energy,' -- WARNING WARNING !!!'
         nbase=nbase+1
         IF (nbase > mxbase ) THEN
           WRITE(STD_OUT,*) 'Error in  setbasis -- too many functions ', nbase,mxbase
           STOP
         ENDIF
         PAW%l(nbase)=l
         PAW%np(nbase)=999
         PAW%nodes(nbase)=currentnode+1
         currentnode=PAW%nodes(nbase)
         write(std_out,*) 'l,nbase,node',l,nbase,currentnode
         PAW%eig(nbase)=energy
         PAW%occ(nbase)=0.d0
         PAW%phi(1:n,nbase)=0.d0
         if (scalarrelativistic) then
            CALL unboundsr(Grid,Pot,n,l,energy,PAW%phi(:,nbase),nodes)
         else if (needvtau) then
            CALL unboundked(Grid,Pot,n,l,energy,PAW%phi(:,nbase),nodes)
         elseif ( Orbit%exctype=='HF') then
            CALL HFunocc(Grid,Orbit,l,energy,Pot%rv,Pot%v0,Pot%v0p, &
&                    PAW%phi(:,nbase),PAW%rng(nbase))
         else
            CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,&
&                 nr,l,energy,PAW%phi(:,nbase),nodes)
         endif
         rat=MAX(ABS(PAW%phi(irc,nbase)),ABS(PAW%phi(irc+1,nbase)))
         rat=DSIGN(rat,PAW%phi(irc,nbase))
         PAW%phi(1:n,nbase)=PAW%phi(1:n,nbase)/rat
         !IF(Orbit%exctype=='HF') PAW%lmbd(:,nbase)=PAW%lmbd(:,nbase)/rat
         WRITE(STD_OUT,'(3i6,1p,2e15.6)') nbase,PAW%np(nbase),l,             &
&               PAW%eig(nbase),PAW%occ(nbase)
         nbl=nbl+1
         ibasis_add=ibasis_add+1
       ENDDO generalizedloop
       !
    ENDDO   ! end lmax loop

    WRITE(std_out,*) 'completed phi basis with ',nbase,' functions '
    PAW%nbase=nbase     ! reset nbase


!  Check density
!    checkden=0
!    do io=1,PAW%nbase
!       checkden=checkden+PAW%occ(io)*(PAW%phi(:,io)**2)
!    enddo   

!   open(1001,file='Checkdeninsetbasis',form='formatted')
!    do i=1,n
!       write(1001,'(1p, 50e16.7)') Grid%r(i),checkden(i),valden(i)
!    enddo   
!   close(1001)
!   deallocate(checkden,valden)
  END SUBROUTINE setbasis

  !**************************************************************************
  !  Program to generate atomic basis functions
  !    Version using Bloechl's form of projector and orthogonalization procedure
  !     At the end of this subroutine, the basis functions and projectors are
  !     orthogonalized with a Gram-Schmidt like scheme
  !**************************************************************************
  SUBROUTINE makebasis_bloechl(Grid,Pot,option)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: POT
    INTEGER,INTENT(IN) :: option

    INTEGER :: n,irc,nr
    INTEGER :: i,j,k,io,ok,lmax
    REAL(8) :: h,rc
    REAL(8), ALLOCATABLE :: tmp(:),VNC(:)
    TYPE(PotentialInfo), TARGET:: PS
    REAL(8), POINTER  :: r(:)

    if(needvtau) then
       write(std_out,*) 'Sorry -- this case has not been programmed yet'
       write(std_out,*)  ' need to modify bsolv routine for extra derivatives '
       stop
    endif   
    n=Grid%n
    h=Grid%h
    r=>Grid%r
    irc=PAW%irc
    nr=irc+20
    rc=PAW%rc
    lmax=PAW%lmax
    PS%nz=0.d0

    ALLOCATE(tmp(n),VNC(n),PS%rv(n),stat=ok)
    IF (ok /= 0 ) THEN
       WRITE(STD_OUT,*) 'Error in allocating  denout in makebasis ',n,ok
       STOP
    ENDIF

    ! find basis functions
    PS%rv=PAW%rveff
    call zeropot(Grid,PS%rv,PS%v0,PS%v0p)
    !write(std_out,*) 'VNC  v0 ', PS%v0,PS%v0p,PS%rv(5)
    CALL formprojectors(Grid,Pot,PS,option)

    DEALLOCATE(tmp,VNC,PS%rv)
    do io=1,PAW%nbase
      PAW%rcio(io)=PAW%rc
    end do

  END SUBROUTINE makebasis_bloechl


  !**************************************************************************
  !  Program to generate atomic basis functions
  !
  !   1) Pseudization of partial waves:
  !        - simple polynom scheme                                            [optps=1]
  !                   r^(l+1).Sum[Ci.r^2i]  0<=i<=4
  !   OR   - ultrasoft polynom scheme                                         [optps=2]
  !                   r^(l+1).{Sum[Ci.r^2i]+Sum[Cj.r^2j]}  0<=i<=3
  !                           3<j adjusted using Fourier filtering
  !   OR   - RRKJ scheme with 2 Bessel functions (PHYS REV B 41,1227 (1990))  [optps=3]
  !
  !   2) Build and orthogonalization of projectors
  !        - Vanderbilt generation method (PHYS REV B 41,7892 (1990))  [optorth=0]
  !   OR   - Gram-Schmidt like sheme                                   [optorth=1]
  !**************************************************************************
  SUBROUTINE makebasis_custom(Grid,Pot,optps,optorth,pdeg,qcut)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    INTEGER, INTENT(IN) :: optps,optorth,pdeg
    REAL(8), INTENT(IN) :: qcut

    INTEGER :: i,j,k,l,io,jo,ok,lmax,nbase,n,irc,irc_vloc,nr,np,thisrc
    INTEGER :: icount,jcount,istart,ifinish,ibase,jbase
    REAL(8) :: choice,rc,xx,yy,gg,g,gp,gpp,al(2),ql(2)
    INTEGER, ALLOCATABLE :: omap(:)
    REAL(8), ALLOCATABLE :: VNC(:),Ci(:),aa(:,:),ai(:,:)
    REAL(8), POINTER  :: r(:)
 
    if (optps<1.or.optps>3) stop 'bug: error calling makebasis_custom routine'
    if (optorth<0.or.optorth>1) stop 'bug: error calling makebasis_custom routine'

    n=Grid%n
    r=>Grid%r
    irc=PAW%irc
    irc_vloc=PAW%irc_vloc
    nbase=PAW%nbase
    lmax=PAW%lmax

    np=5;if (optps==2) np=pdeg+1 ;
    if (optps==1.or.optps==2) allocate(Ci(np))

  ! Set screened local pseudopotential
    allocate(VNC(n),stat=i)
    if (i/=0) stop 'allocation error in makebasis_vanderbilt'
    VNC(2:n)=PAW%rveff(2:n)/r(2:n)
    call extrapolate(Grid,VNC)

  ! Loop on basis elements
    do io=1,nbase
       l=PAW%l(io)

     ! Read matching radius (from input dataset)
       thisrc=FindGridIndex(Grid,input_dataset%basis_func_rc(io))
       thisrc=MIN(thisrc,irc)       ! make sure rc<total rc
       rc=r(thisrc);PAW%rcio(io)=rc
       WRITE(STD_OUT,'(a,3i5,1p,1e15.7)') ' For this wfn: ',io,PAW%np(io),PAW%l(io),PAW%eig(io)
       WRITE(STD_OUT,'(a,f10.7)') '  >>> rc =', rc
       if (thisrc<3.or.thisrc>irc.or. &
&          (optps==1.and.thisrc>n-3).or. &
&          (optps==2.and.thisrc>n-6)) then
          write(std_out,*) 'rc out of range', thisrc,n,irc
          stop
       endif

     ! Find partial wave pseudization
       if (optps==1) then
        call pspolyn(PAW%phi(:,io),Ci,r,l,np,thisrc,n)
       else if (optps==2) then
        call psuspolyn(PAW%phi(:,io),Ci,r,l,np,thisrc,n,qcut)
       else if (optps==3) then
        call psbes(PAW%phi(:,io),al,ql,Grid,l,thisrc,n)
       endif

     ! Compute pseudized partial wave and unnormalized projector
       PAW%tphi(:,io)=PAW%phi(:,io)
       PAW%tp(:,io)=0.d0
       if (optps==1.or.optps==2) then
        do i=1,thisrc-1
         xx=r(i)*r(i)
         PAW%tphi(i,io)=Ci(1)+Ci(2)*xx
         gp=2.d0*Ci(2)
         gpp=2.d0*Ci(2)
         do j=3,np
          PAW%tphi(i,io)=PAW%tphi(i,io)+Ci(j)*xx**(j-1)
          gp=gp+dble(2*j-2)*Ci( j)*xx**(j-2)
          gpp=gpp+dble((2*j-2)*(2*j-3))*Ci(j)*xx**(j-2)
         enddo
         PAW%tphi(i,io)=PAW%tphi(i,io)*r(i)**(l+1)
         PAW%tp(i,io)=(dble(2*(l+1))*gp+gpp)*r(i)**(l+1)+&
&                   (PAW%eig(io)-VNC(i))*PAW%tphi(i,io)
        enddo
       else if (optps==3) then
        PAW%tphi(1,io)=0.d0
        do i=2,thisrc-1
         xx=ql(1)*r(i)
         call jbessel(g,gp,gpp,l,2,xx)
         PAW%tphi(i,io)=al(1)*g*r(i)
         gg=al(1)*(2.d0*ql(1)*gp+ql(1)*xx*gpp)
         xx=ql(2)*r(i)
         call jbessel(g,gp,gpp,l,2,xx)
         PAW%tphi(i,io)=PAW%tphi(i,io)+al(2)*g*r(i)
         gg=gg+al(2)*(2.d0*ql(2)*gp+ql(2)*xx*gpp)
         PAW%tp(i,io)=(PAW%eig(io)-VNC(i)-dble(l*(l+1))/(r(i)**2))*PAW%tphi(i,io)+gg
        enddo
       endif

       if (thisrc<irc_vloc) then
         do i=thisrc,irc_vloc-1
             gpp=Gsecondderiv(Grid,i,PAW%tphi(:,io))
             PAW%tp(i,io)=(PAW%eig(io)-VNC(i)-dble(l*(l+1))/(r(i)**2))*PAW%tphi(i,io)+gpp
         enddo
       endif

     ! If Gram-Schmidt orthogonalization, form normalized projector
       if (optorth==1) then
        xx=overlap(Grid,PAW%tp(:,io),PAW%tphi(:,io),1,max(thisrc,irc_vloc))
        PAW%tp(:,io)=PAW%tp(:,io)/xx
       endif

    enddo  !nbase

    deallocate(VNC);if (optps==1.or.optps==2) deallocate(Ci)

  ! Form orthogonalized projector functions
    do io=1,nbase
       PAW%ophi(:,io)=PAW%phi(:,io)
       PAW%otphi(:,io)=PAW%tphi(:,io)
       PAW%Kop(1,io)=0
       PAW%Kop(2:n,io)=(PAW%eig(io)-Pot%rv(2:n)/Grid%r(2:n))*PAW%phi(2:n,io)
       if (optorth==1) PAW%otp(:,io)=PAW%tp(:,io)
    enddo

  ! First option : VANDERBILT SCHEME
    if (optorth==0) then
     do l=0,lmax
       icount=0
       do io=1,nbase
        if (PAW%l(io)==l) icount=icount+1
       enddo
       write(std_out,*) 'For l = ', l,icount,' basis functions'
       if (icount==0) cycle
       allocate(aa(icount,icount),ai(icount,icount),omap(icount))
       aa=0;icount=0
       do io=1,nbase
        if (PAW%l(io)==l) then
          icount=icount+1
          omap(icount)=io
        endif
       enddo
       do i=1,icount
         io=omap(i)
         do j=1,icount
           jo=omap(j)
           aa(i,j)=overlap(Grid,PAW%otphi(:,io),PAW%tp(:,jo),1,irc)
         enddo
       enddo
       ai=aa;call minverse(ai,icount,icount,icount)

       do i=1,icount
         io=omap(i)
         PAW%ck(io)=ai(i,i)
         PAW%otp(:,io)=0
         do j=1,icount
           jo=omap(j)
           PAW%otp(:,io)=PAW%otp(:,io)+PAW%tp(:,jo)*ai(j,i)
         enddo
       enddo
       deallocate(aa,ai,omap)
     enddo

  ! Second option : GRAM-SCHMIDT SCHEME
    else if (optorth==1) then
     DO l=0,lmax
       icount=0
       do io=1,nbase
        if (PAW%l(io)==l) icount=icount+1
       enddo
       if (icount==0) cycle
       allocate(aa(icount,icount),ai(icount,icount),omap(icount))
       icount=0
       DO io=1,nbase
         IF (PAW%l(io)==l) THEN
           IF  (icount==0) istart=io
           IF  (icount>=0) ifinish=io
           icount=icount+1;omap(icount)=io
         ENDIF
       ENDDO
       DO ibase=istart,ifinish
         DO jbase=istart,ibase
           IF (jbase.LT.ibase) THEN
             xx=overlap(Grid,PAW%otp(:,jbase),PAW%otphi(:,ibase),1,irc)
             yy=overlap(Grid,PAW%otphi(:,jbase),PAW%otp(:,ibase),1,irc)
             PAW%ophi(1:n,ibase)=PAW%ophi(1:n,ibase)-PAW%ophi(1:n,jbase)*xx
             PAW%Kop(1:n,ibase)=PAW%Kop(1:n,ibase)-PAW%Kop(1:n,jbase)*xx
             PAW%otphi(1:n,ibase)=PAW%otphi(1:n,ibase)-PAW%otphi(1:n,jbase)*xx
             PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)-PAW%otp(1:n,jbase)*yy
             aa(ibase-istart+1,jbase-istart+1)=xx
             aa(jbase-istart+1,ibase-istart+1)=xx
           ELSE IF (jbase.EQ.ibase) THEN
             xx=overlap(Grid,PAW%otp(:,jbase),PAW%otphi(:,ibase),1,irc)
             choice=1.d0/SQRT(ABS(xx))
             PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)*DSIGN(choice,xx)
             PAW%otphi(1:n,ibase)=PAW%otphi(1:n,ibase)*choice
             PAW%ophi(1:n,ibase)=PAW%ophi(1:n,ibase)*choice
             PAW%Kop(1:n,ibase)=PAW%Kop(1:n,ibase)*choice
             aa(ibase-istart+1,ibase-istart+1)=xx
           ENDIF
         ENDDO
       ENDDO
       ai=aa;  !call minverse(aa,icount,icount,icount)
       do i=1,icount
        io=omap(i);PAW%ck(io)=ai(i,i)
       enddo
      deallocate(aa,ai,omap)
     ENDDO
    endif

  END SUBROUTINE makebasis_custom

  !**************************************************************************
  !  Program to generate atomic basis functions -- MGGA version
  !     Part 1 -- generate pseudo basis functions only
  !
  !   1) Pseudization of partial waves:
  !        - simple polynom scheme                                            [optps=1]
  !                   r^(l+1).Sum[Ci.r^2i]  0<=i<=4
  !   OR   - ultrasoft polynom scheme                                         [optps=2]
  !                   r^(l+1).{Sum[Ci.r^2i]+Sum[Cj.r^2j]}  0<=i<=3
  !                           3<j adjusted using Fourier filtering
  !   OR   - RRKJ scheme with 2 Bessel functions (PHYS REV B 41,1227 (1990))  [optps=3]
  !
  !   2) Build and orthogonalization of projectors
  !        - Vanderbilt generation method (PHYS REV B 41,7892 (1990))  [optorth=0]
  !   OR   - Gram-Schmidt like sheme                                   [optorth=1]
  !   For this case, only Vanderbilt scheme programmed
  !**************************************************************************
  SUBROUTINE makebasis_custom_vtau1(Grid,Pot,optps,optorth,pdeg,qcut)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    INTEGER, INTENT(IN) :: optps,optorth,pdeg
    REAL(8), INTENT(IN) :: qcut

    INTEGER :: i,j,k,l,io,jo,ok,lmax,nbase,n,irc,irc_vloc,nr,np,thisrc
    INTEGER :: icount,jcount,istart,ifinish,ibase,jbase
    REAL(8) :: choice,rc,xx,yy,gg,g,gp,gpp,al(2),ql(2)
    REAL(8), ALLOCATABLE :: Ci(:)
    REAL(8), POINTER  :: r(:)
    
    write(std_out,*) 'Entering makebasis_custom_vtau1'
    if (optps<1.or.optps>3) stop 'bug: error calling makebasis_custom routine'
    if (optorth<0.or.optorth>1) stop 'bug: error calling makebasis_custom routine'

    n=Grid%n
    r=>Grid%r
    irc=PAW%irc
    nbase=PAW%nbase
    lmax=PAW%lmax

    np=5;if (optps==2) np=pdeg+1 ;
    if (optps==1.or.optps==2) allocate(Ci(np))

  ! Loop on basis elements
    do io=1,nbase
       l=PAW%l(io)

     ! Read matching radius
       thisrc=FindGridIndex(Grid,input_dataset%basis_func_rc(io))
       thisrc=MIN(thisrc,irc)       ! make sure rc<total rc
       rc=r(thisrc);PAW%rcio(io)=rc
       WRITE(std_out,'(a,3i5,1p,1e15.7)') ' For this wfn: ',io,PAW%np(io),PAW%l(io),PAW%eig(io)
       WRITE(std_out,'(a,f10.7)') '  >>> rc =', rc
       if (thisrc<3.or.thisrc>irc.or. &
&          (optps==1.and.thisrc>n-3).or. &
&          (optps==2.and.thisrc>n-6)) then
          write(std_out,*) 'rc out of range', thisrc,n,irc
          stop
       endif

     ! Find partial wave pseudization
       if (optps==1) then
        call pspolyn(PAW%phi(:,io),Ci,r,l,np,thisrc,n)
       else if (optps==2) then
        call psuspolyn(PAW%phi(:,io),Ci,r,l,np,thisrc,n,qcut)
       else if (optps==3) then
        call psbes(PAW%phi(:,io),al,ql,Grid,l,thisrc,n)
       endif

     ! Compute pseudized partial wave and unnormalized projector
       PAW%tphi(:,io)=PAW%phi(:,io)
       if (optps==1.or.optps==2) then
        do i=1,thisrc-1
         xx=r(i)*r(i)
         PAW%tphi(i,io)=Ci(1)+Ci(2)*xx
         gp=2.d0*Ci(2)
         gpp=2.d0*Ci(2)
         do j=3,np
          PAW%tphi(i,io)=PAW%tphi(i,io)+Ci(j)*xx**(j-1)
          gp=gp+dble(2*j-2)*Ci( j)*xx**(j-2)
          gpp=gpp+dble((2*j-2)*(2*j-3))*Ci(j)*xx**(j-2)
         enddo
         PAW%tphi(i,io)=PAW%tphi(i,io)*r(i)**(l+1)
        enddo
       else if (optps==3) then
        PAW%tphi(1,io)=0.d0
        do i=2,thisrc-1
         xx=ql(1)*r(i)
         call jbessel(g,gp,gpp,l,2,xx)
         PAW%tphi(i,io)=al(1)*g*r(i)
         gg=al(1)*(2.d0*ql(1)*gp+ql(1)*xx*gpp)
         xx=ql(2)*r(i)
         call jbessel(g,gp,gpp,l,2,xx)
         PAW%tphi(i,io)=PAW%tphi(i,io)+al(2)*g*r(i)
         gg=gg+al(2)*(2.d0*ql(2)*gp+ql(2)*xx*gpp)
        enddo
       endif

    enddo  !nbase

    if (optps==1.or.optps==2) deallocate(Ci)

  END SUBROUTINE makebasis_custom_vtau1
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! makebasis_modrrkj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE makebasis_modrrkj(Grid,Pot,Orthoindex,success)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    INTEGER, INTENT(IN) :: Orthoindex

    INTEGER :: i,j,k,l,io,jo,ok,lmax,nbase,n,irc,irc_vloc,nr,np,thisrc
    INTEGER :: icount,jcount,istart,ifinish,ibase,jbase,lprev
    INTEGER :: match=5
    REAL(8) :: dk
    REAL(8) :: choice,rc,xx,yy,gg,g,gp,gpp,al(2),ql(2),root1,root2,logderiv
    INTEGER, ALLOCATABLE :: omap(:),rcindex(:)
    REAL(8), ALLOCATABLE :: VNC(:),Ci(:),aa(:,:),ai(:,:),jl(:,:),kv(:)
    REAL(8), ALLOCATABLE :: f(:),rcval(:)
    REAL(8), ALLOCATABLE :: U(:,:),VT(:,:),WORK(:),S(:),X(:,:)
    INTEGER :: LWORK
    REAL(8), POINTER  :: r(:)
    LOGICAL :: success

    success=.false.
    n=Grid%n
    r=>Grid%r
    irc=PAW%irc
    irc_vloc=PAW%irc_vloc
    nbase=PAW%nbase
    lmax=PAW%lmax

    ALLOCATE(rcindex(nbase),rcval(nbase))

  !  Input matching radii for each basis function

      call readmatchradius(Grid,rcindex,rcval)


  ! Set screened local pseudopotential
    allocate(VNC(n),stat=i)
    if (i/=0) stop 'allocation error in make_modrrkj'
    VNC(2:n)=PAW%rveff(2:n)/r(2:n)
    call extrapolate(Grid,VNC)


    allocate(kv(match),jl(n,match),f(n),Ci(match),aa(match,match))
    if (i/=0) stop 'allocation error in make_nrecipe'
  ! Loop on basis elements
    lprev=-1;np=0
    do io=1,nbase

       PAW%Kop(1,io)=0
       PAW%Kop(2:n,io)=(PAW%eig(io)-Pot%rv(2:n)/Grid%r(2:n))*PAW%phi(2:n,io)

       l=PAW%l(io)
       if (l==lprev) then
         np=np+1
       else
         np=1
         lprev=l
       endif

       thisrc=rcindex(io)
       rc=rcval(io);PAW%rcio(io)=rc

    ! Set normalization for PAW%eig(io)>0
      if (PAW%eig(io)>0) then
         PAW%phi(:,io)=PAW%phi(:,io)/PAW%phi(thisrc,io)
      endif
     ! Find logderiv
       logderiv=Gfirstderiv(Grid,thisrc,PAW%phi(:,io))/PAW%phi(thisrc,io)
       write(std_out,*) 'logderiv ', io, l, logderiv

     !  Find   bessel function zeros
        call solvbes(kv,1.d0,0.d0,l,np)
        write(std_out,*) 'Searching for solution ',kv(np)
        if (np==1) then
         root1=0.001d0; root2=kv(1)-0.001d0; xx=0.5d0*(root1+root2)
        else
         root1=kv(np-1)+0.001d0; root2=kv(np)-0.001d0; xx=0.5d0*(root1+root2)
        endif
        icount=0
        Do
          call jbessel(g,gp,gpp,l,2,xx)
          yy=(logderiv-(1.d0/xx+gp/g))/(-1.d0/xx**2+gpp/g-(gp/g)**2)
          if (abs(yy)<1.d-6) then
            write(std_out,*) 'exiting Bessel loop', icount,xx,yy
            exit
          else
            if (xx+yy>root2) then
                 xx=0.5*(xx+root2)
            else if (xx+yy<root1) then
                 xx=0.5*(xx+root1)
            else
                xx=xx+yy
            endif
          endif
          !write(std_out,*) 'iter Bessel', xx,yy
          icount=icount+1
          if (icount>1000) then
            write(std_out,*) 'Giving up on Bessel root'
            return
          endif
        Enddo

        success=.true.
        write(std_out,*) 'Found Bessel root at xx' , xx
        dk=(pi/40)/rc             !  could be made adjustable
        write(std_out,*) 'dk value for this case ', dk
        kv=0
        kv(1)=xx/rc-dk*(match-1)/2
        write(std_out,*) '#  kv      1', kv(1)
        Do i=2,match
          kv(i)=kv(i-1)+dk
          write(std_out,*) '#  kv   ',   i, kv(i)
        enddo

        Do i=1,match
          f(:)=kv(i)*r(:)
          call sphbes(l,n,f)
          jl(:,i)=f(:)*r(:)
        Enddo

      ! Match bessel functions with wavefunctions

        Ci=0; aa=0
        Do j=1,match
           i=match/2-j+1
           Ci(j)=PAW%phi(thisrc-i*1,io)                ! 1 step should vary
           aa(j,1:match)=jl(thisrc-i*1,1:match)
        enddo

        call SolveAXeqBM(match,aa,Ci,match-1)
        write(std_out,*) 'Completed SolveAXeqB with coefficients'
        write(std_out,'(1p,10e15.7)') (Ci(i),i=1,match)

     ! Compute pseudized partial wave and unnormalized projector
       PAW%tphi(:,io)=PAW%phi(:,io);
       PAW%tp(:,io)=0.d0
        do i=1,thisrc-1
         PAW%tphi(i,io)=sum(Ci(1:match)*jl(i,1:match))
         do j=1,match
            PAW%tp(i,io)=PAW%tp(i,io)+&
&             Ci(j)*(-PAW%eig(io)+(kv(j)**2)+VNC(i))*jl(i,j)
         enddo
        enddo

       if (thisrc<=irc_vloc) then
         do i=thisrc,irc_vloc-1
             gpp=Gsecondderiv(Grid,i,PAW%tphi(:,io))
             PAW%tp(i,io)=(-PAW%eig(io)&
&              +VNC(i)+dble(l*(l+1))/(r(i)**2))*PAW%tphi(i,io)  -gpp
         enddo
       endif


    enddo  !nbase

     Deallocate(aa)
!
 Selectcase(Orthoindex)
   Case default     ! Orthoindex=VANDERBILTORTHO
  !! Form orthogonalized projectors --   VANDERBILTORTHO
     write(std_out,*) ' Vanderbilt ortho'
     do l=0,lmax
       icount=0
       do io=1,nbase
        if (PAW%l(io)==l) icount=icount+1
       enddo
       if (icount==0) cycle
       write(std_out,*) 'For l = ', l,icount,' basis functions'
       allocate(aa(icount,icount),ai(icount,icount),omap(icount))
       aa=0;icount=0
       do io=1,nbase
        if (PAW%l(io)==l) then
          icount=icount+1
          omap(icount)=io
        endif
       enddo
       do i=1,icount
         io=omap(i)
         PAW%otphi(:,io)=PAW%tphi(:,io)
         PAW%ophi(:,io)=PAW%phi(:,io)
       enddo
       do i=1,icount
          io=omap(i)
         do j=1,icount
           jo=omap(j)
           aa(i,j)=overlap(Grid,PAW%otphi(:,io),PAW%tp(:,jo),1,irc)
         enddo
       enddo
       ai=aa;call minverse(ai,icount,icount,icount)

       do i=1,icount
         io=omap(i)
         PAW%ck(io)=ai(i,i)
         PAW%otp(:,io)=0
         do j=1,icount
           jo=omap(j)
           PAW%otp(:,io)=PAW%otp(:,io)+PAW%tp(:,jo)*ai(j,i)
         enddo
       enddo

       write(std_out,*) 'Check  otp for l = ', l
       do i = 1, icount
          io=omap(i)
          do j = 1, icount
             jo=omap(j)
             write(std_out,*) 'Overlap i j ', i,j, &
&                    overlap(Grid,PAW%otphi(:,io),PAW%otp(:,jo),1,irc)
          enddo
       enddo
       deallocate(aa,ai,omap)
    enddo

  Case(GRAMSCHMIDTORTHO) ! Form orthonormal projectors --  GRAM-SCHMIDT SCHEME
     write(std_out,*) 'Gramschmidt ortho'
     DO l=0,lmax
        icount=0
        do io=1,nbase
           if (PAW%l(io)==l) icount=icount+1
        enddo
        if (icount==0) cycle
        allocate(aa(icount,icount),ai(icount,icount),omap(icount))
        icount=0
        DO io=1,nbase
          IF (PAW%l(io)==l) THEN
             IF (icount==0) istart=io
             IF (icount>=0) ifinish=io
               icount=icount+1;omap(icount)=io
          ENDIF
       ENDDO
     DO ibase=istart,ifinish
          PAW%otp(:,ibase)=PAW%tp(:,ibase)
          PAW%otphi(:,ibase)=PAW%tphi(:,ibase)
          PAW%ophi(:,ibase)=PAW%phi(:,ibase)
          DO jbase=istart,ibase
             IF (jbase.LT.ibase) THEN
             xx=overlap(Grid,PAW%otp(:,jbase),PAW%otphi(:,ibase),1,irc)
             yy=overlap(Grid,PAW%otphi(:,jbase),PAW%otp(:,ibase),1,irc)
             PAW%ophi(1:n,ibase)=PAW%ophi(1:n,ibase)-PAW%ophi(1:n,jbase)*xx
             PAW%Kop(1:n,ibase)=PAW%Kop(1:n,ibase)-PAW%Kop(1:n,jbase)*xx
             PAW%otphi(1:n,ibase)=PAW%otphi(1:n,ibase)-PAW%otphi(1:n,jbase)*xx
             PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)-PAW%otp(1:n,jbase)*yy
             aa(ibase-istart+1,jbase-istart+1)=xx
             aa(jbase-istart+1,ibase-istart+1)=xx
             ELSE IF (jbase.EQ.ibase) THEN
                   xx=overlap(Grid,PAW%otp(:,jbase),PAW%otphi(:,ibase),1,irc)
                   !choice=1.d0/SQRT(ABS(xx))
                   !PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)*DSIGN(choice,xx)
                   !PAW%otphi(1:n,ibase)=PAW%otphi(1:n,ibase)*choice
                   !PAW%ophi(1:n,ibase)=PAW%ophi(1:n,ibase)*choice
                   !PAW%Kop(1:n,ibase)=PAW%Kop(1:n,ibase)*choice
                   PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)/xx
                   aa(ibase-istart+1,ibase-istart+1)=xx
             ENDIF
          ENDDO
       ENDDO
       ai=aa;
       do i=1,icount
        io=omap(i);PAW%ck(io)=ai(i,i)
       enddo
       deallocate(aa,ai,omap)
     ENDDO

  Case (SVDORTHO)    ! SVD construction
     write(std_out,*) 'SVD ortho'
     DO l=0,lmax
        icount=0
        do io=1,nbase
           if (PAW%l(io)==l) icount=icount+1
        enddo
        if (icount==0) cycle
        allocate(aa(icount,icount),ai(icount,icount),omap(icount),X(n,icount))
        icount=0
        DO io=1,nbase
          IF (PAW%l(io)==l) THEN
             IF (icount==0) istart=io
             IF (icount>=0) ifinish=io
               icount=icount+1;omap(icount)=io
          ENDIF
        ENDDO
        do i=1,icount
          io=omap(i)
          PAW%otphi(:,io)=0
          PAW%ophi(:,io)=0
          PAW%otp(:,io)=0
          do j=1,icount
            jo=omap(j)
            aa(i,j)=overlap(Grid,PAW%tphi(:,io),PAW%tp(:,jo),1,irc)
          enddo
        enddo
        LWORK=MAX(200,icount,icount)
        ALLOCATE(U(icount,icount),VT(icount,icount),WORK(LWORK),S(icount))
        ai=aa;
        CALL DGESVD('A','A',icount,icount,ai,icount,S,&
&          U,icount,VT,icount,WORK,LWORK,k)

        write(std_out,*) 'For l = ', l , 'Completed SVD'
        write(std_out,*) 'S = ', S(:)

        Do i=1,icount
           io=omap(i)
           X(:,i)=PAW%Kop(:,io);
        Enddo
        do j=1,icount
           jo=omap(j)
           PAW%ophi(:,jo)=0
           PAW%otphi(:,jo)=0
           PAW%otp(:,jo)=0
           PAW%Kop(:,jo)=0
           do i=1,icount
              io=omap(i)
              PAW%ophi(:,jo)=PAW%ophi(:,jo)+PAW%phi(:,io)*U(i,j)
              PAW%otphi(:,jo)=PAW%otphi(:,jo)+PAW%tphi(:,io)*U(i,j)
              PAW%Kop(:,jo)=PAW%Kop(:,jo)+X(:,i)*U(i,j)
              PAW%otp(:,jo)=PAW%otp(:,jo)+PAW%tp(:,io)*VT(j,i)/S(j)
           enddo
        enddo

         do i=1,icount
            io=omap(i)
            do j=1,icount
               jo=omap(j)
               write(std_out,*) 'Overlap ', i,j, &
&                 overlap(Grid,PAW%otphi(:,io),PAW%otp(:,jo),1,irc)
            enddo
         enddo

      deallocate(aa,ai,omap,U,VT,S,WORK,X)

    ENDDO

  End Select

  Deallocate(VNC,kv,jl,f,Ci,rcindex,rcval)
  END SUBROUTINE makebasis_modrrkj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! makebasis_modrrkj_vtau1
!     version for MGGA case -- Part 1 to determine pseudobasis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE makebasis_modrrkj_vtau1(Grid,Pot,Orthoindex,success)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PotentialInfo), INTENT(IN) :: Pot
    INTEGER, INTENT(IN) :: Orthoindex

    INTEGER :: i,j,k,l,io,jo,ok,lmax,nbase,n,irc,irc_vloc,nr,np,thisrc
    INTEGER :: icount,jcount,istart,ifinish,ibase,jbase,lprev
    INTEGER :: match=5
    REAL(8) :: dk
    REAL(8) :: choice,rc,xx,yy,gg,g,gp,gpp,al(2),ql(2),root1,root2,logderiv
    INTEGER, ALLOCATABLE :: omap(:),rcindex(:)
    REAL(8), ALLOCATABLE :: VNC(:),Ci(:),aa(:,:),ai(:,:),jl(:,:),kv(:)
    REAL(8), ALLOCATABLE :: f(:),rcval(:)
    REAL(8), ALLOCATABLE :: U(:,:),VT(:,:),WORK(:),S(:),X(:,:)
    INTEGER :: LWORK
    REAL(8), POINTER  :: r(:)
    LOGICAL :: success

    write(std_out,*) 'Entering makebasis_modrrkj_vtau1'
    success=.false.
    n=Grid%n
    r=>Grid%r
    irc=PAW%irc
    nbase=PAW%nbase
    lmax=PAW%lmax

    ALLOCATE(rcindex(nbase),rcval(nbase))

   !Input matching radii for each basis function
    call readmatchradius(Grid,rcindex,rcval)

    allocate(kv(match),jl(n,match),f(n),Ci(match),aa(match,match))

  ! Loop on basis elements
    lprev=-1;np=0
    do io=1,nbase

       PAW%Kop(1,io)=0
       PAW%Kop(2:n,io)=(PAW%eig(io)-Pot%rv(2:n)/Grid%r(2:n))*PAW%phi(2:n,io)

       l=PAW%l(io)
       if (l==lprev) then
         np=np+1
       else
         np=1
         lprev=l
       endif

       thisrc=rcindex(io)
       rc=rcval(io);PAW%rcio(io)=rc

    ! Set normalization for PAW%eig(io)>0
      if (PAW%eig(io)>0) then
         PAW%phi(:,io)=PAW%phi(:,io)/PAW%phi(thisrc,io)
      endif
     ! Find logderiv
       logderiv=Gfirstderiv(Grid,thisrc,PAW%phi(:,io))/PAW%phi(thisrc,io)
       write(std_out,*) 'logderiv ', io, l, logderiv

     !  Find   bessel function zeros
        call solvbes(kv,1.d0,0.d0,l,np)
        write(std_out,*) 'Searching for solution ',kv(np)
        if (np==1) then
         root1=0.001d0; root2=kv(1)-0.001d0; xx=0.5d0*(root1+root2)
        else
         root1=kv(np-1)+0.001d0; root2=kv(np)-0.001d0; xx=0.5d0*(root1+root2)
        endif
        icount=0
        Do
          call jbessel(g,gp,gpp,l,2,xx)
          yy=(logderiv-(1.d0/xx+gp/g))/(-1.d0/xx**2+gpp/g-(gp/g)**2)
          if (abs(yy)<1.d-6) then
            write(std_out,*) 'exiting Bessel loop', icount,xx,yy
            exit
          else
            if (xx+yy>root2) then
                 xx=0.5*(xx+root2)
            else if (xx+yy<root1) then
                 xx=0.5*(xx+root1)
            else
                xx=xx+yy
            endif
          endif
          !write(std_out,*) 'iter Bessel', xx,yy
          icount=icount+1
          if (icount>1000) then
            write(std_out,*) 'Giving up on Bessel root'
            return
          endif
        Enddo

        success=.true.
        write(std_out,*) 'Found Bessel root at xx' , xx
        dk=(pi/40)/rc             !  could be made adjustable
        write(std_out,*) 'dk value for this case ', dk
        kv=0
        kv(1)=xx/rc-dk*(match-1)/2
        write(std_out,*) '#  kv      1', kv(1)
        Do i=2,match
          kv(i)=kv(i-1)+dk
          write(std_out,*) '#  kv   ',   i, kv(i)
        enddo

        Do i=1,match
          f(:)=kv(i)*r(:)
          call sphbes(l,n,f)
          jl(:,i)=f(:)*r(:)
        Enddo

      ! Match bessel functions with wavefunctions

        Ci=0; aa=0
        Do j=1,match
           i=match-j
           Ci(j)=PAW%phi(thisrc-i,io)                ! 1 step should vary
           aa(j,1:match)=jl(thisrc-i,1:match)
        enddo

        call SolveAXeqBM(match,aa,Ci,match-1)
        write(std_out,*) 'Completed SolveAXeqB with coefficients'
        write(std_out,'(1p,10e15.7)') (Ci(i),i=1,match)

     ! Compute pseudized partial wave 
       PAW%tphi(:,io)=PAW%phi(:,io);
        do i=1,thisrc-1
         PAW%tphi(i,io)=sum(Ci(1:match)*jl(i,1:match))
        enddo

    enddo  !nbase

    deallocate(kv,jl,f,Ci,aa)
  END SUBROUTINE makebasis_modrrkj_vtau1

  SUBROUTINE readmatchradius(Grid,rcindex,rcval)
      TYPE(GridInfo), INTENT(IN) :: Grid
      INTEGER, INTENT(INOUT) :: rcindex(:)
      REAL(8), INTENT(INOUT) :: rcval(:)

      INTEGER :: io,nbase,n,irc,lmax,thisrc
      REAL(8) :: rc

    nbase=PAW%nbase
    irc=PAW%irc  

  ! Loop on basis elements
    rcindex=0;rcval=0
    do io=1,nbase

      ! Read matching radius (from input dataset)
      thisrc=FindGridIndex(Grid,input_dataset%basis_func_rc(io))
      thisrc=MIN(thisrc,irc)       ! make sure rc<total rc
      rc=Grid%r(thisrc)
      WRITE(STD_OUT,'(a,3i5,1p,1e15.7)') ' For this wfn: ',io,PAW%np(io),PAW%l(io),PAW%eig(io)
      WRITE(STD_OUT,'(a,f10.7)') '  >>> rc =', rc
      if (thisrc<3.or.thisrc>irc) then
        write(std_out,*) 'rc out of range', thisrc,n,irc
        stop
      endif
      rcindex(io)=thisrc;rcval(io)=rc
      write(std_out,'(" For io = ", i5," rc = ", i5,f10.5)') io,rcindex(io),rcval(io)

   ENDDO
  END SUBROUTINE      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  calculate_tvtau    for MGGA case
!     Assume valence pseudo wavefunctions known
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calculate_tvtau(Grid)
    TYPE(GridInfo), INTENT(IN) :: Grid

    INTEGER :: i,j,k,l
    REAL(8), allocatable :: dp(:),ddp(:),vxc(:),tvxc(:),locald(:),localtd(:)
    REAL(8), allocatable :: Ktvxc(:),Kd(:)
    REAL(8) :: fac,exc,texc,sum,tsum
    REAL(8), parameter :: small=1.d-5

    allocate(dp(Grid%n),ddp(Grid%n),vxc(Grid%n),tvxc(Grid%n))
    allocate(locald(Grid%n),localtd(Grid%n),Ktvxc(Grid%n),Kd(Grid%n))

    locald=PAW%core
    localtd=PAW%tcore
    PAW%valetau=0.d0
    PAW%tvaletau=0.d0

    do i=1,PAW%nbase
      write(std_out,*) 'tvtau -- ', i,PAW%l(i),PAW%occ(i),PAW%eig(i)
      if (PAW%occ(i).gt.small) then
         locald=locald+PAW%occ(i)*(PAW%phi(:,i)**2)
         localtd=localtd+PAW%occ(i)*(PAW%tphi(:,i)**2)
         CALL taufromwfn(dp,Grid,PAW%phi(:,i),PAW%l(i))
         CALL taufromwfn(ddp,Grid,PAW%tphi(:,i),PAW%l(i))
         PAW%valetau=PAW%valetau+PAW%occ(i)*dp
         PAW%tvaletau=PAW%tvaletau+PAW%occ(i)*ddp
      endif   
    enddo
     dp=PAW%coretau+PAW%valetau
     ddp=PAW%tcoretau+PAW%tvaletau
     CALL exch(Grid,locald,vxc,sum,exc,&
&       tau=dp,vtau=PAW%vtau)

     CALL exch(Grid,localtd,tvxc,tsum,texc,&
&       tau=ddp,vtau=PAW%tvtau)

    write(std_out,*) 'tvtau exc texc ', exc, texc

! Kresse form    
    Kd=locald-PAW%core-localtd+PAW%tcore
    sum=integrator(Grid,Kd)

    write(std_out,*) 'compensation charge in Ktvtau ', sum
    Kd=localtd+sum*PAW%hatden

     CALL exch(Grid,Kd,Ktvxc,tsum,texc,&
&       tau=ddp,vtau=PAW%Ktvtau)



  open(1001,file='pseudovtau',form='formatted')
  write(1001,*)'#r  den    tildeden   Kd    tau    tildetau   rvxc    trvxc    Krtvxc    vtau     tvtau     Ktvtau'
     do i=1,Grid%n
        write(1001,'(1p, 50e17.7)') Grid%r(i),locald(i),localtd(i),Kd(i) &
&       ,dp(i),ddp(i),vxc(i),tvxc(i),Ktvxc(i),PAW%vtau(i),PAW%tvtau(i)  &
&       ,PAW%Ktvtau(i)
    enddo
 close(1001)

  open(1001,file='coredenandtau',form='formatted')
     do i=1,Grid%n
        write(1001,'(1p, 50e17.7)') Grid%r(i),PAW%core(i),PAW%tcore(i)&
&          ,PAW%coretau(i),PAW%tcoretau(i)
    enddo
 close(1001)


    deallocate(dp,ddp,vxc,tvxc,locald,localtd,Kd,Ktvxc)
  END SUBROUTINE calculate_tvtau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Smoothtvtau(Grid,PAW,range)
!      In the range  of radii  0 .. range,
!       consistently modify PAW%tvtau and PAW%Ktvtau to quadradic form
!         so that derivatives  PAW%dtvtaudr and PAW%Kdtvtaudr are well behaved
!         Call after calculate_tvtau
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Smoothtvtau(Grid,PAW,range)
    TYPE(GridInfo), INTENT(IN) :: Grid
    TYPE(PseudoInfo), INTENT(INOUT) :: PAW
    REAL(8) , INTENT(IN) :: range

    INTEGER :: i,j,k,l,nrin
    REAL(8), allocatable :: original(:),aux(:),Koriginal(:),Kaux(:)
    REAL(8), allocatable :: Ktvxc(:),Kd(:)
    REAL(8) :: fac,exc,texc,sum,tsum,rr,one,two,diff,Kdiff,ddiff,Kddiff,r1,r2,r3,r4,r5
    REAL(8), parameter :: small=1.d-5
    INTEGER :: times=-1
    CHARACTER(4) :: stuff

    allocate(original(Grid%n),aux(Grid%n),Koriginal(Grid%n),Kaux(Grid%n))

    times=times+1

    j=0
    do i=1,Grid%n
       if (Grid%r(i)>range) exit
       j=j+1
    enddo
    nrin=j
     write(std_out,*) 'in smoothetvtau -- times ', times
     write(std_out,*) 'in smoothtvtau --  range, Grid%r(nrin) ' , nrin,  range, Grid%r(nrin)
     r1=(Grid%r(nrin-2)/range)**2
     r2=(Grid%r(nrin-1)/range)**2
     r3=(Grid%r(nrin)/range)**2
     r4=(Grid%r(nrin+1)/range)**2
     r5=(Grid%r(nrin+2)/range)**2
     one=(PAW%tvtau(1)+PAW%tvtau(2)+PAW%tvtau(3)+PAW%tvtau(4)+PAW%tvtau(5))/5
     two=(PAW%tvtau(nrin-2)-one)/r1
     two=two+(PAW%tvtau(nrin-1)-one)/r2
     two=two+(PAW%tvtau(nrin)-one)/r3
     two=two+(PAW%tvtau(nrin+1)-one)/r4
     two=two+(PAW%tvtau(nrin+2)-one)/r5
     two=two/5
     write(std_out,*) 'in smoothtvtau --  tvtau -- one two  ', one,two
     original(:)=PAW%tvtau(:)
     call derivative(Grid,PAW%tvtau,aux)
     PAW%dtvtaudr(:)=aux(:)
     diff=0.d0; ddiff=0.d0
     do i=1,nrin
        PAW%tvtau(i)=one+two*(Grid%r(i)/range)**2
        diff=diff+(PAW%tvtau(i)-original(i))**2
        PAW%dtvtaudr(i)=2*two*(Grid%r(i)/range**2)
        ddiff=ddiff+(PAW%dtvtaudr(i)-aux(i))**2
     enddo
     write(std_out,*) 'in smoothtvtau --  tvtau -- diff  ddiff ', diff,ddiff
     one=(PAW%Ktvtau(1)+PAW%Ktvtau(2)+PAW%Ktvtau(3)+PAW%Ktvtau(4)+PAW%Ktvtau(5))/5
     two=(PAW%Ktvtau(nrin-2)-one)/r1
     two=two+(PAW%Ktvtau(nrin-1)-one)/r2
     two=two+(PAW%Ktvtau(nrin)-one)/r3
     two=two+(PAW%Ktvtau(nrin+1)-one)/r4
     two=two+(PAW%Ktvtau(nrin+2)-one)/r5
     two=two/5
     write(std_out,*) 'in smoothtvtau --  Ktvtau -- one two  ', one,two
     Koriginal(:)=PAW%Ktvtau(:)
     call derivative(Grid,PAW%Ktvtau,Kaux)
     PAW%Kdtvtaudr(:)=Kaux(:)
     Kdiff=0.d0; Kddiff=0.d0
     do i=1,nrin
        PAW%Ktvtau(i)=one+two*(Grid%r(i)/range)**2
        Kdiff=Kdiff+(PAW%Ktvtau(i)-Koriginal(i))**2
        PAW%Kdtvtaudr(i)=2*two*(Grid%r(i)/range**2)
        Kddiff=Kddiff+(PAW%Kdtvtaudr(i)-Kaux(i))**2
     enddo
     write(std_out,*) 'in smoothtvtau --  Ktvtau -- Kdiff  Kddiff ', Kdiff,Kddiff

     call mkname(times,stuff)
     open(1000,file='checktvtau'//TRIM(stuff),form='formatted')
     write(1000,*) 'r        tvtau      original       dtvtaudr       aux '
     do i=1,nrin+1
        write(1000,'(1p,9e15.7)') Grid%r(i),PAW%tvtau(i),original(i),PAW%dtvtaudr(i),aux(i) &
&             ,PAW%Ktvtau(i),Koriginal(i),PAW%Kdtvtaudr(i),Kaux(i)
     enddo
     close(1000)

    deallocate(original,aux,Koriginal,Kaux)

  END SUBROUTINE Smoothtvtau
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! makeprojectors_vtau
!     Designed for MGGA case
!     Needs PAW%rveff, PAW%tphi
!     Uses Vanderbilt scheme to develop projectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE makeprojectors_vtau(Grid)
    TYPE(GridInfo), INTENT(IN) :: Grid

    INTEGER :: i,j,k,l,io,jo,ok,lmax,nbase,n,irc,irc_vloc,nr,np,thisrc
    INTEGER :: icount,jcount,istart,ifinish,ibase,jbase,lprev
    INTEGER :: match=5
    REAL(8) :: choice,rc,xx,yy,gg,g,gp,gpp,al(2),ql(2),root1,root2,logderiv
    INTEGER, ALLOCATABLE :: omap(:)
    REAL(8), ALLOCATABLE :: VNC(:),aa(:,:),ai(:,:),aux(:)
    REAL(8), ALLOCATABLE :: f(:),d(:),dd(:)
    INTEGER :: LWORK
    REAL(8), POINTER  :: r(:)
    LOGICAL :: success

    n=Grid%n
    r=>Grid%r
    irc=PAW%irc
    irc_vloc=PAW%irc_vloc
    nbase=PAW%nbase
    lmax=PAW%lmax


  ! Set screened local pseudopotential
    allocate(VNC(n),aux(n),f(n),d(n),dd(n),stat=i)
    if (i/=0) stop 'allocation error in make_modrrkj'
    VNC(2:n)=PAW%rveff(2:n)/r(2:n)
    call extrapolate(Grid,VNC)
    !#call derivative(Grid,PAW%tvtau,aux)
 
    write(std_out,*) 'In program makeprojectors_vtau '; call flush_unit(std_out)
! Form unorthogonalized projectors
    do io=1,nbase
       f=PAW%tphi(:,io)
       l=PAW%l(io)
       xx=l*(l+1)
       call derivative(Grid,f,d)
       call derivative(Grid,d,dd)
       dd(2:n)=dd(2:n)-xx*f(2:n)/(Grid%r(2:n)**2)
       d(2:n)=d(2:n)-f(2:n)/Grid%r(2:n)
       call extrapolate(Grid,d)
       call extrapolate(Grid,dd)
       PAW%tp(:,io)=(PAW%eig(io)-VNC(:))*f(:)+&
&                 (1.d0+PAW%tvtau(:))*dd(:)+PAW%dtvtaudr(:)*d(:)
!       do i=1,n
!          write(100+io,'(1p,50e16.7)') Grid%r(i),f(i),d(i),dd(i),PAW%tp(i,io)
!       enddo   
    enddo   

!    stop
  !! Form orthogonalized projectors --   VANDERBILTORTHO
     write(std_out,*) ' Vanderbilt ortho'
     do l=0,lmax
       icount=0
       do io=1,nbase
        if (PAW%l(io)==l) icount=icount+1
       enddo
       if (icount==0) cycle
       write(std_out,*) 'For l = ', l,icount,' basis functions'
       allocate(aa(icount,icount),ai(icount,icount),omap(icount))
       aa=0;icount=0
       do io=1,nbase
        if (PAW%l(io)==l) then
          icount=icount+1
          omap(icount)=io
        endif
       enddo
       do i=1,icount
         io=omap(i)
         PAW%otphi(:,io)=PAW%tphi(:,io)
         PAW%ophi(:,io)=PAW%phi(:,io)
       enddo
       do i=1,icount
          io=omap(i)
         do j=1,icount
           jo=omap(j)
           aa(i,j)=overlap(Grid,PAW%otphi(:,io),PAW%tp(:,jo),1,irc)
         enddo
       enddo
       ai=aa;call minverse(ai,icount,icount,icount)

       do i=1,icount
         io=omap(i)
         PAW%ck(io)=ai(i,i)
         PAW%otp(:,io)=0
         do j=1,icount
           jo=omap(j)
           PAW%otp(:,io)=PAW%otp(:,io)+PAW%tp(:,jo)*ai(j,i)
         enddo
       enddo

       write(std_out,*) 'Check  otp for l = ', l
       do i = 1, icount
          io=omap(i)
          do j = 1, icount
             jo=omap(j)
             write(std_out,*) 'Overlap i j ', i,j, &
&                    overlap(Grid,PAW%otphi(:,io),PAW%otp(:,jo),1,irc)
          enddo
       enddo
       deallocate(aa,ai,omap)
    enddo

    open(1001,file='nonorthoproj',form='formatted')
      do i=1,irc+100
         write(1001,'(1p,50e16.7)') Grid%r(i),(PAW%tp(i,io),io=1,nbase)
      enddo   
    close(1001)  
    open(1001,file='orthoproj',form='formatted')
      do i=1,irc+100
         write(1001,'(1p,50e16.7)') Grid%r(i),(PAW%otp(i,io),io=1,nbase)
      enddo   
    close(1001)  
    do io=1,nbase
!wrong       write(std_out,*) 'For basis ', io, ' projector error ',SUM(abs(PAW%tp(irc_vloc:n,io)))
       write(std_out,*) 'For basis ', io, ' projector error ',SUM(abs(PAW%tp(irc:n,io)))
!wrong       PAW%tp(irc_vloc+1:n,io)=0.d0
       PAW%tp(irc+1:n,io)=0.d0
!wong       PAW%otp(irc_vloc+1:n,io)=0.d0
       PAW%otp(irc+1:n,io)=0.d0
    enddo
  deallocate(VNC,aux,f,d,dd)
  END SUBROUTINE makeprojectors_vtau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !modified version of SUBROUTINE makebasis_custom which was
  !           written by Marc Torrent from earlier version by NAWH
  !           for the case the vloc==0 (KS only; not HF)
  !           PAW%core and PAW%tcore already loaded
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE makebasis_V_setvloc(Grid,Pot,PAW)
      Type(GridInfo), INTENT(IN) :: Grid
      Type(PotentialInfo), INTENT(IN) :: Pot
      Type(PseudoInfo), INTENT(INOUT) :: PAW

      INTEGER ::  i,j,k,l,n,ib,jb,io,jo,nbase,irc,np,thisrc,lmax,ni
      INTEGER :: icount
      INTEGER, ALLOCATABLE :: omap(:)
      REAL(8) :: rc,gp,gpp,xx,occ,stuff
      REAL(8), POINTER :: r(:)
      REAL(8), ALLOCATABLE :: Ci(:),dum(:),tdum(:),hat(:),aa(:,:),ai(:,:)
      REAL(8), ALLOCATABLE :: dpdr(:),pdr(:)
      REAL(8), PARAMETER :: threshold=1.d-6

      if(needvtau) then
           write(std_out,*) 'In makebasis_V_setvloc -- needvtau not setup'
           stop
      endif   

      if (TRIM(PAW%exctype)=='HF'.or.TRIM(PAW%exctype)=='EXX') then
         write(std_out,*) 'makebasis_V_setvloc is not designed for ', PAW%exctype
         stop
      endif
      n=Grid%n
      r=>Grid%r
      irc=PAW%irc
      lmax=PAW%lmax
      nbase=PAW%nbase
      np=5
      allocate(Ci(np),dum(n),tdum(n),hat(n))
      PAW%rveff=PAW%vloc*Grid%r

  ! Loop on basis elements
      do ib=1,nbase
         l=PAW%l(ib)
     ! Read matching radius (from input dataset)
         thisrc=FindGridIndex(Grid,input_dataset%basis_func_rc(ib))
         thisrc=MIN(thisrc,irc)       ! make sure rc<total rc
         rc=r(thisrc);PAW%rcio(io)=rc
         WRITE(STD_OUT,'(a,3i5,1p,1e15.7)') ' For this wfn: ',ib,PAW%np(ib),PAW%l(ib),PAW%eig(ib)
         WRITE(STD_OUT,'(a,f10.7)') '  >>> rc =', rc
         if (thisrc<3.or.thisrc>irc.or.thisrc>n-3) then
            write(std_out,*) 'rc out of range', thisrc,n,irc
            stop
         endif
         call pspolyn(PAW%phi(:,ib),Ci,r,l,np,thisrc,n)
     ! Compute pseudized partial wave and part of  projector
         PAW%tphi(:,ib)=PAW%phi(:,ib)
         PAW%tp(:,ib)=0.d0 ; PAW%Kop(:,ib)=0
         do i=1,thisrc-1
            xx=r(i)*r(i)
            PAW%tphi(i,ib)=Ci(1)+Ci(2)*xx
            gp=2.d0*Ci(2)
            gpp=2.d0*Ci(2)
            do j=3,np
               PAW%tphi(i,ib)=PAW%tphi(i,ib)+Ci(j)*xx**(j-1)
               gp=gp+dble(2*j-2)*Ci( j)*xx**(j-2)
               gpp=gpp+dble((2*j-2)*(2*j-3))*Ci(j)*xx**(j-2)
            enddo
            PAW%tphi(i,ib)=PAW%tphi(i,ib)*r(i)**(l+1)
            PAW%tp(i,ib)=-(dble(2*(l+1))*gp+gpp)*r(i)**(l+1)  !kinetic pt.
            PAW%Kop(i,ib)=PAW%tp(i,ib)
          enddo
          do i=thisrc,n
             gpp=Gsecondderiv(Grid,i,PAW%tphi(:,ib))
             PAW%Kop(i,ib)=-gpp+dble(l*(l+1))/(r(i)**2)*PAW%tphi(i,ib)
             if (i<irc) PAW%tp(i,ib)=PAW%Kop(i,ib)
          enddo
       enddo     !nbase

    do io=1,PAW%OCCWFN%norbit
       if(PAW%valencemap(io)>0) then
          ib=PAW%valencemap(io)
          PAW%TOCCWFN%wfn(:,io)=PAW%tphi(:,ib)
          !PAW%OCCWFN%wfn(:,io)=PAW%phi(:,ib)    !presumably also true
       endif
    enddo

  ! Find rveff
   PAW%den=0.d0; PAW%tden=0.d0
   allocate(dpdr(n),pdr(n))
   do io=1,PAW%OCCWFN%norbit
      if (.not.PAW%OCCwfn%iscore(io)) then
         occ=PAW%OCCwfn%occ(io)
          PAW%den=PAW%den+occ*PAW%OCCwfn%wfn(:,io)**2
          PAW%tden=PAW%tden+occ*PAW%TOCCwfn%wfn(:,io)**2
      endif
   enddo
   dum=PAW%core+PAW%den-PAW%tcore-PAW%tden
   xx=-Pot%nz+integrator(Grid,dum)
   write(std_out,*) 'Checking charge for pseudo scheme ', xx,integrator(Grid,tdum)
   PAW%rveff=PAW%rveff+xx*PAW%hatpot
   tdum=PAW%tcore+PAW%tden
   call poisson(Grid,gp,tdum,dum,stuff)    ! Coulomb from tcore and tden
   write(std_out,*) 'After Poisson ', gp,stuff
   PAW%rveff=PAW%rveff+dum
      CALL exch(Grid,tdum,hat,gp,stuff)
   PAW%rveff=PAW%rveff+hat
   do i=1,n
      write(901,'(1p,20e15.7)') Grid%r(i),PAW%rveff(i),PAW%AErefrv(i),dum(i),&
&              xx*PAW%hatpot(i),hat(i)
   enddo
   tdum=0; tdum(2:n)=PAW%rveff(2:n)/r(2:n)
  ! complete projector functions  ; at this point PAW%tp stores -K*tilde{wfn}
   do ib=1,nbase
      PAW%tp(1:irc-1,ib)=PAW%tp(1:irc-1,ib)+&
&              (tdum(1:irc-1)-PAW%eig(ib))*PAW%tphi(1:irc-1,ib)
   enddo

     ! Form orthogonalized projector functions
    do ib=1,nbase
       PAW%ophi(:,ib)=PAW%phi(:,ib)
       PAW%otphi(:,ib)=PAW%tphi(:,ib)
    enddo

    do l=0,lmax
       icount=0
       do ib=1,nbase
          if (PAW%l(ib)==l) icount=icount+1
       enddo
       if (icount==0) cycle
       write(std_out,*) 'For l = ', l,icount,' basis functions'
       allocate(aa(icount,icount),ai(icount,icount),omap(icount))
       aa=0;icount=0
       do ib=1,nbase
          if (PAW%l(ib)==l) then
             icount=icount+1
             omap(icount)=ib
          endif
       enddo
       do i=1,icount
          ib=omap(i)
          do j=1,icount
             jb=omap(j)
             aa(i,j)=overlap(Grid,PAW%otphi(:,ib),PAW%tp(:,jb),1,irc)
          enddo
       enddo
       ai=aa;call minverse(ai,icount,icount,icount)

       do i=1,icount
          ib=omap(i)
          PAW%ck(ib)=ai(i,i)
          PAW%otp(:,ib)=0
          do j=1,icount
             jb=omap(j)
             PAW%otp(:,ib)=PAW%otp(:,ib)+PAW%tp(:,jb)*ai(j,i)
          enddo
       enddo
       deallocate(aa,ai,omap)
     enddo

    do i=1,n
       write(801,'(1p,50e15.7)') r(i),(PAW%otp(i,ib),PAW%tp(i,ib),ib=1,nbase)
    enddo

   deallocate(Ci,dum,tdum,hat)
 END SUBROUTINE makebasis_V_setvloc

    !*************************************************************************
    !  program to generate projector functions for Blochl's paw formalism
    !    starting with smooth functions
    !    for every basis function phi, choose smooth function with
    !     form for r<rc: (r**(l+1))*sum(n)*(cn*(r**(2*n)))
    !       where, rc==r(iiirc)
    !  phi(i,io), eval(io), tocc(io)  original basis function, energy, occupancy
    !  tphi(i,io) smooth basis function corresponding to phi(i,io)
    !  tp(i,io)=constant*normalized Harmonic Oscillator function
    !  tp == 0 for r>r(iiirc)
    !  functions defined to be identically zero for r>r(iiirc))
    !*************************************************************************
    SUBROUTINE formprojectors(Grid,Pot,PS,option)
      TYPE(GridInfo),  INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot,PS
      INTEGER,INTENT(IN) :: option

      INTEGER :: nbase,lmax,l,io,irc,wantednodes,nb,n,icount
      INTEGER ::  istart,ifinish,ibase,jbase
      REAL(8) :: h,xx,yy,choice,rc
      INTEGER,allocatable :: irc_io(:)

      n=Grid%n
      h=Grid%h
      lmax=PAW%lmax
      nbase=PAW%nbase
      irc=PAW%irc
      !
      !   form projector functions for each l
      !

      !Read rc from input dataset
      if (option==1) then
       allocate(irc_io(nbase))
       do io=1,nbase
        irc_io(io)=FindGridIndex(Grid,input_dataset%basis_func_rc(io))
        rc=Grid%r(irc_io(io))
        if(irc_io(io)>PAW%irc) then
         write(std_out,*) 'rc out of range', irc_io(io),n,PAW%irc
         stop
        endif
       enddo
      endif

      DO l=0,lmax
         icount=0
         DO io=1,nbase
            IF (PAW%l(io)==l) THEN
               IF  (icount==0) istart=io
               IF  (icount >= 0) ifinish=io
               icount=icount+1
               wantednodes=icount-1
               ! form unorthonormalized projector functions tp
               WRITE(STD_OUT,*) '******* projector for l = ',l
               if (option==1) then
                CALL  bsolv(Grid,PS,io,wantednodes,irc_io(io))
               else
                CALL  bsolv(Grid,PS,io,wantednodes)
               endif
               PAW%ophi(:,io)=PAW%phi(:,io)
               PAW%otphi(:,io)=PAW%tphi(:,io)
               PAW%otp(:,io)=PAW%tp(:,io)
               PAW%Kop(1,io)=0
               PAW%Kop(2:n,io)=(PAW%eig(io)-Pot%rv(2:n)/Grid%r(2:n))&
&                       *PAW%phi(2:n,io)
            ENDIF
         ENDDO
         !write(std_out,*) 'orthnormalization'
         !write(std_out,*) 'start orthogonalization',istart,ifinish
         DO ibase=istart,ifinish
            DO jbase=istart,ibase
               IF (jbase.LT.ibase) THEN
                  xx=overlap(Grid,PAW%otp(:,jbase),PAW%otphi(:,ibase),1,irc)
                  yy=overlap(Grid,PAW%otphi(:,jbase),PAW%otp(:,ibase),1,irc)
                  !write(std_out,*) 'before',jbase,ibase,xx,yy
                  PAW%ophi(1:n,ibase)=PAW%ophi(1:n,ibase)-PAW%ophi(1:n,jbase)*xx
                  PAW%Kop(1:n,ibase)=PAW%Kop(1:n,ibase)-PAW%Kop(1:n,jbase)*xx
                  PAW%otphi(1:n,ibase)=PAW%otphi(1:n,ibase)-PAW%otphi(1:n,jbase)*xx
                  PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)-PAW%otp(1:n,jbase)*yy
               ELSE IF (jbase.EQ.ibase) THEN
                  xx=overlap(Grid,PAW%otp(:,ibase),PAW%otphi(:,ibase),1,irc)
                !write(std_out,*) 'before',jbase,ibase,xx
                  choice=1.d0/SQRT(ABS(xx))
                  PAW%otp(1:n,ibase)=PAW%otp(1:n,ibase)*DSIGN(choice,xx)
                  PAW%otphi(1:n,ibase)=PAW%otphi(1:n,ibase)*choice
                  PAW%ophi(1:n,ibase)=PAW%ophi(1:n,ibase)*choice
                  PAW%Kop(1:n,ibase)=PAW%Kop(1:n,ibase)*choice
               ENDIF
            ENDDO
         ENDDO

      ENDDO

      if (option==1) deallocate(irc_io)

    END SUBROUTINE formprojectors
    !*************************************************************************
    !  on input tphi=phi
    !  on output tphi recalculated for r<nrc*h
    !  l = angular momentum quantum number
    !  nodes = the number of desired nodes in tp and tphi
    !  shape function
    !*************************************************************************
    SUBROUTINE bsolv(Grid,Pot,io,wantednodes,irc_io)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      INTEGER, INTENT(IN) :: io,wantednodes
      INTEGER, OPTIONAL :: irc_io

      REAL(8), ALLOCATABLE :: f(:),chi(:),fakerv(:)
      REAL(8), ALLOCATABLE :: tphi(:),tp(:)
      REAL(8), POINTER :: r(:),rv(:)
      REAL(8) :: h,en,v0,v0p,ol,h2,eh2,veh2,xnorm,rc_io
      INTEGER ::  l,n,num,nrc,i,j,k,ii,jj,iter,irc,nodes,ok
      REAL(8) :: cpmin,cpmax,zeroval
      REAL(8) :: cp,cph2,del,deriv,tderiv
      REAL(8), PARAMETER :: small=1.d-10,step=1.d0
      INTEGER, PARAMETER :: mxiter=1000

      n=Grid%n
      h=Grid%h
      r=>Grid%r
      rv=>Pot%rv
      if (present(irc_io)) then
       irc=irc_io
      else
       irc=PAW%irc
      endif
      ALLOCATE(fakerv(n),f(n),chi(n),tphi(n),tp(n),stat=ok)
      IF (ok /= 0) THEN
         WRITE(STD_OUT,*) 'Error in bsolv allocation ', n,ok
         STOP
      ENDIF
      tphi=PAW%phi(:,io)
      en=PAW%eig(io)
      l=PAW%l(io)
      deriv=(Gfirstderiv(Grid,irc,tphi))/tphi(irc) ! log derivative
      tderiv=0    ! initial log derive
      cp=0         ! initial projector constant
      del=1000     ! initial error
      cpmin=-100
      cpmax=100
      if (present(irc_io)) then
       rc_io=r(irc_io)
       f(1)=1.d0;f(irc_io:n)=0.d0
       DO i=2,irc_io-1
        f(i)=(SIN(pi*r(i)/rc_io)/(pi*r(i)/rc_io))**2
       ENDDO
      else
       DO i=1,n
          f(i)=PAW%projshape(i)
       ENDDO
      endif

      WRITE(STD_OUT,*) 'in bsolv -- l, en, n',l,en,wantednodes

      iter=0

      DO

         iter=iter+1

         WRITE(STD_OUT,*) 'bsolv iter cp',iter,cp

         fakerv=rv-cp*f*Grid%r
         call zeropot(Grid,fakerv,v0,v0p)
         ! initialize chi
         chi=0
         chi(2)=wfninit(0.d0,l,v0,v0p,en,Grid%r(2))
         zeroval=0
         if (l==1) zeroval=2

         call forward_numerov(Grid,l,irc+5,en,fakerv,zeroval,chi,nodes)

         WRITE(STD_OUT,'("iter nodes cp cpmin cpmax",2i5,1p,3e15.7)')&
&             iter,nodes,cp,cpmin,cpmax
         IF (nodes.EQ.wantednodes) THEN
            tderiv=(Gfirstderiv(Grid,irc,chi))/chi(irc) ! log derivative
            tp=0
            tp(1:irc)=chi(1:irc)/chi(irc)
            chi(1:irc)=f(1:irc)*(tp(1:irc))**2
            xnorm=integrator(Grid,chi(1:irc),1,irc)
            del=(tderiv-deriv)/xnorm
            WRITE(STD_OUT,*) 'iter nodes del', iter,nodes,del

            IF (ABS(del).LT.small) EXIT

            IF (iter.GE.mxiter) THEN
               WRITE(STD_OUT,*)' terminating projector',iter
               STOP
            ENDIF

            IF (ABS(del).GT.step) del=DSIGN(step,del)
            cp=cp+del
            IF (cp.GT.cpmax) cp=cpmax-ranx()*step
            IF (cp.LT.cpmin) cp=cpmin+ranx()*step

         ELSE IF(nodes.GT.wantednodes) THEN
            cpmax=cp
            cp=cpmax-ranx()*step
         ELSE IF(nodes.LT.wantednodes) THEN
            cpmin=cp
            cp=cpmin+ranx()*step
         ENDIF

      ENDDO

      tphi(1:irc-1)=tphi(irc)*tp(1:irc-1)
      tphi(irc:n)=PAW%phi(irc:n,io)
      tp(1:irc)=f(1:irc)*tphi(1:irc)
      chi(1:irc)=tphi(1:irc)*tp(1:irc)
      xnorm=integrator(Grid,chi(1:irc),1,irc)
      WRITE(STD_OUT,*) 'normalization for projector l,n=',l,nodes,xnorm
      tp(1:irc)=tp(1:irc)/xnorm
      tp(irc+1:n)=0

      PAW%tphi(:,io)=tphi
      PAW%tp(:,io)=tp
      PAW%ck(io)=cp

      WRITE(STD_OUT,*) 'completed bsolv',io,cp

      DEALLOCATE(fakerv,f,chi,tphi,tp)
    END SUBROUTINE bsolv

    !***********************************************************************
    !  program to calculate Fourier transform of paw product tp*tphi
    !   in order to get an idea of the convergence
    !   and output them to files tprod.l
    !***********************************************************************
    SUBROUTINE ftprod(Grid)
      TYPE(GridInfo), INTENT(IN) :: Grid

      INTEGER, PARAMETER :: nq=200
      REAL(8), PARAMETER :: qmax=15.d0

      REAL(8), ALLOCATABLE :: q(:),dum(:),dum1(:)
      REAL(8), POINTER :: r(:)
      REAL(8) :: h,dq,tphij
      INTEGER :: i,ib,l,iq,n,irc
      CHARACTER(4) flnm
      !
      WRITE(STD_OUT,*) 'calculating Fourier transforms of tp*tphi products  ',&
&          'For bound states only '

      n=Grid%n
      h=Grid%h
      r=>Grid%r
      irc=PAW%irc
      ALLOCATE(q(nq),dum(n),dum1(n),stat=i)
      IF (i/=0) THEN
         WRITE(STD_OUT,*) 'Error in allocating space in ftprod',n,nq,i
         STOP
      ENDIF

      dq=qmax/nq
      DO ib=1,PAW%nbase
         IF (PAW%eig(ib)< 0.d0) THEN
            l=PAW%l(ib)
            CALL mkname(ib,flnm)
            OPEN(55,file='tprod.'//TRIM(flnm),form='formatted')
            DO iq=1,nq
               q(iq)=iq*dq
               DO i=1,n
                  dum(i)=q(iq)*r(i)
               ENDDO
               CALL sphbes(l,n,dum)
               DO i=1,n
                  dum1(i)=dum(i)*r(i)*PAW%otphi(i,ib)
                  IF (i.LE.irc) dum(i)=dum(i)*r(i)*PAW%otp(i,ib)
               ENDDO
               tphij=integrator(Grid,dum1(1:n))*&
&                    integrator(Grid,dum(1:irc),1,irc)

               WRITE(55,'(1p,2e16.7)') q(iq),tphij
            ENDDO
            CLOSE(55)
         ENDIF  !(e<=0)
      ENDDO  ! ib
      DEALLOCATE(q,dum,dum1)
    END SUBROUTINE ftprod
    !***********************************************************************
    !  program to calculate Fourier transform of hatpot functions
    !   in order to get an idea of the convergence
    !   and output them to files hatpot.l
    !***********************************************************************
    SUBROUTINE fthatpot(Grid)
      TYPE(GridInfo), INTENT(IN) :: Grid

      INTEGER, PARAMETER :: nq=200
      REAL(8), PARAMETER :: qmax=15.d0

      REAL(8), ALLOCATABLE :: q(:),dum(:),dum1(:),dum2(:),dum3(:),arg(:)
      REAL(8), POINTER :: r(:)
      REAL(8) :: h,dq,fthatp,fthatd,fthattd
      INTEGER :: i,ib,l,iq,n,irc,ll
      CHARACTER(4) flnm
      !
      WRITE(STD_OUT,*) 'calculating Fourier transforms of hatpot for  each l '

      n=Grid%n
      h=Grid%h
      r=>Grid%r
      irc=PAW%irc
      ALLOCATE(q(nq),dum(n),dum1(n),dum2(n),dum3(n),arg(n),stat=i)
      IF (i/=0) THEN
         WRITE(STD_OUT,*) 'Error in allocating space in fthatpot',n,nq,i
         STOP
      ENDIF

      dq=qmax/nq
      ll=2*MAXVAL(PAW%l(:))
      DO l=0,ll
         CALL mkname(l,flnm)
         OPEN(55,file='fthatpot.'//TRIM(flnm),form='formatted')
         CALL hatL(Grid,PAW,l,dum2)
         DO iq=1,nq
            q(iq)=iq*dq
            DO i=1,n
               dum(i)=q(iq)*r(i)
            ENDDO
            CALL sphbes(l,n,dum)
            DO i=1,n
               dum3(i)=dum(i)*dum2(i)
            ENDDO
            fthatd=integrator(Grid,dum3(1:n))
            fthatp=8*PI*fthatd/q(iq)**2

            fthattd=0
            IF(l==0) THEN
               arg(1:n)=dum(1:n)*PAW%tden(1:n)
               fthattd=integrator(Grid,arg)
            ENDIF

            WRITE(55,'(1p,5e16.7)') q(iq),fthatp,fthatd,fthatp*fthatd,fthatp*fthattd
         ENDDO
         CLOSE(55)
      ENDDO  ! l
      DEALLOCATE(q,dum,dum1,dum2,dum3,arg)
    END SUBROUTINE fthatpot

    !***********************************************************************
    !  program to calculate Fourier transform of tphi*q^2
    !   in order to get an idea of the convergence of kinetic energy
    !   and output them to files ftkin.l
    !***********************************************************************
    SUBROUTINE ftkin(Grid)
      TYPE(GridInfo), INTENT(IN) :: Grid

      INTEGER, PARAMETER :: nq=200
      REAL(8), PARAMETER :: qmax=15.d0

      REAL(8), ALLOCATABLE :: q(:),dum(:),dum1(:)
      REAL(8), POINTER :: r(:)
      REAL(8) :: h,dq,kin
      INTEGER :: i,ib,l,iq,n,irc
      CHARACTER(4) flnm
      !
      WRITE(STD_OUT,*) 'calculating Fourier transforms of tphi  ',&
&          'For bound states only '

      n=Grid%n
      h=Grid%h
      r=>Grid%r
      ALLOCATE(q(nq),dum(n),dum1(n),stat=i)
      IF (i/=0) THEN
         WRITE(STD_OUT,*) 'Error in allocating space in ftkin',n,nq,i
         STOP
      ENDIF

      dq=qmax/nq
      DO ib=1,PAW%nbase
         IF (PAW%eig(ib)<=0.d0) THEN
            l=PAW%l(ib)
            CALL mkname(ib,flnm)
            OPEN(55,file='ftkin.'//TRIM(flnm),form='formatted')
            DO iq=1,nq
               q(iq)=iq*dq
               DO i=1,n
                  dum(i)=q(iq)*r(i)
               ENDDO
               CALL sphbes(l,n,dum)
               DO i=1,n
                  dum1(i)=dum(i)*r(i)*PAW%tphi(i,ib)
               ENDDO
               kin=integrator(Grid,dum1(1:n))*(q(iq))
               kin=kin*kin

               WRITE(55,'(1p,2e16.7)') q(iq),kin
            ENDDO
            CLOSE(55)
         ENDIF  !(e<=0)
      ENDDO  ! ib
      DEALLOCATE(q,dum,dum1)
    END SUBROUTINE ftkin

    !***********************************************************************
    !  program to calculate Fourier transform of vloc and  tden
    !   in order to get an idea of their convergence
    !   and output them to files ftvloc
    !***********************************************************************
    SUBROUTINE ftvloc(Grid)
      TYPE(GridInfo), INTENT(IN) :: Grid

      INTEGER, PARAMETER :: nq=200
      REAL(8), PARAMETER :: qmax=15.d0

      REAL(8), ALLOCATABLE :: q(:),dum(:),dum1(:)
      REAL(8), POINTER :: r(:)
      REAL(8) :: h,dq,vloc,tden
      INTEGER :: i,ib,l,iq,n,irc
      CHARACTER(4) flnm
      !
      WRITE(STD_OUT,*) 'calculating Fourier transforms of vloc and tden'

      n=Grid%n
      h=Grid%h
      irc=PAW%irc
      r=>Grid%r
      ALLOCATE(q(nq),dum(n),dum1(n),stat=i)
      IF (i/=0) THEN
         WRITE(STD_OUT,*) 'Error in allocating space in ftvloc',n,nq,i
         STOP
      ENDIF

      OPEN(55,file='ftvloc',form='formatted')
      dq=qmax/nq
      l=0
      DO iq=1,nq
         q(iq)=iq*dq
         DO i=1,n
            dum(i)=q(iq)*r(i)
         ENDDO
         CALL sphbes(l,n,dum)
         DO i=1,n
            dum1(i)=dum(i)*PAW%vloc(i)*(r(i)**2)
            dum(i)=dum(i)*PAW%tden(i)
         ENDDO
         vloc=integrator(Grid,dum1(1:n),1,irc)
         tden=integrator(Grid,dum(1:n))

         WRITE(55,'(1p,4e16.7)') q(iq),vloc,tden,vloc*tden
      ENDDO
      CLOSE(55)

      DEALLOCATE(q,dum,dum1)
    END SUBROUTINE ftvloc


    !***************************************************************************
    !  pgm to solve separable radial schroedinger equation
    !    at energy 'energy' and at angular momentum l
    !
    !    with smooth potential rveff/r, given in uniform mesh of n points
    !   r=i*h, i=1,...n-1 ;assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
    !                               p((n+1)*h)=0
    !
    !  uses Noumerov algorithm
    !
    !  For l=0,1 corrections are needed to approximate wfn(r=0)
    !     These depend upon:
    !         e0 (current guess of energy eigenvalue)
    !         l,nz==0
    !         v(0) == v0 electronic potential at r=0
    !         v'(0) == v0p derivative of electronic potential at r=0
    !
    ! also returns node == number of nodes for calculated state
    !
    !  proj == projector functions
    !  hij and qij == hamiltonianian and overlap matrix elements
    !***************************************************************************
    SUBROUTINE unboundsep(Grid,Pot,PAW,nr,l,energy,wfn,node)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: l,nr
      REAL(8), INTENT(IN) :: energy
      REAL(8), INTENT(INOUT) :: wfn(:)
      INTEGER, INTENT(INOUT) :: node

      INTEGER :: n,ia,ib,ic,nbase,icount,jcount,lcount,irc
      REAL(8) :: summ,h,scale,zeroval
      REAL(8), ALLOCATABLE :: y(:,:),b(:),a(:,:)

      n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc
      !write(std_out,*) 'Entering unboundsep with l energy = ', l, energy

      IF (nr<irc) THEN
         WRITE(STD_OUT,*) 'error in unboundsep -- nr < irc', nr,irc
         STOP
      ENDIF
      !
      ! initialize wfn
      wfn=0
      wfn(2)=wfninit(0.d0,l,Pot%v0,Pot%v0p,energy,Grid%r(2))
      zeroval=0
      if (l==1) zeroval=2

      call forward_numerov(Grid,l,nr,energy,Pot%rv,zeroval,wfn,node)

      ALLOCATE(y(nr,nbase),b(nbase),stat=ib)
      IF (ib/=0) THEN
         WRITE(STD_OUT,*) 'Error in unboundsep allocation',  nr,nbase,ib
         STOP
      ENDIF

      lcount=0
      DO ib=1,nbase
         IF (l==PAW%l(ib)) THEN
            lcount=lcount+1
            CALL inhomogeneous_numerov(Grid,l,nr,energy,&
&              PAW%rveff,PAW%otp(:,ib),y(:,lcount))
         ENDIF
      ENDDO

      IF(lcount>0) THEN
         ALLOCATE(a(lcount,lcount),stat=ib)
         IF (ib/=0) THEN
            WRITE(STD_OUT,*) 'Error in unboundsep allocation',  lcount,ib
            STOP
         ENDIF

         icount=0
         DO ib=1,nbase
            summ=0.d0
            IF (l==PAW%l(ib)) THEN
               icount=icount+1
               DO ic=1,nbase
                  IF (l==PAW%l(ic)) THEN
                     summ=summ+(PAW%dij(ib,ic)-energy*PAW%oij(ib,ic))*&
&                         overlap(Grid,PAW%otp(:,ic),wfn,1,irc)
                  ENDIF
               ENDDO
               b(icount)=-summ
            ENDIF
         ENDDO
         !
         icount=0
         DO ia=1,nbase
            IF (l==PAW%l(ia)) THEN
               icount=icount+1
               jcount=0
               DO ib=1,nbase
                  IF (l==PAW%l(ib)) THEN
                     jcount=jcount+1
                     summ=0.d0
                     IF (ia.EQ.ib) summ=1.d0
                     DO ic=1,nbase
                        IF (l==PAW%l(ic)) THEN
                           summ=summ+(PAW%dij(ia,ic)-energy*PAW%oij(ia,ic))*&
&                          overlap(Grid,PAW%otp(:,ic),y(:,jcount),1,irc)
                        ENDIF
                     ENDDO
                     a(icount,jcount)=summ
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         !
         CALL linsol(a,b,lcount,lcount,lcount,nbase)

         icount=0
         DO ib=1,nbase
            IF(l==PAw%l(ib)) THEN
               icount=icount+1
               wfn(1:nr)=wfn(1:nr)+b(icount)*y(1:nr,icount)
            ENDIF
         ENDDO
         DEALLOCATE(a)
      ENDIF
      !
      ! normalize to unity within integration range
      !
      scale=1.d0/sepnorm(Grid,PAW,nr,l,wfn)
      IF (scale.LE.0.d0) THEN
         WRITE(STD_OUT,*) 'warning -- negative norm for l=',l
         scale=-scale
         IF (scale.EQ.0.d0) scale=1.d0
      ENDIF
      scale=DSIGN(SQRT(scale),wfn(nr-2))
      wfn(1:nr)=wfn(1:nr)*scale
      DEALLOCATE(b,y)
    END SUBROUTINE unboundsep

    !***************************************************************************
    !  Modified for meta gga  and extra terms
    !
    !  pgm to solve separable radial schroedinger equation
    !    at energy 'energy' and at angular momentum l
    !
    !    with smooth potential rveff/r, given in uniform mesh of n points
    !   r=i*h, i=1,...n-1 ;assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
    !                               p((n+1)*h)=0
    !
    !  uses inhomocfdsol
    !
    !  For l=0,1 corrections are needed to approximate wfn(r=0)
    !     These depend upon:
    !         e0 (current guess of energy eigenvalue)
    !         l,nz==0
    !         v(0) == v0 electronic potential at r=0
    !         v'(0) == v0p derivative of electronic potential at r=0
    !
    ! also returns node == number of nodes for calculated state
    !
    !  proj == projector functions
    !  hij and qij == hamiltonianian and overlap matrix elements
    !***************************************************************************
    SUBROUTINE unboundsepked(Grid,Pot,PAW,nr,l,energy,wfn,node)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: l,nr
      REAL(8), INTENT(IN) :: energy
      REAL(8), INTENT(INOUT) :: wfn(:)
      INTEGER, INTENT(INOUT) :: node

      INTEGER :: n,ia,ib,ic,nbase,icount,jcount,lcount,irc
      REAL(8) :: summ,h,scale,zeroval
      REAL(8), ALLOCATABLE :: y(:,:),b(:),a(:,:)

      n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc

      IF (nr<irc) THEN
         WRITE(STD_OUT,*) 'error in unboundsepked -- nr < irc', nr,irc
         STOP
      ENDIF
      !

      call unboundked(Grid,Pot,nr,l,energy,wfn,node)

      ALLOCATE(y(nr,nbase),b(nbase),stat=ib)
      IF (ib/=0) THEN
         WRITE(STD_OUT,*) 'Error in unboundsepked allocation',  nr,nbase,ib
         STOP
      ENDIF

      lcount=0
      DO ib=1,nbase
         IF (l==PAW%l(ib)) THEN
            lcount=lcount+1
            CALL unboundkedinhomo(Grid,Pot,PAW%otp(1:nr,ib),nr,l,energy,&
&              y(1:nr,lcount))
         ENDIF
      ENDDO

      IF(lcount>0) THEN
         ALLOCATE(a(lcount,lcount),stat=ib)
         IF (ib/=0) THEN
            WRITE(STD_OUT,*) 'Error in unboundsepked allocation',  lcount,ib
            STOP
         ENDIF

         icount=0
         DO ib=1,nbase
            summ=0.d0
            IF (l==PAW%l(ib)) THEN
               icount=icount+1
               DO ic=1,nbase
                  IF (l==PAW%l(ic)) THEN
                     summ=summ+(PAW%dij(ib,ic)-energy*PAW%oij(ib,ic))*&
&                         overlap(Grid,PAW%otp(:,ic),wfn,1,irc)
                  ENDIF
               ENDDO
               b(icount)=-summ
            ENDIF
         ENDDO
         !
         icount=0
         DO ia=1,nbase
            IF (l==PAW%l(ia)) THEN
               icount=icount+1
               jcount=0
               DO ib=1,nbase
                  IF (l==PAW%l(ib)) THEN
                     jcount=jcount+1
                     summ=0.d0
                     IF (ia.EQ.ib) summ=1.d0
                     DO ic=1,nbase
                        IF (l==PAW%l(ic)) THEN
                           summ=summ+(PAW%dij(ia,ic)-energy*PAW%oij(ia,ic))*&
&                          overlap(Grid,PAW%otp(:,ic),y(:,jcount),1,irc)
                        ENDIF
                     ENDDO
                     a(icount,jcount)=summ
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         !
         CALL linsol(a,b,lcount,lcount,lcount,nbase)

         icount=0
         DO ib=1,nbase
            IF(l==PAw%l(ib)) THEN
               icount=icount+1
               wfn(1:nr)=wfn(1:nr)+b(icount)*y(1:nr,icount)
            ENDIF
         ENDDO
         DEALLOCATE(a)
      ENDIF
      !
      ! normalize to unity within integration range
      !
      scale=1.d0/sepnorm(Grid,PAW,nr,l,wfn)
      IF (scale.LE.0.d0) THEN
         WRITE(STD_OUT,*) 'warning -- negative norm for l=',l
         scale=-scale
         IF (scale.EQ.0.d0) scale=1.d0
      ENDIF
      scale=DSIGN(SQRT(scale),wfn(nr-2))
      wfn(1:nr)=wfn(1:nr)*scale
      DEALLOCATE(b,y)
    END SUBROUTINE unboundsepked

    !***************************************************************************
    !  pgm to solve separable radial schroedinger equation
    !    for bound state near energy 'energy' and at angular momentum l
    !
    !    with smooth potential rveff/r, given in uniform mesh of n points
    !   r=i*h, i=1,...n-1 ;assuming p(r)=C*r**(l+1)*polynomial(r) for r==0;
    !                               p((n+1)*h)=0
    !
    !  uses Noumerov algorithm
    !
    !  For l=0,1 corrections are needed to approximate wfn(r=0)
    !     These depend upon:
    !         e0 (current guess of energy eigenvalue)
    !         l,nz==0
    !         v(0) == v0 electronic potential at r=0
    !         v'(0) == v0p derivative of electronic potential at r=0
    !
    ! also returns node == number of nodes for calculated state
    !
    !  proj == projector functions
    !  hij and qij == hamiltonianian and overlap matrix elements
    !***************************************************************************
    SUBROUTINE boundsep(Grid,Pot,PAW,l,node,energy,emin,emax,wfn)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: l,node
      REAL(8), INTENT(INOUT) :: energy,emin,emax
      REAL(8), INTENT(INOUT) :: wfn(:)

      INTEGER, PARAMETER :: niter=100
      INTEGER :: n,i,j,ia,ib,ic,nbase,icount,jcount,lcount
      INTEGER :: match,mmatch,irc,node1,iter
      REAL(8), PARAMETER :: convre=1.d-10,err=1.d-9
      REAL(8) :: summ,h,scale,best,energy0,dele,x,rin,rout,zeroval,ppp
      REAL(8), ALLOCATABLE :: p1(:),u(:),y(:,:),b(:),a(:,:)

      n=Grid%n; h=Grid%h; nbase=PAW%nbase;  irc=PAW%irc
      !
      ALLOCATE(p1(n),u(n),y(n,nbase),b(nbase),stat=ib)
      IF (ib/=0) THEN
         WRITE(STD_OUT,*) 'Error in boundsep allocation',  n,nbase,ib
         STOP
      ENDIF
      !WRITE(STD_OUT,*) 'in boundsep with ', l,node,energy,emin,emax
      lcount=0
      DO ib=1,nbase
         IF (l==PAW%l(ib)) THEN
            lcount=lcount+1
         ENDIF
      ENDDO
      !
      IF(lcount>0) THEN
         ALLOCATE(a(lcount,lcount),stat=ib)
         IF (ib/=0) THEN
            WRITE(STD_OUT,*) 'Error in boundsep allocation',  lcount,ib
            STOP
         ENDIF
      ENDIF

      IF (emax.GT.0.d0) emax=0.d0
      best=1.d10
      IF (energy.LT.emin) energy=emin+err
      IF (energy.GT.emax) energy=emax
      energy0=energy

      iter=0


      DO
         match=MIN(irc+10,n-10)
         x=0.5d0*(Pot%rv(n)/Grid%r(n)+Pot%rv(n-1)/Grid%r(n-1))+&
&               l*(l+1)/(Grid%r(n)**2)
         ppp=SQRT(ABS(x-energy))
         p1=0
         p1(n)=1
         p1(n-1)=exp(-ppp*(Grid%r(n-1)-Grid%r(n)))

         !write(std_out,*) 'before backward', n,p1(n-1),p1(n)
         !write(std_out,*) 'before backward', x,ppp,exp(-ppp*(Grid%r(n-1)-Grid%r(n)))
         !write(std_out,*) 'x,energy', x,energy,ABS(x-energy),SQRT(ABS(x-energy))
         !call flush_unit(std_out)

         CALL backward_numerov(Grid,l,match-5,energy,Pot%rv,p1)
         rin=Gfirstderiv(Grid,match,p1)/p1(match)
         mmatch=match+1
         !WRITE(STD_OUT,*) 'match, rin ' ,match,rin
         !
         !  perform outward integration until match point -- it is assumed
         !   that projector functions proj are zero for r>r(match)
         !
         ! initialize u
         u=0
         u(2)=wfninit(0.d0,l,Pot%v0,Pot%v0p,energy,Grid%r(2))
         zeroval=0
         if (l==1) zeroval=2

         call forward_numerov(Grid,l,mmatch+5,energy&
&                 ,Pot%rv,zeroval,u,node1)
         !
         lcount=0
         DO ib=1,nbase
            IF (l==PAW%l(ib)) THEN
               lcount=lcount+1
               CALL inhomogeneous_numerov(Grid,l,mmatch+5,energy,&
&                  Pot%rv,PAW%otp(:,ib),y(:,lcount))
            ENDIF
         ENDDO
         !
         IF(lcount>0) THEN

            icount=0
            DO ib=1,nbase
               summ=0.d0
               IF (l==PAW%l(ib)) THEN
                  icount=icount+1
                  DO ic=1,nbase
                     IF (l==PAW%l(ic)) THEN
                        summ=summ+(PAW%dij(ib,ic)-energy*PAW%oij(ib,ic))*&
&                            overlap(Grid,PAW%otp(:,ic),u,1,irc)
                     ENDIF
                  ENDDO
                  b(icount)=-summ
               ENDIF
            ENDDO
            !
            icount=0
            DO ia=1,nbase
               IF (l==PAW%l(ia)) THEN
                  icount=icount+1
                  jcount=0
                  DO ib=1,nbase
                     IF (l==PAW%l(ib)) THEN
                        jcount=jcount+1
                        summ=0.d0
                        IF (ia.EQ.ib) summ=1.d0
                        DO ic=1,nbase
                           IF (l==PAW%l(ic)) THEN
                              summ=summ+(PAW%dij(ia,ic)-energy*PAW%oij(ia,ic))*&
&                             overlap(Grid,PAW%otp(:,ic),y(:,jcount),1,irc)
                           ENDIF
                        ENDDO
                        a(icount,jcount)=summ
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            !
            CALL linsol(a,b,lcount,lcount,lcount,nbase)

            wfn=0
            wfn(1:mmatch+5)=u(1:mmatch+5)
            icount=0
            DO ib=1,nbase
               IF(l==PAw%l(ib)) THEN
                  icount=icount+1
                  wfn(1:mmatch+5)=wfn(1:mmatch+5)+b(icount)*y(1:mmatch+5,icount)
               ENDIF
            ENDDO
         ENDIF

         rout=Gfirstderiv(Grid,match,wfn)/wfn(match)
         !WRITE(STD_OUT,'("node,match,rin,rout",2i8,1p,2e15.7)') node1,match,rin,rout
         !   -- estimate correction
         node1=0
         wfn(:)=wfn(:)/wfn(match)
         DO j=3,match
               IF (wfn(j)*wfn(j-1).LT.0.d0) node1=node1+1
         ENDDO
         !WRITE(STD_OUT,*) 'actual number of nodes', node1

!This test is obsolete: pseudo-WFs do not have to be orthogonal
!          IF (node1<node) THEN
!             ! too few nodes -- raise energy
!             emin=MAX(energy+err,emin)
!             energy=emax-(emax-energy)*ranx()
!             !WRITE(STD_OUT,*) 'too few nodes -- energy raised', energy,emin,emax
!          ELSEIF (node1>node) THEN
!             ! too many nodes -- lower energy
!             emax=MIN(energy-err,emax)
!             energy=emin+(energy-emin)*ranx()
!             !WRITE(STD_OUT,*) 'too many nodes -- energy lowered', energy,emin,emax
!             !do i=1,mmatch
!             !write(200+iter,'(1p,7e15.7)')Grid%r(i),wfn(i)
!             !enddo
!          ELSEIF (node1==node) THEN
            DO j=match,n
               wfn(j)=p1(j)/p1(match)
            ENDDO
            !  normalization
            scale=1.d0/sepnorm(Grid,PAW,n,l,wfn)
            dele=(rout-rin)*scale
            !WRITE(STD_OUT,*) 'dele,scale',dele,scale
            scale=SQRT(scale)
            wfn=scale*wfn

            !do i=1,n
            ! write(100+iter,'(1p,2E15.7)') Grid%r(i),wfn(i)
            !enddo

            x=ABS(dele)
            IF (x.LT.best) THEN
               energy0=energy
               best=x
            ENDIF
            IF (ABS(dele).LE.convre) EXIT
            energy=energy+dele
            !WRITE(STD_OUT,*) 'next energy' , energy,dele
            ! if energy is out of range, pick random energy in correct range
            IF (emin-energy.GT.convre.OR.energy-emax.GT.convre) THEN
               energy=emin+(emax-emin)*ranx()+err
            !   WRITE(STD_OUT,*) 'energy out of range -- reranged --', energy
            ENDIF
!          ENDIF
         iter=iter+1
         !WRITE(STD_OUT,*) 'Energy for next iteration ', iter,energy
         IF (iter.GT.niter) THEN
            WRITE(STD_OUT,*) 'no convergence in boundsep',l,dele,energy
            WRITE(STD_OUT,*) ' best guess of eig, dele = ',energy0,best
            STOP
         ENDIF
      ENDDO
      !
      ! normalize to unity within integration range
      !
      CALL filter(n,wfn,1.d-11)
      scale=1.d0/sepnorm(Grid,PAW,n,l,wfn)
      IF (scale.LE.0.d0) THEN
         WRITE(STD_OUT,*) 'warning -- negative norm for l=',l
         scale=-scale
         IF (scale.EQ.0.d0) scale=1.d0
      ENDIF
      scale=DSIGN(SQRT(scale),wfn(n-2))
      wfn(1:n)=wfn(1:n)*scale
      !WRITE(STD_OUT,*) 'exiting boundsep with energy ', l,energy
      DEALLOCATE(a,b,y,u,p1)
    END SUBROUTINE boundsep

    !*********************************************************
    !  program to to transform smooth wavefunction to all-electron
    !     wavefunction within Blochl's paw formalism
    !   otp == projector function
    !   odphi == ophi - otphi
    !*********************************************************
    SUBROUTINE PStoAE(Grid,PAW,nr,l,tpsi,psi)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(IN) :: PAW
      INTEGER, INTENT(IN) :: nr,l
      REAL(8), INTENT(IN) :: tpsi(:)
      REAL(8), INTENT(INOUT) :: psi(:)

      INTEGER :: n,ib,nbase,irc
      REAL(8) :: h,scale

      n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc

      IF (nr<irc) THEN
         WRITE(STD_OUT,*) 'error in PStoAE  -- nr < irc', nr,irc
         STOP
      ENDIF

      psi(1:nr)=tpsi(1:nr)
      DO ib=1,nbase
         IF (l==PAW%l(ib)) psi(1:nr)=psi(1:nr)+&
&             (PAW%ophi(1:nr,ib)-PAW%otphi(1:nr,ib))*&
&             overlap(Grid,PAW%otp(:,ib),tpsi,1,irc)
      ENDDO
    END SUBROUTINE PStoAE

    SUBROUTINE Set_PAW_MatrixElements(Grid,PAW,ifen)
  ! Calculate PAW matrix elements from reference configuration before
  !       unscreening/rescreening
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER, INTENT(IN) :: ifen

      INTEGER :: nbase,l,ib,ic,io
      REAL(8) :: x,y,occ
      TYPE(OrbitInfo), POINTER :: PSO
      REAL(8), allocatable :: wij(:,:)

      PAW%oij=0
      PAW%dij=0
      PAW%wij=0
      nbase=PAW%nbase
      PSO=>PAW%TOCCWFN
      ALLOCATE (wij(nbase,nbase))

      do ib=1,nbase
         do ic=1,nbase
            IF (PAW%l(ib)==PAW%l(ic)) THEN
               CALL dqij(Grid,PAW,ib,ic,PAW%oij(ib,ic))
               CALL altdtij(Grid,PAW,ib,ic,PAW%dij(ib,ic))
               CALL avij(Grid,PAW,ib,ic,x)
               PAW%dij(ib,ic)=PAW%dij(ib,ic)+x
             write(std_out,'(2i5, 1p,3e17.8)') ib,ic,PAW%oij(ib,ic),PAW%dij(ib,ic),x
             Endif
         enddo
      enddo

      wij=0
      Do io=1,PSO%norbit
         l=PSO%l(io)
         if (PSO%occ(io)>1.d-8.and..NOT.PSO%iscore(io)) then
             CALL calcwij(Grid,PAW,l,PSO%occ(io),PSO%wfn(:,io),wij)
         endif
      enddo

      do ib=1,nbase
         do ic=1,nbase
             PAW%wij(ib,ic)=wij(ib,ic)
        enddo
      enddo

      DEALLOCATE(wij)

      if (.not.Check_overlap_of_projectors(Grid,PAW,ifen)) then
         write(std_out,*) "The overlap operator has at least one negative eigenvalue!"
         write(std_out,*) "It might be not positive definite..."
         write(std_out,*) "Program is stopping."
         write(std_out,*) "This probably means that your projectors are too similar"
         write(std_out,*) "or that your PAW basis is incomplete."
         write(std_out,*) "Advice: try to change your input parameters (f.i. : reference energies)."
         stop
      end if

  END SUBROUTINE Set_PAW_MatrixElements

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  FUNCTION Check_overlap_of_projectors(Grid,PAW,ifen)
  !     function written by Francois Jollet and Marc Torrent 
  !       Approximate assessment of "completeness" of PAW basis set
  !       measured by the eigenvalues of the 1+O operator int
  !       the "basis" of the projectors     May 2019
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  FUNCTION Check_overlap_of_projectors(Grid,PAW,ifen)
  ! Check that the overlap operator is positive definite
  !  in the "basis" of projectors
      LOGICAL :: Check_overlap_of_projectors
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      INTEGER, INTENT(IN) :: ifen

      REAL(8),PARAMETER :: tol=1.d-10

      INTEGER :: n,nbase,i,j,k,m,info
      REAL(8) :: ovlp_ik,ovlp_jm 
      REAL(8),allocatable :: ovlp(:,:),dum(:)
!      REAL(8),allocatable :: wr(:),wi(:),work(:),vl(:,:),vr(:,:)
      REAL(8), allocatable :: eig(:),work(:)

      Check_overlap_of_projectors=.true.

      n=PAW%irc ; nbase=PAW%nbase
      allocate(ovlp(nbase,nbase))
      ovlp(:,:)=0.d0

      allocate(dum(n))
      do k=1,nbase
        do m=1,nbase
          if (PAW%l(k)==PAW%l(m)) then
            dum(1:n)=PAW%otp(1:n,k)*PAW%otp(1:n,m)
            ovlp(m,k)=integrator(Grid,dum(1:n),1,n)

            do i=1,nbase
              if (PAW%l(i)==PAW%l(k)) then
                dum(1:n)=PAW%otp(1:n,i)*PAW%otp(1:n,k)
                ovlp_ik=integrator(Grid,dum(1:n),1,n)

                do j=1,nbase
                  if (PAW%l(j)==PAW%l(m)) then
                    dum(1:n)=PAW%otp(1:n,j)*PAW%otp(1:n,m)
                    ovlp_jm=integrator(Grid,dum(1:n),1,n)

                    ovlp(m,k)=ovlp(m,k)+ovlp_ik*PAW%oij(i,j)*ovlp_jm

                  end if ! l_j=l_m
                end do   ! Loop j
              end if     ! l_i=l_k
            end do       ! Loop i
          end if         ! l_k=l_m
        end do           ! Loop m
      end do             ! Loop k
      deallocate(dum)

!      allocate(wr(nbase),wi(nbase),work(4*nbase))
!      allocate(vl(nbase,nbase),vr(nbase,nbase))
!      call dgeev('N','N',nbase,ovlp,nbase,wr,wi,vl,nbase,vr,nbase,work,4*nbase,info)

      allocate(eig(nbase),work(4*nbase))
      call DSYEV('N','U',nbase,ovlp,nbase,eig,work,4*nbase,info)
      write(std_out,'(" Completed diagonalization of ovlp with info = ", i8)')info
      write(ifen,'(" Completed diagonalization of ovlp with info = ", i8)')info
       if (info==0) then
         write(std_out,*) " "
         write(std_out,*) "Eigenvalues of overlap operator (in the basis of projectors):"
         write(ifen,*) " "
         write(ifen,*) "Eigenvalues of overlap operator (in the basis of projectors):"
         do i=1,nbase
            write(std_out,'( i5,5x,1p,1e17.8)') i,eig(i)
            write(ifen,'( i5,5x,1p,1e17.8)') i,eig(i)
            if (eig(i)<tol) Check_overlap_of_projectors=.false.
         enddo
        else
          write(std_out,*) 'Stopping due to failure of ovlp diagonalization'
          write(ifen,*) 'Stopping due to failure of ovlp diagonalization'
          stop
       endif       
      write(std_out,*) " "
      write(ifen,*) " "

!      do i=1,nbase
!        write(std_out,*) i,wr(i)
!        if (wr(i)<tol) Check_overlap_of_projectors=.false.
!      end do

      deallocate(eig,work)
      deallocate(ovlp)

   END FUNCTION Check_overlap_of_projectors

    !************************************************************************
    !  program to calculate logerivatives of paw wavefunctions
    !   and to compare them with all electron wavefunctions
    !  optionally, wavefunctions are written to a file
    !  Assumes prior call to SUBROUTINE Set_PAW_MatrixElements(Grid,PAW)
    ! 7/2018 program modified by Casey Brock to output atan(logderiv)
    !      values, with pi jumps taken into account
    !   
    !************************************************************************
    SUBROUTINE logderiv(Grid,Pot,PAW)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW

      TYPE(PotentialInfo) :: PS
      INTEGER :: n,l,ie,ke,nbase,ib,ic,nr,i,nodes,mbase,irc,lng
      REAL(8), PARAMETER :: e0=-5.d0 ! Where to print WFn in file
      REAL(8) :: de,h,x,dwdr,dcwdr,scale,energy
      REAL(8), ALLOCATABLE :: psi(:),tpsi(:),ttpsi(:)
      CHARACTER(4)  :: flnm
      REAL(8) :: y_ae, x_ae, y_ps, x_ps 
      REAL(8) :: cumshift_ae, cumshift_ps
      REAL(8) :: atan_ae, atan_ps, atan_prev_ae, atan_prev_ps 
      REAL(8) :: slope_ae, slope_ps, slope_prev_ae, slope_prev_ps
      REAL(8) :: atan_shifted_ae, atan_shifted_ps
      INTEGER :: cnt_ae, cnt_ps


      n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc;  nr=irc+10

      ALLOCATE(psi(nr),tpsi(nr),ttpsi(nr),PS%rv(n),PS%vtau(n),PS%rvx(n),stat=ie)
      IF (ie/=0) THEN
         WRITE(STD_OUT,*) 'Error in logderiv allocation',n,ie
         STOP
      ENDIF

      ! load  PS
      PS%rv=PAW%rveff ; PS%nz=0.d0 ; PS%rvx=PAW%rveff
      if(needvtau) PS%vtau=PAW%tvtau
      call zeropot(Grid,PS%rv,PS%v0,PS%v0p)
      !
      !   calculate logderivatives at irc
      !
      WRITE(STD_OUT,*) 'calculating log derivatives at irc',Grid%r(irc)
      !
      de=(maxlogderiv-minlogderiv)/dble(nlogderiv-1)
      ke=1+anint((e0-minlogderiv)/de)

      mbase=nbase
      DO l=0,PAW%lmax+1
         CALL mkname(l,flnm)
         OPEN(56,file='logderiv.'//TRIM(flnm),form='formatted')

! initialize variables for atan calculation
         cnt_ae = 1
         slope_ae = 0.d0
         cumshift_ae = 0.d0 
         cnt_ps = 1
         slope_ps = 0.d0
         cumshift_ps = 0.d0 

         DO ie=1,nlogderiv
            energy=minlogderiv+de*(ie-1)
            psi=0;tpsi=0;ttpsi=0
            if (scalarrelativistic) then
               CALL unboundsr(Grid,Pot,nr,l,energy,psi,nodes)
            else if (needvtau) then
               CALL unboundked(Grid,Pot,nr,l,energy,psi,nodes)
            elseif (TRIM(PAW%exctype)=='HF') then
               CALL HFunocc(Grid,PAW%OCCWFN,l,energy,Pot%rv,Pot%v0,Pot%v0p,&
&                       psi,lng)
            else
               CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,energy,psi,nodes)
            endif
            if(needvtau) then
                CALL unboundsepked(Grid,PS,PAW,nr,l,energy,tpsi,nodes)    
            else
                CALL unboundsep(Grid,PS,PAW,nr,l,energy,tpsi,nodes)
            endif
            CALL PStoAE(Grid,PAW,nr,l,tpsi,ttpsi)
            !
!            dwdr=Gfirstderiv(Grid,irc,psi)/psi(irc)
!            dcwdr=Gfirstderiv(Grid,irc,ttpsi)/ttpsi(irc)
             !!! CALCULATE ATANS 
             y_ae = Gfirstderiv(Grid,irc,psi)
             x_ae = psi(irc)
             y_ps = Gfirstderiv(Grid,irc,ttpsi)
             x_ps = ttpsi(irc)
 
             dwdr = y_ae/x_ae
             dcwdr = y_ps/x_ps
             
! calculate atan and 
! phase shift by n*pi if necessary so the atan curve is continous
            atan_ae = atan2(y_ae, x_ae)
            IF (ie>1) THEN
               slope_ae = atan_ae - atan_prev_ae
               CALL phase_unwrap(ie, y_ae, x_ae, atan_ae, atan_prev_ae, slope_ae, slope_prev_ae, cumshift_ae, cnt_ae)
            ENDIF
            atan_shifted_ae = atan_ae + cumshift_ae
            slope_prev_ae = slope_ae
            atan_prev_ae = atan_ae
            atan_ps = atan2(y_ps, x_ps)
            IF (ie>1) THEN
               slope_ps = atan_ps - atan_prev_ps
               CALL phase_unwrap(ie, y_ps, x_ps, atan_ps, atan_prev_ps, slope_ps, slope_prev_ps, cumshift_ps, cnt_ps)
            ENDIF
            atan_shifted_ps = atan_ps + cumshift_ps
            slope_prev_ps = slope_ps
            atan_prev_ps = atan_ps

            WRITE(56,'(1p,5e12.4)') energy,dwdr,dcwdr,atan_shifted_ae,&
&                 atan_shifted_ps
            IF (ie.EQ.ke) THEN
               mbase=mbase+1
               CALL mkname(mbase,flnm)
               OPEN(57,file='wfn'//TRIM(flnm),form='formatted')
               WRITE(57,*) '# l=',l,'energy=',energy
               !
               ! form converted wavefunction and rescale exact wavefunction
               !
               scale=ttpsi(irc)/psi(irc)
               DO i=1,nr
                  WRITE(57, '(1p,5e12.4)') Grid%r(i),tpsi(i),ttpsi(i),psi(i)*scale
               ENDDO
               CLOSE(57)
            ENDIF
         ENDDO !ie
         CLOSE(56)

      ENDDO !l

      DEALLOCATE(psi,tpsi,PS%rv,PS%vtau)
    END SUBROUTINE logderiv



    !  Assumes prior call to SUBROUTINE calculate_tvtau
    !  which now fills PAW%tden and PAW%ttau
    SUBROUTINE FindVlocfromVeff(Grid,Orbit,PAW)
      TYPE(GridInfo), INTENT(INOUT) :: Grid
      TYPE(OrbitInfo), INTENT(IN) :: Orbit
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW

      REAL(8), POINTER  :: r(:)
      REAL(8) :: h,qeff,tq,rat,q00,ecoul,etxc,eexc,occ,fac
      INTEGER :: n,i,irc,io,nbase,ib,ic
      REAL(8), allocatable :: d(:),dv(:),dvx(:),v(:),vv(:),vxB(:),vxK(:)
      REAL(8), allocatable :: t(:),vt(:),vthat(:)

      CALL FillHat(Grid,PAW)

      if (Vlocalindex == SETVLOC) then
         write(std_out,*) 'Vloc == VlocCoef*shapefunc  '
         return
      endif

      n=Grid%n;irc=PAW%irc
      nbase=PAW%nbase
      h=Grid%h ; r=>Grid%r
      irc=max(PAW%irc,PAW%irc_shap,PAW%irc_vloc,PAW%irc_core)

! Recalculate den and tau      
      PAW%den=0.d0;PAW%tden=0.d0
      PAW%valetau=0.d0;PAW%tvaletau=0.d0
      allocate(d(n))
      do io=1,PAW%OCCWFN%norbit
         if (.not.PAW%OCCwfn%iscore(io)) then
            occ=PAW%OCCwfn%occ(io)
            PAW%den=PAW%den+occ*PAW%OCCwfn%wfn(:,io)**2
            PAW%tden=PAW%tden+occ*PAW%TOCCwfn%wfn(:,io)**2
            CALL taufromwfn(d,Grid,PAW%OCCwfn%wfn(:,io),PAW%OCCwfn%l(io))
            PAW%valetau=PAW%valetau+occ*d
            CALL taufromwfn(d,Grid,PAW%TOCCwfn%wfn(:,io),PAW%TOCCwfn%l(io))
            PAW%tvaletau=PAW%tvaletau+occ*d
          endif
      enddo
      deallocate(d)

      allocate(d(n),dv(n),dvx(n),vxB(n),vxK(n),STAT=i)
      allocate(t(n),vt(n),vthat(n))       
      if (i/=0) stop 'Error (1) in allocating  arrays in findvlocfromveff'

      d=PAW%tden+PAW%tcore
      call poisson(Grid,tq,d,dv,ecoul)

      d=PAW%core-PAW%tcore
      q00=0.5d0*PAW%AErefrv(1)+integrator(Grid,d)
      write(std_out,*) 'nucleus and core q00 ', q00
      do ib=1,nbase
         do ic=1,nbase
            q00=q00+PAW%wij(ib,ic)*PAW%oij(ib,ic)
         enddo
      enddo
      write(std_out,*) 'Complete q00 ', q00

      dv=dv+q00*PAW%hatpot
      If (TRIM(PAW%exctype)=='EXX') THEN
        CALL EXX_pseudoVx(Grid,Orbit,PAW,dvx)
        PAW%trvx=dvx
        dv=dv+dvx
      ELSEIf (TRIM(PAW%exctype)=='EXXKLI') THEN
        write(std_out,*) 'before EXX_pseudo'; call flush_unit(std_out)
        CALL EXXKLI_pseudoVx(Grid,PAW,dvx)
        PAW%trvx=dvx
        dv=dv+dvx
      ELSEIf (TRIM(PAW%exctype)=='HF') THEN
        dvx=0.d0
        PAW%trvx=0.d0
      ELSE
        d=PAW%tden+PAW%tcore
        t=PAW%tcoretau+PAW%tvaletau
        CALL exch(Grid,d,dvx,etxc,eexc,tau=t,vtau=vt)
        PAW%trvx=dvx
        if (needvtau) PAW%tvtau=vt
        dv=dv+dvx
      endif

      rat=0.d0
      DO i=2,n
         PAW%vloc(i)=(PAW%rveff(i)-dv(i))/r(i)
         IF (i>=irc) rat=rat+ABS(PAW%vloc(i))
      ENDDO
      WRITE(STD_OUT,*) 'Error in vloc -- ',rat
      call extrapolate(Grid,PAW%vloc)

!!! Check rveff and unscreening potential
      open(1001,file='Vloccheck',form='formatted')
      write(1001,'(a)') '# r     rveff           vscreen         vloc      vhat'
      do i=1,n
         write(1001,'(1p,10e15.7)') Grid%r(i),PAW%rveff(i),dv(i),PAW%vloc(i),q00*PAW%hatpot(i)
      enddo    
      close(1001)

      PAW%vloc(irc:n)=0.d0


!!!!!!!!!!!This part does not work for HF!!!!
!     Construct ionic local potential for abinit from screened pseudopotential
!     in addition to ionic unscreening include hathat density
!     in exchange-correlation functional
!     also generate version without hathat density in vxc

      PAW%abinitvloc=0.d0
      PAW%abinitnohat=0.d0
      allocate(v(n),vv(n),STAT=i)
      if (i/=0) stop 'Error (2) in allocating  arrays in findvlocfromveff'

      d=PAW%den-PAW%tden
      tq=integrator(Grid,d,1,irc)
      write(std_out,*) ' abinit tq = ', tq

!     Compute VH(tDEN+hatDEN)
      d=PAW%tden+tq*PAW%hatden
      CALL poisson(Grid,q00,d,v,rat)

!     Compute Vxc(tcore+tDEN)
      d=PAW%tcore+PAW%tden
        t=PAW%tcoretau+PAW%tvaletau
        CALL exch(Grid,d,v,etxc,eexc,tau=t,vtau=vt)

!     Compute Vxc(tcore+tDEN+hatDEN)
      d=PAW%tcore+PAW%tden+tq*PAW%hatden
        t=PAW%tcoretau+PAW%tvaletau
        CALL exch(Grid,d,vv,etxc,eexc,tau=t,vtau=vthat)

!     Store Vxc(tcore+tDEN)-Vxc(tcore+tDEN+hatDEN)
      do i=2,n
        PAW%abinitvloc(i)=(v(i)-vv(i))/Grid%r(i)
      enddo
      call extrapolate(Grid,PAW%abinitvloc)

!     Compute VH(tcore+(Nc-tNc-Z)*hatDEN)
!     (For consistency with ABINIT, need to integrate Poisson
!      equation up to coretailpoints (not the whole mesh))
      d=PAW%core-PAW%tcore
      qeff=0.5d0*PAW%AErefrv(1)+integrator(Grid,d,1,PAW%irc_core)
      d=PAW%tcore+qeff*PAW%hatden
      Grid%n=PAW%coretailpoints
      CALL poisson_marc(Grid,q00,d(1:PAW%coretailpoints),vv(1:PAW%coretailpoints),rat)
      Grid%n=n

!     Compute ionic potential following Bloechl formalism
      do i=2,PAW%coretailpoints
        PAW%abinitnohat(i)=PAW%vloc(i)+vv(i)/r(i)  ! in Rydberg units
      enddo
      call extrapolate(Grid,PAW%abinitnohat)
      PAW%abinitnohat(PAW%coretailpoints+1:n)=PAW%vloc(PAW%coretailpoints+1:n) &
&                     +vv(PAW%coretailpoints)/Grid%r(PAW%coretailpoints+1:n)

!     Compute ionic potential following Kresse formalism
      do i=1,n
        PAW%abinitvloc(i)=PAW%abinitvloc(i)+PAW%abinitnohat(i)  ! in Rydberg units
      enddo

! check if PAW%tcore+tq*PAW%hatden is positive      
      PAW%poscorenhat=.true.
      PAW%nhatv=tq*PAW%hatden
      do i=1,irc
         if ((PAW%tcore(i)+PAW%nhatv(i))<-1.d-8) then
          PAW%poscorenhat=.false.
          exit       
         endif
      enddo   

      open(123,file='compare.abinit', form='formatted')
      write(123,*) '#r    veff      vloc(B)        vloc(K)        vloc(Ab) vtau(Kresse)'
      do i=2,n
         write(123,'(1p,6e16.7)') r(i),PAW%rveff(i)/r(i),&
&             PAW%vloc(i),PAW%abinitvloc(i),PAW%abinitnohat(i),vthat(i)
      enddo
      close(123)

!!!!   Unscreen Kresse pseudopotential
      d=PAW%den-PAW%tden
      tq=integrator(Grid,d,1,irc)
      write(std_out,*) ' abinit tq = ', tq

!     Compute VH(tDEN+hatDEN)
      d=PAW%tden+tq*PAW%hatden
      CALL poisson(Grid,q00,d,v,rat)

!     Compute Vxc(tcore+tDEN+hatDEN)
      d=PAW%tcore+PAW%tden+tq*PAW%hatden
        t=PAW%tcoretau+PAW%tvaletau
        CALL exch(Grid,d,vv,etxc,eexc,tau=t,vtau=vthat)

     PAW%Kunscreen=PAW%Krveff-v-vv
      open(123,file='Kressepot.dat',form='formatted')
      write(123,*)'#r   AErefrv      Krveff      Har    vxc    Unscreened  Ktvt'
      do i=1,n
        write(123,'(1p,8e16.7)') r(i),PAW%AErefrv(i),PAW%Krveff(i), &
&        v(i),vv(i),PAW%Kunscreen(i),vthat(i),PAW%Ktvtau(i)
     enddo
     close(123)     

     write(std_out,*) 'Resmoothing tvtau and Ktvtau'
     PAW%tvtau=vt
     PAW%Ktvtau=vthat
     Call Smoothtvtau(Grid,PAW,0.001d0)
      open(123,file='Resmoothvtau.dat',form='formatted')
      write(123,*)'#r    Blochl       Blochlsmooth         Kresse     Kressesmooth'
      do i=1,n
        write(123,'(1p,8e16.7)') r(i),vt(i),PAW%tvtau(i),vthat(i),PAW%Ktvtau(i)
      enddo
      close(123)     


!  Recalculate  PAW%abinitvloc (For Kresse)  and PAW%abinitnohat  (For Blochl)    
!!!!  Redo potentials
!  \Delta Qval
      d=PAW%den-PAW%tden
      tq=integrator(Grid,d,1,irc)
      write(std_out,*) ' Delta Qval = ', tq

!     Compute VH(tDEN+hatDEN)
      d=PAW%tden+tq*PAW%hatden
      CALL poisson(Grid,q00,d,v,rat)
      write(std_out,*) ' Completed Poisson with q00 = ', q00
      write(std_out,*) ' Completed Poisson with v(n) = ', v(n)

!  Compute Blochl exc      
        d=PAW%tden+PAW%tcore
        t=PAW%tcoretau+PAW%tvaletau
        CALL exch(Grid,d,vxB,etxc,eexc,tau=t,vtau=vt)

!     Compute Kresse   exc Vxc(tcore+tDEN+hatDEN)
      d=PAW%tcore+PAW%tden+tq*PAW%hatden
        t=PAW%tcoretau+PAW%tvaletau
        CALL exch(Grid,d,vxK,etxc,eexc,tau=t,vtau=vthat)

      PAW%abinitvloc=0.d0; PAW%abinitnohat=0
      PAW%abinitnohat(2:n)=(PAW%rveff(2:n)-v(2:n)-vxB(2:n))/r(2:n)  
      call extrapolate(Grid,PAW%abinitnohat)
      PAW%abinitvloc(2:n)=(PAW%rveff(2:n)-v(2:n)-vxK(2:n))/r(2:n)  
      call extrapolate(Grid,PAW%abinitvloc)

      open(123,file='newcompare.abinit', form='formatted')
      write(123,*) '#r    veff    rvH        Blochl Pseudo   vxc  vtau      Kresse Pseudo   vxc  vtau '
      do i=2,n
         write(123,'(1p,9e16.7)') r(i),PAW%rveff(i)/r(i), v(i),&
&             PAW%abinitnohat(i),vxB(i),vt(i),PAW%abinitvloc(i),vxK(i),vthat(i)
      enddo
      close(123)


! Reassess poscorenhat      
! check if PAW%tcore+PAW%tden+tq*PAW%hatden is positive      
      PAW%poscorenhat=.true.
      open(123,file='checkpositivity',form='formatted')
      do i=1,n
         rat=PAW%tcore(i)+PAW%tden(i)+tq*PAW%hatden(i)
         write(123,'(1p, 5e16.7)') r(i),PAW%tcore(i),PAW%tden(i),tq*PAW%hatden(i),rat 
         if ((PAW%tcore(i)+PAW%tden(i)+tq*PAW%hatden(i))<-1.d-8) then
          PAW%poscorenhat=.false.
          exit       
         endif
      enddo   
      close(123)


      deallocate(v,vv,vxB,vxK)
      deallocate(d,dv,dvx,vt,vthat)

    END SUBROUTINE FindVlocfromVeff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   SCFPAW
!     It is assumed that the new configuration involves only occupation
!         changes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE SCFPAW(Grid,PAW)
     Type(GridInfo), INTENT(IN) :: Grid
     Type(PseudoInfo), INTENT(INOUT) :: PAW

     INTEGER :: i,j,k,l,n,nvorbit,nbase,io,ib,jo,jb,ip,nfix
     INTEGER :: loop
     INTEGER :: mxloop=1000
     REAL(8) :: err0=1.d-7,mix=0.5d0
     TYPE(OrbitInfo), POINTER :: AEO,PSO
     INTEGER :: firsttime=0,norbit_mod,np(5)
     REAL(8) :: xocc,err
     LOGICAL :: success
     CHARACTER(4) :: stuff
     INTEGER, ALLOCATABLE :: tmap(:)
     INTEGER, ALLOCATABLE :: orbit_mod_l(:),orbit_mod_n(:),orbit_mod_k(:)
     REAL(8), ALLOCATABLE :: orbit_mod_occ(:)

     AEO=>PAW%OCCWFN;     PSO=>PAW%TOCCWFN
     WRITE(STD_OUT,*) 'Current occupancies:'
     IF (.NOT.diracrelativistic) THEN
       WRITE(STD_OUT,*) ' n l   occupancy        energy    '
     ELSE
       WRITE(STD_OUT,*) ' n l kap   occupancy        energy    '
     ENDIF
     DO io=1,PSO%norbit
         IF (.NOT.PSO%iscore(io)) THEN
             IF (.NOT.diracrelativistic) THEN
               WRITE(STD_OUT,'(i2,1x,i2,4x,1p,2e15.7)') PSO%np(io),&
&                   PSO%l(io),PSO%occ(io), PSO%eig(io)
             ELSE
               WRITE(STD_OUT,'(i2,1x,i2,2x,i2,4x,1p,2e15.7)') PSO%np(io),&
&                   PSO%l(io),PSO%kappa(io),PSO%occ(io), PSO%eig(io)
             ENDIF
             if (firsttime==0) then
                ib=PAW%valencemap(io)
                write(std_out,*) 'Setting pseudo orbital ', io,ib
                PSO%wfn(:,io)=PAW%tphi(:,ib)
             endif
         ENDIF
     ENDDO
     firsttime=firsttime+1

     !Read modified occupations from standard input
     np(1)=AEO%nps;np(2)=AEO%npp;np(3)=AEO%npd;np(4)=AEO%npf;np(5)=AEO%npg
     CALL input_dataset_read_occ(norbit_mod,orbit_mod_l,orbit_mod_n,orbit_mod_k,&
&                                orbit_mod_occ,np,diracrelativistic)

     DO jo=1,norbit_mod
       nfix=-100
       DO io=1,PSO%norbit
          IF (orbit_mod_n(jo)==PSO%np(io).AND.orbit_mod_l(jo)==PSO%l(io).AND.&
&             ((.NOT.diracrelativistic).OR.orbit_mod_k(jo)==PSO%kappa(io))) THEN
             IF (PSO%iscore(io)) THEN
                 write(std_out,*) 'Core orbitals cannot be changed',&
                 &  orbit_mod_n(jo),orbit_mod_l(jo),orbit_mod_occ(jo)
                 stop
             endif
             nfix=io
             EXIT
          ENDIF
       ENDDO
       IF (nfix.LE.0) THEN
          WRITE(STD_OUT,*) 'error in occupations -- ip,l,xocc', &
&             orbit_mod_n(jo),orbit_mod_l(jo),orbit_mod_occ(jo),nfix
          STOP
       ENDIF
       PSO%occ(nfix)=orbit_mod_occ(jo)
       AEO%occ(nfix)=orbit_mod_occ(jo)
     ENDDO
     DEALLOCATE(orbit_mod_l)
     DEALLOCATE(orbit_mod_n)
     DEALLOCATE(orbit_mod_k)
     DEALLOCATE(orbit_mod_occ)

     WRITE(STD_OUT,*) 'New configuration:'
     DO io=1,PSO%norbit
        If (.NOT.PSO%iscore(io)) then
           IF (.NOT.diracrelativistic) THEN
             WRITE(STD_OUT,'(i2,1x,i2,4x,1p,2e15.7)') PSO%np(io),&
&                  PSO%l(io),PSO%occ(io)
           ELSE
             WRITE(STD_OUT,'(i2,1x,i2,2x,i2,4x,1p,2e15.7)') PSO%np(io),&
&                  PSO%l(io),PSO%kappa(io),PSO%occ(io)
           END IF
        Endif
     Enddo

     success=.FALSE.
     do loop=1,mxloop
        if (TRIM(PSO%exctype)=='HF') then
           call PAWIter_HF(Grid,PAW,mix*0.01,err0,err,success)
        else
           call PAWIter_LDA(Grid,PAW,mix,err0,err,success)
        endif
        write(std_out,*)  '--Results for Iter -- ', loop
        IF (.NOT.diracrelativistic) THEN
          write(std_out,*)  '  n   l   occupancy          energy    '
          do io=1,PSO%norbit
             If (.NOT.PSO%iscore(io)) then
                WRITE(STD_OUT,'(i2,1x,i2,4x,1p,2e15.7)') PSO%np(io),&
&                      PSO%l(io),PSO%occ(io),PSO%eig(io)
             endif
          enddo
        ELSE
          write(std_out,*)  '  n   l kap   occupancy          energy    '
          do io=1,PSO%norbit
             If (.NOT.PSO%iscore(io)) then
                WRITE(STD_OUT,'(i2,1x,i2,2x,i2,4x,1p,2e15.7)') PSO%np(io),&
&                      PSO%l(io),PSO%kappa(io),PSO%occ(io),PSO%eig(io)
             endif
          enddo
        END IF
        IF (success) then
           WRITE(STD_OUT,*) ' PS wfn iteration converged ', loop
           write(std_out,*)  '--Results for Iter -- ', loop
           IF (.NOT.diracrelativistic) THEN
             write(std_out,*)  '  n   l   occupancy          energy    '
             do io=1,PSO%norbit
                If (.NOT.PSO%iscore(io)) then
                   WRITE(STD_OUT,'(i2,1x,i2,4x,1p,2e15.7)') PSO%np(io),&
&                         PSO%l(io),PSO%occ(io),PSO%eig(io)
                endif
             enddo
           ELSE
             write(std_out,*)  '  n   l kap   occupancy          energy    '
             do io=1,PSO%norbit
                If (.NOT.PSO%iscore(io)) then
                   WRITE(STD_OUT,'(i2,1x,i2,2x,i2,4x,1p,2e15.7)') PSO%np(io),&
&                         PSO%l(io),PSO%kappa(io),PSO%occ(io),PSO%eig(io)
                endif
             enddo
           ENDIF
           exit
        ENDIF
     enddo

     Call mkname(firsttime,stuff)
     allocate(tmap(PSO%norbit))
     ip=0; tmap=0
     do io=1,PSO%norbit
        if (.not.PSO%iscore(io)) then
           ip=ip+1
           tmap(ip)=io
           call PStoAE(Grid,PAW,Grid%n,PSO%l(io),PSO%wfn(:,io),&
&                 PAW%OCCWFN%wfn(:,io))
        endif
      enddo
      OPEN(unit=1001,file='PAWwfn.'//TRIM(stuff),form='formatted')
      do i=1,Grid%n
         write(1001,'(1p,60e15.7)') Grid%r(i),(PSO%wfn(i,tmap(k)),&
&                   PAW%OCCWFN%wfn(i,tmap(k)),k=1,ip)
      enddo
      CLOSE(1001)

     deallocate(tmap)

   END SUBROUTINE SCFPAW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PAWIter_LDA(Grid,PAW,err0,err,success)
!     On input PAW%TOCCWFN  contains initial guess of smooth wfn's
!        On output the guess is updated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PAWIter_LDA(Grid,PAW,mix,err0,err,success)
      Type(GridInfo), INTENT(IN) :: Grid
      Type(PseudoInfo), INTENT(INOUT) :: PAW
      REAL(8) , INTENT(IN) :: mix,err0
      REAL(8) , INTENT(OUT) :: err
      LOGICAL , INTENT(OUT) :: success

      INTEGER :: i,j,k,l,n,io,jo,ib,jb,kb,lb,irc,nbase,nocc
      INTEGER, allocatable :: tmap(:)
      REAL(8) , ALLOCATABLE :: arg(:),rhs(:),rv(:),aden(:),v1(:),v2(:),o(:,:)
      REAL(8) , ALLOCATABLE :: t(:),vt(:)
      Type(OrbitInfo), POINTER :: PSO
      Type(OrbitInfo) :: tmpOrbit
      REAL(8) :: occ,x,y,q,v0term,en,val,mix1
      INTEGER :: fcount=0

      success=.false.
      n=Grid%n; nbase=PAW%nbase;  irc=PAW%irc
      ALLOCATE(arg(n),rhs(n),rv(n),aden(n),v1(n),v2(n),&
&           tmap(PAW%OCCWFN%norbit),o(PAW%OCCWFN%norbit,nbase))
      ALLOCATE(t(n),vt(n))

      PSO=>PAW%TOCCWFN
      ! orthonormalize
      nocc=0
      do io=1,PSO%norbit
         if (.not.PSO%iscore(io).and.io>1) then
            do jo=1,io-1
               if (PSO%l(jo)==PSO%l(io).and..not.PSO%iscore(jo)) then
                     call genOrthog(Grid,PAW,n,PSO%l(io),&
&                        PSO%wfn(:,io),PSO%wfn(:,jo))
                     write(std_out,*) 'orthog', io,jo
               endif
            enddo
         x=genoverlap(Grid,PAW,n,PSO%l(io),PSO%wfn(:,io),PSO%wfn(:,io))
         PSO%wfn(:,io)=PSO%wfn(:,io)/sqrt(x)
         CALL ADJUSTSIGN(PSO%wfn(:,io),3)
         write(std_out,*) 'normalize ', io, x
         nocc=nocc+1
         tmap(nocc)=io
         endif
      enddo

     !  prepare overlap terms
      o=0
      do k=1,nocc
         io=tmap(k); l=PSO%l(io)
         do ib=1,PAW%nbase
            if (l==PAW%l(ib)) then
               o(io,ib)=overlap(Grid,PSO%wfn(:,io),PAW%otp(:,ib),1,irc)
                write(std_out,'("<p|psi> ", 2i5,1p,e15.7)') io,ib,o(io,ib)
            endif
         enddo
      enddo

      Call CopyOrbit(PSO,tmpOrbit)

      rv=PAW%rtVf; PAW%tkin=0;   PAW%wij=0
      PAW%tden=0;rhs=0;aden=0;x=0
      do k=1,nocc
         io=tmap(k) ; l=PSO%l(io)
            occ=PSO%occ(io)
            PAW%tden=PAW%tden+occ*(PSO%wfn(:,io))**2
            call kinetic_ij(Grid,PSO%wfn(:,io),PSO%wfn(:,io),l,y)
            PAW%tkin= PAW%tkin + occ*y
            PSO%eig(io)=y
            write(std_out,*) 'Kinetic ',io,y
            do ib=1,PAW%nbase
               do jb=1,PAW%nbase
                  if (PAW%l(ib)==l.and.PAW%l(jb)==l) then
                      x=x+occ*o(io,ib)*o(io,jb)*PAW%mLij(ib,jb,1)
                      PAW%wij(ib,jb)=PAW%wij(ib,jb)+occ*o(io,ib)*o(io,jb)
                   endif
               enddo
            enddo
      enddo

        aden=PAW%tden+x*PAW%g(:,1)
        arg=0; arg(2:n)=(aden(2:n)*PAW%hatpot(2:n))/Grid%r(2:n)
        v0term=integrator(Grid,arg)
        write(std_out,*) 'v0term, aden', v0term,integrator(Grid,aden)

        call poisson(Grid,q,aden,rhs,x)
        write(std_out,*) 'PAW poisson ', q,x
        PAW%tvale=x
        arg=0; arg(2:n)=PAW%rtVf(2:n)/Grid%r(2:n)
        PAW%tion=overlap(Grid,arg,PAW%tden)
        rv=rv+rhs    ! vion + valence-Hartree
        arg=PAW%tden+PAW%tcore
        t=PAW%tcoretau+PAW%tvaletau
        call exch(Grid,arg,rhs,x,y,tau=t,vtau=vt)
        PAW%txc=y
        rv=rv+rhs    ! + vxc

        do k=1,nocc
           io=tmap(k);l=PSO%l(io)
           arg=rv*(PSO%wfn(:,io)**2);
           arg(1)=0; arg(2:n)=arg(2:n)/Grid%r(2:n)
           x=integrator(Grid,arg)
           write(std_out,*) ' potential term ',io, x
           PSO%eig(io)=PSO%eig(io)+x
        enddo

        PAW%dij=0; PAW%Ea=0
        do ib=1,PAW%nbase
           do jb=1,PAW%nbase
              if (PAW%l(ib)==PAW%l(jb)) then
                 PAW%dij(ib,jb)=PAW%dij(ib,jb) + PAW%Kij(ib,jb) + &
&                        PAW%VFij(ib,jb)+PAW%mLij(ib,jb,1)*v0term
                 PAW%Ea=PAW%Ea + PAW%wij(ib,jb)*(PAW%Kij(ib,jb) + &
&                        PAW%VFij(ib,jb))
                 x=0;    !  accumulate Hartree term
                 do kb=1,PAW%nbase
                    do lb=1,PAW%nbase
                       if (PAW%l(kb)==PAW%l(lb)) then
                          x=x+PAW%wij(kb,lb)*PAW%DR(ib,jb,kb,lb,1)
                       endif
                    enddo
                 enddo
                 PAW%dij(ib,jb)=PAW%dij(ib,jb)+x
                 PAW%Ea=PAW%Ea + 0.5d0*PAW%wij(ib,jb)*x
              endif
           enddo
        enddo

      !  exchange-correlation part
      arg=PAW%core;   rhs=PAW%tcore
      do ib=1,PAW%nbase
         do jb=1,PAW%nbase
            if (PAW%l(ib)==PAW%l(jb)) then
               arg=arg+PAW%wij(ib,jb)*PAW%ophi(:,ib)*PAW%ophi(:,jb)
               rhs=rhs+PAW%wij(ib,jb)*PAW%otphi(:,ib)*PAW%otphi(:,jb)
            endif
         enddo
      enddo

      Write(STD_OUT,*) 'Before EXC ', PAW%Ea
      irc=PAW%irc
         t=PAW%coretau+PAW%valetau
         call exch(Grid,arg,v1,x,y,fin=irc,tau=t,vtau=vt)
      PAW%Ea=PAW%Ea+y  ; write(std_out,*) 'AE EXC ' ,y
         t=PAW%tcoretau+PAW%tvaletau
         call exch(Grid,arg,v2,x,y,fin=irc,tau=t,vtau=vt)
      PAW%Ea=PAW%Ea-y  ; write(std_out,*) 'PS EXC ' ,y

      PAW%Etotal=PAW%tkin+PAW%tion+PAW%tvale+PAW%txc+PAW%Ea
      Write(STD_OUT,*) '*******Total energy*********', PAW%Etotal
      Write(STD_OUT,*) 'PAW%tkin                    ',PAW%tkin
      Write(STD_OUT,*) 'PAW%tion                    ',PAW%tion
      Write(STD_OUT,*) 'PAW%tvale                   ',PAW%tvale
      Write(STD_OUT,*) 'PAW%txc                     ',PAW%txc
      Write(STD_OUT,*) 'PAW%Ea                      ',PAW%Ea

      do ib=1,PAW%nbase
         do jb=1,PAW%nbase
            if (PAW%l(ib)==PAW%l(jb)) then
               arg=PAW%ophi(:,ib)*PAW%ophi(:,jb)*v1(:) &
&                   - PAW%otphi(:,ib)*PAW%otphi(:,jb)*v2(:)
               arg(2:n)=arg(2:n)/Grid%r(2:n)
               call extrapolate(Grid,arg)
               PAW%dij(ib,jb)=PAW%dij(ib,jb)+integrator(Grid,arg,1,irc)
            endif
          enddo
      enddo


     ! complete estimate of PSO%eig(io)
     do k=1,nocc
        io=tmap(k); l=PSO%l(io)
         write(std_out,*) 'Eig before dij ', PSO%eig(io)
           do ib=1,PAW%nbase
              do jb=1,PAW%nbase
                 if (PAW%l(ib)==l.and.PAW%l(jb)==l) then
                    PSO%eig(io)=PSO%eig(io)+PAW%dij(ib,jb)*o(io,ib)*o(io,jb)
                 endif
              enddo
           enddo
         write(std_out,*) 'Eig after dij ', PSO%eig(io)
      enddo

    ! Solve inhomogeneous diffeq. and store result in tmpOrbit
    if(needvtau) then
        write(std_out,*) 'PAWIter_LDA --orbital solutions not quite correct for mgga'
    endif    
    err=0;
    do k=1,nocc
       io=tmap(k); l=PSO%l(io)
          en=PSO%eig(io)
          rhs=0
          do ib=1,PAW%nbase
             do jb=1,PAW%nbase
                if (PAW%l(ib)==l.and.PAW%l(jb)==l) then
                   rhs=rhs-PAW%otp(:,ib)*(en*PAW%oij(ib,jb)-&
&                        PAW%dij(ib,jb))*o(io,jb)
                endif
             enddo
          enddo
          tmpOrbit%wfn(:,io)=0
            !  note rhs is -(what we expect)
          CALL inhomo_bound_numerov(Grid,l,n,en,rv,rhs,tmpOrbit%wfn(:,io))
          arg=(PSO%wfn(:,io)-tmpOrbit%wfn(:,io))**2
          err=err+PSO%occ(io)*Integrator(Grid,arg)
          !do i=1,Grid%n
          !   write(800+k,'(1p,20e15.7)') Grid%r(i),PSO%wfn(i,io),&
          !&     tmpOrbit%wfn(i,io),rhs(i),rv(i)
          !enddo
    enddo

    write(std_out,*) 'PAWIter ', fcount,err
       ! update wfn if tolerance not satisfied
       IF (err>err0) THEN
          mix1=MAX(mix,mix/err)
          val=(1.d0-mix1)
          WRITE(STD_OUT,*) 'mixing wfns ', val
          DO k=1,nocc
             io=tmap(k)
                PSO%wfn(:,io)=val*PSO%wfn(:,io)+mix1*tmpOrbit%wfn(:,io)
          ENDDO
       ELSE
          success=.TRUE.
       ENDIF
     fcount=fcount+1
     DEALLOCATE(arg,rhs,rv,aden,v1,v2,tmap,o)
     Deallocate(t,vt)
     Call DestroyOrbit(tmpOrbit)
  END SUBROUTINE PAWIter_LDA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  exploreparms reads ndata sets of PAW parameters and calculates
!    the corresponding logderivs   It is intended that this routine
!    would be used to choose the parameters with the best logderivs.
!    It should be called after the program has first run a "normal"
!    calculation.   exploreparms recalculates the wfni (i=1,2,..) and
!    logderivl (l=0,1,2,..) values but does not generate the full
!    PAW dataset.   Upon analyzing the logderiv results, the atompaw
!    program should be run in the "normal" mode to generate the best
!    dataset.  exploreparams should only be called once.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE exploreparms(Grid,Pot,FC,Orbit,PAW)
     Type(GridInfo), INTENT(IN) :: Grid
     Type(PotentialInfo), INTENT(IN) :: Pot
     Type(FCInfo), INTENT(IN) :: FC
     Type (OrbitInfo), INTENT(IN) :: Orbit
     Type(PseudoInfo), INTENT(INOUT) :: PAW

     INTEGER :: ndata,i,j,k,l,nrcs,ircs,ictrc,iprev
     INTEGER :: ifen=9
     CHARACTER(4) ::fdata
     CHARACTER(132) :: inputline,keyword
     REAL(8), allocatable ::  logderiverror(:,:)
     REAL(8) :: thisrc
     REAL(8) :: EBEGIN=-10, EEND=10       ! default logderiv range (Ry)
     LOGICAL :: success
     TYPE(input_dataset_t) :: backup_dataset

     Type rcresults
         INTEGER :: beginindex, endindex
         REAL(8) :: rc
     End type rcresults   

     Type (rcresults), POINTER :: rcsummary(:)

     success=.true.
     write(std_out,*) 'Input the number of PAW parameter sets in this run'
     READ(STD_IN,'(a)') inputline
     CALL eliminate_comment(inputline)
     read(inputline,*) ndata,nrcs
     if (ndata>9999) then
          write(std_out,*) 'Error : ndata must be <= 9999', ndata
          stop
     endif        
     if (nrcs<1.or.nrcs>9999) then
          write(std_out,*) 'Error : nrcs must be <= 9999 and >=1', nrcs
          stop
     endif        

     call UpperCase(inputline)
     i=0; i=INDEX(inputline,'LOGDERIVERANGE')
       if (i>0) then
          read(unit=inputline(i+14:80),fmt=*,err=111,end=111,iostat=i) &
&             EBEGIN , EEND
 111      continue            
        endif
               
     allocate(rcsummary(nrcs))
     allocate(logderiverror(6,ndata))         ! 6 represents the max l+1
     logderiverror=9.d20

     write(std_out,*) 'Begin explore runs'
     OPEN(20,file='EXPLORERESULTS',form='formatted')
     OPEN(21,file='EXPLORESUMMARY',form='formatted')
     write(20,'("#dataset","    rc       ",6(2x,i3,10x))') (l,l=0,PAW%lmax+1)
     write(21,'(" Logderiv errors based on energy range", 2f12.2)') &
&        EBEGIN, EEND

!    Save input dataset
     CALL input_dataset_copy(input_dataset,backup_dataset)

     ircs=0; thisrc=-1
     do i=1,ndata
        write(std_out,*) '===================== #',i,'=================='
        call mkname(i,fdata)
        OPEN(ifen, file='EXPLOREout.'//TRIM(fdata), form='formatted')


        !Read new basis parameters
        CALL input_dataset_read(echofile='EXPLOREIN.'//TRIM(fdata),&
&            read_global_data=.false.,read_elec_data=.false.,&
&            read_coreval_data=.false.,read_basis_data=.true.)

        Call InitPAW(PAW,Grid,FCOrbit,input_dataset%lmax)
        CALL SetPAWOptions1(ifen,Grid)
        CALL setbasis(Grid,FCPot,FCOrbit)
        Call setcoretail(Grid,FC%coreden)
        Call setttau(Grid,FC%coreden,FC%coretau)
        If (TRIM(FCorbit%exctype)=='HF'.or.TRIM(FCorbit%exctype)=='EXXKLI') PAW%tcore=0
        If (TRIM(FCorbit%exctype)=='EXXKLI') Call fixtcorewfn(Grid,PAW)
        Call SetPAWOptions2(ifen,Grid,FCOrbit,FCPot,success)
        Call Report_PseudobasisRP(Grid,PAW,ifen,fdata)
        if (success) then
           Call Set_PAW_MatrixElements(Grid,PAW,ifen)
           CALL EXPLORElogderiv(Grid,FCPot,PAW,fdata,EBEGIN,EEND,&
&           logderiverror(:,i))
        endif
        write(20,'(i5,2x,f12.5,1p,6e15.7)')&
&             i,PAW%rc,(logderiverror(l+1,i),l=0,PAW%lmax+1)
        close(ifen)
        Call DestroyPAW(PAW)
        if (abs(PAW%rc-thisrc)>1.d-10 ) then
           ircs=ircs+1      
           rcsummary(ircs)%beginindex=i
           rcsummary(ircs)%endindex=i
           thisrc=PAW%rc
           rcsummary(ircs)%rc=thisrc
        else
           rcsummary(ircs)%rc=thisrc
           rcsummary(ircs)%endindex=i
        endif        
           write(std_out,*) 'test ', i,ircs, rcsummary(ircs)%beginindex,&
&             rcsummary(ircs)%endindex          

     enddo

!    Restore input dataset
     CALL input_dataset_copy(backup_dataset,input_dataset)
    
     Write(STD_OUT,*) 'Results for minimum logderiverror'
     Write(21,*) 'Results for minimum logderiverror'
     Do ircs=1,nrcs
        write(std_out,'( "=== Rc = ", f20.6, " ====")') rcsummary(ircs)%rc
        write(21,'( "=== Rc = ", f20.6, " ====")') rcsummary(ircs)%rc
        j=rcsummary(ircs)%beginindex
        k=rcsummary(ircs)%endindex
        Do l=0,PAW%lmax+1
           i=MINLOC(logderiverror(l+1,j:k),1)+j-1
           Write(STD_OUT,'(" l =", i5,2x, i6, 1pe15.7)') l,i,logderiverror(l+1,i)
           Write(21,'(" l =", i5,2x, i6, 1pe15.7)') l,i,logderiverror(l+1,i)
        enddo  
     enddo   

     CLOSE(20); CLOSE(21)
     Deallocate(logderiverror,rcsummary)

 END SUBROUTINE exploreparms


    !************************************************************************
    !  program to calculate logerivatives of paw wavefunctions
    !   and to compare them with all electron wavefunctions
    !  optionally, wavefunctions are written to a file
    !  Assumes prior call to SUBROUTINE Set_PAW_MatrixElements(Grid,PAW)
    !    Version for exploreparms
    !************************************************************************
    SUBROUTINE EXPLORElogderiv(Grid,Pot,PAW,label,EBEGIN,EEND,lderror)
      TYPE(GridInfo), INTENT(IN) :: Grid
      TYPE(PotentialInfo), INTENT(IN) :: Pot
      TYPE(PseudoInfo), INTENT(INOUT) :: PAW
      CHARACTER(4), INTENT(IN) :: label
      REAL(8), INTENT(IN) :: EBEGIN,EEND
      REAL(8), INTENT(INOUT) :: lderror(:)

      TYPE(PotentialInfo) :: PS
      INTEGER :: n,l,ie,nbase,ib,ic,nr,i,nodes,mbase,irc,lng,ne
      REAL(8), PARAMETER :: e0=-5.d0,de=0.05d0
      REAL(8) :: h,x,dwdr,dcwdr,scale,energy
      REAL(8), ALLOCATABLE :: psi(:),tpsi(:),ttpsi(:)
      CHARACTER(4)  :: flnm
      LOGICAL :: OK

      n=Grid%n; h=Grid%h; nbase=PAW%nbase; irc=PAW%irc;  nr=irc+10
      ne=(EEND-EBEGIN)/de + 1

      ALLOCATE(psi(nr),tpsi(nr),ttpsi(nr),PS%rv(n),stat=ie)
      IF (ie/=0) THEN
         WRITE(STD_OUT,*) 'Error in logderiv allocation',n,ie
         STOP
      ENDIF

      ! load  PS
      PS%rv=PAW%rveff ; PS%nz=0.d0
      call zeropot(Grid,PS%rv,PS%v0,PS%v0p)
      !
      !   calculate logderivatives at irc
      !
      WRITE(STD_OUT,*) 'calculating log derivatives at irc',Grid%r(irc)
      !
      mbase=nbase; irc=PAW%irc
      write(std_out,*) 'Nodes counted to radius ', Grid%r(irc)
      DO l=0,PAW%lmax+1
         OK=.true.
         If (l<=PAW%lmax) then
         Basislist: do ib=1,nbase
             if (PAW%l(ib)==l) then
                nodes=countnodes(2,irc,PAW%ophi(:,ib))
                if (nodes/=PAW%nodes(ib).and.PAW%eig(ib)>0) then
                   lderror(l+1)=9.d20 
                   write(std_out,*) 'Warning node problems for case ', l,ib,&
&                     nodes,PAW%nodes(ib) 
                   OK=.false.
                   exit Basislist
                endif 
            endif 
        enddo Basislist
        endif
        if (OK) then
            CALL mkname(l,flnm)
            OPEN(56,file=TRIM(label)//'.logderiv.'//TRIM(flnm),form='formatted')
            lderror(l+1)=0

          DO ie=1,ne
             energy=EBEGIN+de*(ie-1)
             psi=0;tpsi=0;ttpsi=0
             if (scalarrelativistic) then
                CALL unboundsr(Grid,Pot,nr,l,energy,psi,nodes)
             elseif (TRIM(PAW%exctype)=='HF') then
                CALL HFunocc(Grid,PAW%OCCWFN,l,energy,Pot%rv,Pot%v0,Pot%v0p,&
                      psi,lng)    
             else
               CALL unboundsch(Grid,Pot%rv,Pot%v0,Pot%v0p,nr,l,energy,psi,nodes)
             endif
               CALL unboundsep(Grid,PS,PAW,nr,l,energy,tpsi,nodes)
               CALL PStoAE(Grid,PAW,nr,l,tpsi,ttpsi)
            !
               dwdr=Gfirstderiv(Grid,irc,psi)/psi(irc)
               dcwdr=Gfirstderiv(Grid,irc,ttpsi)/ttpsi(irc)

               WRITE(56,'(1p,5e12.4)') energy,dwdr,dcwdr ;
               lderror(l+1)=lderror(l+1)+abs(atan(dwdr)-atan(dcwdr))
          ENDDO !ie
         CLOSE(56)
       ENDIF   !OK

      ENDDO !l

      DEALLOCATE(psi,tpsi,ttpsi,PS%rv)
    END SUBROUTINE EXPLORElogderiv

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !   Repeated version of output for used with EXPLORElogderiv
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE Report_PseudobasisRP(Grid,PAW,ifen,fdata)
     Type(GridInfo), INTENT(IN)  :: Grid
     Type(PseudoInfo), INTENT(IN) :: PAW
     INTEGER, INTENT(IN) :: ifen
     CHARACTER(4) :: fdata

     INTEGER, parameter :: ifout=15
     INTEGER :: i,j,io,nbase,irc,icount,n
     INTEGER, ALLOCATABLE :: mapp(:)
     CHARACTER (len=4) :: flnm

     nbase=PAW%nbase;irc=PAW%irc;n=Grid%n
     WRITE(ifen,'(/"Number of basis functions ",i5)') nbase
     WRITE(ifen,*)'No.   n    l      Energy         Cp coeff         Occ'

     DO io=1,nbase
        WRITE(ifen,'(3i5,1p,3e15.7)') io,PAW%np(io),PAW%l(io),PAW%eig(io),&
&       PAW%ck(io),PAW%occ(io)
        CALL mkname(io,flnm)
        OPEN(ifout,file=TRIM(fdata)//'.wfn'//TRIM(flnm),form='formatted')
        WRITE(ifout,*) '# l=',PAW%l(io),'basis function with energy  ',&
&            PAW%eig(io)
          DO i=1,irc+50
             WRITE(ifout,'(1p,5e12.4)') Grid%r(i),PAW%ophi(i,io),&
&                   PAW%otphi(i,io),PAW%otp(i,io)
          ENDDO
       CLOSE(ifout)
    ENDDO

    ! also write "raw" wavefunctions
     DO io=1,nbase
        CALL mkname(io,flnm)
        OPEN(ifout,file=TRIM(fdata)//'.wfn00'//TRIM(flnm),form='formatted')
        WRITE(ifout,*) '# l=',PAW%l(io),'basis function with energy  ',&
&            PAW%eig(io)
          DO i=1,irc+50
             WRITE(ifout,'(1p,5e12.4)') Grid%r(i),PAW%phi(i,io),&
&                   PAW%tphi(i,io),PAW%tp(i,io)
          ENDDO
       CLOSE(ifout)
    ENDDO

    allocate(mapp(PAW%OCCWFN%norbit))
    mapp=0
    icount=0
    do io=1,PAW%OCCWFN%norbit
       if(PAW%OCCWFN%iscore(io)) then
       else
         icount=icount+1
         mapp(icount)=io
       endif
    enddo
    OPEN(ifout,file='OCCWFN',form='formatted')
    WRITE(ifout,'("#            ",50i30)') (mapp(j),j=1,icount)
    do i=1,n
       write(ifout,'(1p,51e15.7)') Grid%r(i),(PAW%OCCWFN%wfn(i,mapp(j)),&
&               PAW%TOCCWFN%wfn(i,mapp(j)),j=1,icount)
    enddo
    close(ifout)

    deallocate(mapp)

  END SUBROUTINE Report_PseudobasisRP


!************************************************************************
! 7/2018 Program written by Casey Brock from Vanderbilt U.
!  Unwraps phase so it is continuous (no jumps of pi)
!  The algorithm was designed by Alan Tackett.
!
!  x, y: psi and psi' used to calculate atan
!  atan_curr, atan_prev: value of atan, and values from prev. iteration
!     These should always be unshifted.
!  s2, s1: current (s2) and previous (s1) values of atan slope
!  cumshift: an integer multiple of pi, the cumulative shift applied
!     to atan
!  cnt: the number of previous points that were shifted
!************************************************************************
SUBROUTINE phase_unwrap(ie, y, x, atan_curr, atan_prev, s2, s1, cumshift, cnt)
   INTEGER, INTENT(IN) :: ie
   REAL(8), INTENT(IN) :: x, y, s1, s2
   REAL(8), INTENT(IN) :: atan_curr, atan_prev
   REAL(8), INTENT(INOUT) :: cumshift
   INTEGER, INTENT(INOUT) :: cnt

   cnt = cnt + 1

! if (x<0) AND (slope changes sign)
   IF ((x<0) .AND. (atan_curr*atan_prev < 0)) THEN
      cnt = 0
      cumshift = cumshift - sign(pi, y)
! if (slope changes sign) AND (last two points weren't shifted) 
!    AND (absolute change in slope > 0.01)
   ELSE IF ((s1*s2 < 0) .AND. (cnt>2) .AND. (ABS(s2-s1)>0.01)) THEN
      cnt = 0
      cumshift = cumshift + sign(pi, s1)
   ENDIF
END SUBROUTINE phase_unwrap

END MODULE pseudo
