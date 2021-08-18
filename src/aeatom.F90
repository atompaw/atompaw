!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This module contains the following active subroutines:
!     SCFatom_Init,SCFatom,Orbit_Init,NC_Init,SC_Init,FC_Init,Potential_Init,
!       Set_Valence
!  The following subroutines need revision to be useful:
!       dump_aeatom,load_aeatom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

MODULE AEatom

  USE io_tools
  USE atomdata
  USE excor
  USE exx_mod
  USE general_mod
  USE globalmath
  USE gridmod
  USE hf_mod
  USE ldagga_mod
  USE input_dataset_mod
  
  IMPLICIT NONE

  REAL(8), PARAMETER, PRIVATE :: linh=0.0025d0,logh=0.020d0,lor00=1.d-5
  REAL(8), PARAMETER, PRIVATE :: small0=1.d-13,rimix=0.5d0,worst=1.d-8
  INTEGER, PARAMETER, PRIVATE :: niter=1000,mxloop=500
  REAL(8), PARAMETER, PRIVATE :: conv1=4.d13,conv2=3.d13,conv3=2.d13,conv4=1.d13
  REAL(8), PRIVATE :: electrons

  TYPE (PotentialInfo) ,TARGET :: AEPot
  TYPE (OrbitInfo) ,TARGET :: AEOrbit
  TYPE (SCFInfo) ,TARGET :: AESCF

  TYPE (PotentialInfo) ,TARGET :: FCPot
  TYPE (OrbitInfo) ,TARGET :: FCOrbit
  TYPE (SCFInfo) ,TARGET :: FCSCF
  TYPE(FCInfo),TARGET :: FC

  TYPE (Gridinfo) ,TARGET :: Grid

  TYPE (PotentialInfo) ,PRIVATE, POINTER :: PotPtr
  TYPE (OrbitInfo) ,PRIVATE, POINTER :: OrbitPtr
  TYPE (SCFInfo) ,PRIVATE, POINTER :: SCFPtr

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SCFatom_Init                      !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SCFatom_Init()

    !  program to calculate the self-consistent density functional
    !    atom ground state for atom with atomic number nz
    !    for self-consistent potential rv

    REAL(8) :: xocc,qf,small,zeff,hval
    REAL(8) :: qcal,rescale,nzeff,h,r0
    INTEGER :: icount,i,j,k,it,start,np,ierr
    INTEGER :: is,ip,id,jf,ig,io,l,nfix,ir,kappa
    INTEGER :: nps,npp,npd,npf,npg
    INTEGER, ALLOCATABLE :: nl(:,:)
    
!   Various initializations
    AEPot%v0=0.d0;AEPot%v0p=0.d0;AEPot%Nv0=0;AEPot%Nv0p=0 
    frozencorecalculation=.FALSE.;setupfrozencore=.false.
    gaussianshapefunction=.FALSE.;besselshapefunction=.FALSE.

!   Atomic symbol and atomic number (from input dataset)
    AEPot%sym=input_dataset%atomic_symbol
    AEPot%nz=input_dataset%atomic_charge
    AEPot%zz=AEPot%nz

!   Relativistic, point-nucleus, HF data (from input dataset)
    scalarrelativistic=input_dataset%scalarrelativistic
    diracrelativistic=input_dataset%diracrelativistic
    HFpostprocess=input_dataset%HFpostprocess
    finitenucleus=input_dataset%finitenucleus
    AEPot%finitenucleusmodel=input_dataset%finitenucleusmodel

!   Grid data (from input dataset)
    IF (TRIM(input_dataset%gridkey)=='LINEAR') THEN
      hval=input_dataset%gridmatch/(input_dataset%gridpoints-1)
      CALL InitGrid(Grid,hval,input_dataset%gridrange)
    ELSEIF (TRIM(input_dataset%gridkey)=='LOGGRID') THEN
      hval=logh
      CALL findh(AEPot%zz,input_dataset%gridmatch,input_dataset%gridpoints,hval,r0)
      CALL InitGrid(Grid,hval,input_dataset%gridrange,r0=r0)
    ELSEIF (TRIM(input_dataset%gridkey)=='LOGGRID4') THEN
      hval=logh
      CALL findh_given_r0(AEPot%zz,input_dataset%gridmatch,lor00,&
&                         input_dataset%gridpoints,hval)
      CALL InitGrid(Grid,hval,input_dataset%gridrange,r0=lor00/AEPot%zz)
    ENDIF
     
    CALL InitPot(AEPot,Grid%n)
    CALL Get_Nuclearpotential(Grid,AEPot)

!   Logderiv data
    minlogderiv=input_dataset%minlogderiv
    maxlogderiv=input_dataset%maxlogderiv
    nlogderiv=input_dataset%nlogderiv

!   Bound state solver
    BDsolve=input_dataset%BDsolve

!   Exchange-correlation/HF keyword (from input dataset)
    exctype=input_dataset%exctype
    localizedcoreexchange=input_dataset%localizedcoreexchange
    ColleSalvetti=.FALSE.
    IF (TRIM(input_dataset%exctype)=='EXX'.or.&
&       TRIM(input_dataset%exctype)=='EXXKLI'.or.&
&       TRIM(input_dataset%exctype)=='EXXOCC') THEN
      CALL EXX_Input_Settings(input_dataset%fixed_zero,input_dataset%fixed_zero_index)
    ELSE IF (TRIM(input_dataset%exctype)=='EXXCS') THEN
      CALL EXX_Input_Settings(input_dataset%fixed_zero,input_dataset%fixed_zero_index)
      ColleSalvetti=.TRUE.
    ELSE IF (TRIM(input_dataset%exctype)=='HF'.or.&
&            TRIM(input_dataset%exctype)=='HFV') THEN
      CALL initexch
    ELSE
      CALL initexch
    ENDIF

!   Electronic configuration of atom
    WRITE(STD_OUT,'(a,f6.2)') ' Calculation for atomic number = ',AEPot%zz
    call flush_unit(std_out)

    nps=input_dataset%np(1)
    npp=input_dataset%np(2)
    npd=input_dataset%np(3)
    npf=input_dataset%np(4)
    npg=input_dataset%np(5)
    WRITE(STD_OUT,'(5i4)') nps,npp,npd,npf,npg

    i=MAX(nps,npp,npd,npf,npg)
    j=nps
    IF(npp>0) j=j+npp-1
    IF(npd>0) j=j+npd-2
    IF(npf>0) j=j+npf-3
    IF(npg>0) j=j+npg-4

    IF (diracrelativistic) j=j+j   !  need more orbitals
    CALL InitOrbit(AEOrbit,j,Grid%n,exctype)
    AEOrbit%nps=nps;AEOrbit%npp=npp;AEOrbit%npd=npd
    AEOrbit%npf=npf;AEOrbit%npg=npg
    ALLOCATE(nl(i,j));nl=0

    icount=0
    IF (.not.diracrelativistic) THEN
      IF (nps.GT.0) THEN
        DO is=1,nps
           icount=icount+1
           nl(is,1)=icount
           AEOrbit%occ(icount)=2.d0
           AEOrbit%np(icount)=is
           AEOrbit%l(icount)=0
        ENDDO
      ENDIF
      IF (npp.GT.1) THEN
        DO ip=2,npp
           icount=icount+1
           nl(ip,2)=icount
           AEOrbit%occ(icount)=6.d0
           AEOrbit%np(icount)=ip
           AEOrbit%l(icount)=1
        ENDDO
      ENDIF
      IF (npd.GT.2) THEN
        DO id=3,npd
           icount=icount+1
           nl(id,3)=icount
           AEOrbit%occ(icount)=10.d0
           AEOrbit%np(icount)=id
           AEOrbit%l(icount)=2
        ENDDO
      ENDIF
      IF (npf.GT.3) THEN
        DO jf=4,npf
           icount=icount+1
           nl(jf,4)=icount
           AEOrbit%occ(icount)=14.d0
           AEOrbit%np(icount)=jf
           AEOrbit%l(icount)=3
        ENDDO
      ENDIF
      IF(npg.GT.4) THEN
        DO ig=5,npg
           icount=icount+1
           nl(ig,5)=icount
           AEOrbit%occ(icount)=18.d0
           AEOrbit%np(icount)=ig
           AEOrbit%l(icount)=4
        ENDDO
      ENDIF
      AEOrbit%norbit=icount
      AEOrbit%nps=nps
      AEOrbit%npp=npp
      AEOrbit%npd=npd
      AEOrbit%npf=npf
      AEOrbit%npg=npg
      IF (AEOrbit%norbit/=input_dataset%norbit) STOP 'Inconsistent number of orbitals!'
      WRITE(STD_OUT,*) AEOrbit%norbit, ' orbitals will be calculated'
 
      WRITE(STD_OUT,*)' Below are listed the default occupations '
      WRITE(STD_OUT,"(' n  l     occupancy')")
      DO io=1,AEOrbit%norbit
        WRITE(STD_OUT,'(i2,1x,i2,4x,1p,1e15.7)') &
  &          AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
      ENDDO

      !Corrected occupations (from input dataset)
      DO io=1,input_dataset%norbit_mod
        l=input_dataset%orbit_mod_l(io)
        ip=input_dataset%orbit_mod_n(io)
        xocc=input_dataset%orbit_mod_occ(io)
        nfix=nl(ip,l+1)
        IF (nfix<=0.OR.nfix>AEOrbit%norbit) THEN
          WRITE(STD_OUT,*) 'error in occupations -- ip,l,xocc: ',ip,l,xocc,nfix,AEOrbit%norbit
          STOP
        ENDIF
        AEOrbit%occ(nfix)=xocc
      END DO

      WRITE(STD_OUT,*) ' Corrected occupations are: '
      WRITE(STD_OUT,"(' n  l     occupancy')")
      electrons=0.d0
      DO io=1,AEOrbit%norbit
         WRITE(STD_OUT,'(i2,1x,i2,4x,1p,1e15.7)')  &
  &           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
         electrons=electrons+AEOrbit%occ(io)
      ENDDO
    ENDIF ! scalarrelativistic

    IF (diracrelativistic) THEN
      icount=0        
      deallocate(nl)
      i=MAX(nps,npp,npd,npf,npg)
      allocate(nl(i,-5:5))
      nl=0
      IF (nps.GT.0) THEN
         DO is=1,nps
            icount=icount+1
            nl(is,-1)=icount
            AEOrbit%occ(icount)=2.d0
            AEOrbit%np(icount)=is
            AEOrbit%l(icount)=0
            AEOrbit%kappa(icount)=-1
         ENDDO
      ENDIF
      IF (npp.GT.1) THEN
         DO ip=2,npp
            icount=icount+1
            nl(ip,1)=icount
            AEOrbit%occ(icount)=2.d0
            AEOrbit%np(icount)=ip
            AEOrbit%l(icount)=1
            AEOrbit%kappa(icount)=1
         ENDDO
         DO ip=2,npp
            icount=icount+1
            nl(ip,-2)=icount
            AEOrbit%occ(icount)=4.d0
            AEOrbit%np(icount)=ip
            AEOrbit%l(icount)=1
            AEOrbit%kappa(icount)=-2
         ENDDO
      ENDIF
      IF (npd.GT.2) THEN
         DO id=3,npd
            icount=icount+1
            nl(id,2)=icount
            AEOrbit%occ(icount)=4.d0
            AEOrbit%np(icount)=id
            AEOrbit%l(icount)=2
            AEOrbit%kappa(icount)=2
         ENDDO
         DO id=3,npd
            icount=icount+1
            nl(id,-3)=icount
            AEOrbit%occ(icount)=6.d0
            AEOrbit%np(icount)=id
            AEOrbit%l(icount)=2
            AEOrbit%kappa(icount)=-3
         ENDDO
      ENDIF
      IF (npf.GT.3) THEN
         DO jf=4,npf
            icount=icount+1
            nl(jf,3)=icount
            AEOrbit%occ(icount)=6.d0
            AEOrbit%np(icount)=jf
            AEOrbit%l(icount)=3
            AEOrbit%kappa(icount)=3
         ENDDO
         DO jf=4,npf
            icount=icount+1
            nl(jf,-4)=icount
            AEOrbit%occ(icount)=8.d0
            AEOrbit%np(icount)=jf
            AEOrbit%l(icount)=3
            AEOrbit%kappa(icount)=-4
         ENDDO
      ENDIF
      IF(npg.GT.4) THEN
         DO ig=5,npg
            icount=icount+1
            nl(ig,4)=icount
            AEOrbit%occ(icount)=8.d0
            AEOrbit%np(icount)=ig
            AEOrbit%l(icount)=4
            AEOrbit%kappa(icount)=4
         ENDDO
         DO ig=5,npg
            icount=icount+1
            nl(ig,-5)=icount
            AEOrbit%occ(icount)=10.d0
            AEOrbit%np(icount)=ig
            AEOrbit%l(icount)=4
            AEOrbit%kappa(icount)=-5
         ENDDO
      ENDIF
      AEOrbit%norbit=icount
      AEOrbit%nps=nps
      AEOrbit%npp=npp
      AEOrbit%npd=npd
      AEOrbit%npf=npf
      AEOrbit%npg=npg
      IF (AEOrbit%norbit/=input_dataset%norbit) STOP 'Inconsistent number of orbitals!'
      WRITE(STD_OUT,*) AEOrbit%norbit, ' orbitals will be calculated'

      WRITE(STD_OUT,*)' Below are listed the default occupations '
      WRITE(STD_OUT,"(' n  l kappa     occupancy')")
      DO io=1,AEOrbit%norbit
         WRITE(STD_OUT,'(i2,1x,i2,3x,i2,4x,1p,1e15.7)') &
  &           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%kappa(io),AEOrbit%occ(io)
      ENDDO

      !Corrected occupations (from input dataset)
      DO io=1,input_dataset%norbit_mod
        l=input_dataset%orbit_mod_l(io)
        ip=input_dataset%orbit_mod_n(io)
        kappa=input_dataset%orbit_mod_k(io)
        xocc=input_dataset%orbit_mod_occ(io)
        nfix=nl(ip,kappa)
        IF (nfix<=0.OR.nfix>AEOrbit%norbit) THEN
          WRITE(STD_OUT,*) 'error in occupations -- ip,l,kappa,xocc: ',ip,l,kappa,xocc,nfix,AEOrbit%norbit
          STOP
        ENDIF
         AEOrbit%occ(nfix)=xocc
      END DO

      WRITE(STD_OUT,*) ' Corrected occupations are: '
      WRITE(STD_OUT,"(' n  l  kappa   occupancy')")
      electrons=0.d0
      DO io=1,AEOrbit%norbit
         WRITE(STD_OUT,'(i2,1x,i2,3x,i2,4x,1p,1e15.7)')  &
  &           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%kappa(io),AEOrbit%occ(io)
         electrons=electrons+AEOrbit%occ(io)
      ENDDO
    ENDIF ! diracrelativistic
 
    AEPot%q=electrons
    qf=AEPot%nz-electrons
    WRITE(STD_OUT,*)
    WRITE(STD_OUT,*) 'nuclear charge    = ', AEPot%nz
    WRITE(STD_OUT,*) 'electronic charge = ', electrons
    WRITE(STD_OUT,*) 'net charge        = ', qf

    CALL InitSCF(AESCF)
    IF (scalarrelativistic) CALL Allocate_Scalar_Relativistic(Grid)
    IF (diracrelativistic)  CALL Allocate_Dirac_Relativistic(Grid)

    write(std_out,*) 'Finish SCFatom_Init' ; call flush_unit(std_out)

    DEALLOCATE(nl)

  END SUBROUTINE SCFatom_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SCFatom                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SCFatom(scftype,lotsofoutput,skip_reading)
    IMPLICIT NONE
    CHARACTER(2),INTENT(IN)::scftype
    LOGICAL,INTENT(IN) :: lotsofoutput
    LOGICAL,INTENT(IN),OPTIONAL :: skip_reading
    TYPE(SCFInfo) :: SCFPP

    !  program to perform self-consistency loop

    IF(scftype=='AE') THEN
       !------------AE calculation--------------
       frozencorecalculation=.FALSE.
       setupfrozencore=.FALSE.
       PotPtr=>AEPot
       OrbitPtr=>AEOrbit
       SCFPtr=>AESCF
       CALL Orbit_Init(OrbitPtr,PotPtr)
       CALL Potential_Init(OrbitPtr,PotPtr)

    ELSEIF(scftype=='NC') THEN
       !------------New Configuration Calculation---------
       frozencorecalculation=.FALSE.
       setupfrozencore=.FALSE.
       PotPtr=>AEPot
       OrbitPtr=>AEOrbit
       SCFPtr=>AESCF
       IF (PRESENT(skip_reading)) THEN
         CALL NC_Init(OrbitPtr,PotPtr,skip_reading=skip_reading)
       ELSE
         CALL NC_Init(OrbitPtr,PotPtr)
       ENDIF
       CALL Potential_Init(OrbitPtr,PotPtr)

    ELSEIF(scftype=='SC') THEN
       !------------Set Core and Valence in current config.---------
       frozencorecalculation=.TRUE.
       setupfrozencore=.TRUE.
       CALL CopyPot(AEPot,FCPot)
       CALL CopyOrbit(AEOrbit,FCOrbit)
       CALL CopySCF(AESCF,FCSCF)
       PotPtr=>FCPot
       OrbitPtr=>FCOrbit
       SCFPtr=>FCSCF
       CALL SC_Init(OrbitPtr,PotPtr,SCFPtr)

    ELSEIF(scftype=='FC') THEN
       !------------Frozen Core Calculation---------
       frozencorecalculation=.TRUE.
       setupfrozencore=.FALSE.
       PotPtr=>FCPot
       OrbitPtr=>FCOrbit
       SCFPtr=>FCSCF
       CALL FC_Init(OrbitPtr,PotPtr)
       CALL Potential_Init(OrbitPtr,PotPtr)

    ENDIF

    IF (TRIM(OrbitPtr%exctype)=='EXX'.OR. &
&       TRIM(OrbitPtr%exctype)=='EXXKLI'.OR. &
        TRIM(OrbitPtr%exctype)=='EXXCS') THEN
       CALL EXX_SCF(scftype,lotsofoutput,Grid,OrbitPtr,PotPtr,FC,SCFPtr)

    ELSE IF (TRIM(OrbitPtr%exctype)=='EXXOCC') THEN
       CALL EXXOCC_SCF(scftype,lotsofoutput,Grid,OrbitPtr,PotPtr,FC,SCFPtr)

    ELSE IF (TRIM(OrbitPtr%exctype)=='HF'.or.TRIM(OrbitPtr%exctype)=='HFV') THEN
       write(std_out,*) 'Just before HF_SCF'; call flush_unit(std_out)
       CALL HF_SCF(scftype,lotsofoutput,Grid,OrbitPtr,PotPtr,FC,SCFPtr)
       write(std_out,*) 'Just after HF_SCF'; call flush_unit(std_out)

    ELSE
      !LDA or GGA
      CALL LDAGGA_SCF(scftype,lotsofoutput,Grid,OrbitPtr,PotPtr,FC,SCFPtr)
    ENDIF

    If (HFpostprocess) then
       write(std_out,*) "******PostProcessing HF*********"
       CALL InitSCF(SCFPP)
       call hf_energy_only(Grid,OrbitPtr,PotPtr,SCFPP)
    endif
  END SUBROUTINE SCFatom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Orbit_Init
  !!       From nuclear charge -- generate hydrogenic-like initial wfns
  !!          and densities --
  !!          fill AEOrbit%wfn, AEOrbit%eig, and AEOrbit%den and AEOrbit%q
  !!          also AEOrbit%otau and AEOrbit%tau
  !!          Note that both den and tau need to be divided by 4 \pi r^2
  !!          Note that tau is only correct for non-relativistic case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Orbit_Init(Orbit,Pot)
    IMPLICIT NONE
    TYPE(PotentialInfo),INTENT(INOUT) :: Pot
    TYPE(OrbitInfo),INTENT(INOUT) :: Orbit
    INTEGER  :: i,io,ir,ip,l,nfix,np,kappa
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,xocc,fac
    REAL(8),allocatable :: dpdr(:),pbr(:)
    INTEGER :: initialconfig=0

       IF (initialconfig/=0) STOP 'Error in aeatom -- Orbit_Init already called'

       !  calculate initial charge density from hydrogen-like functions
       !  also initial energies
       zeff=Pot%nz
       If  (.not.diracrelativistic) then
       DO io=1,Orbit%norbit
          np=Orbit%np(io)
          l=Orbit%l(io)
          xocc=Orbit%occ(io)
          Orbit%eig(io)=-(zeff/(np))**2
          WRITE(STD_OUT,*) io,np,l,xocc,Orbit%eig(io)
          DO ir=1,Grid%n
             Orbit%wfn(ir,io)=hwfn(zeff,np,l,Grid%r(ir))
             IF (ABS(Orbit%wfn(ir,io))<machine_zero) Orbit%wfn(ir,io)=0.d0
          ENDDO
          zeff=zeff-0.5d0*xocc
          zeff=MAX(zeff,1.d0)
       ENDDO
       endif
       If  (diracrelativistic) then
       DO io=1,Orbit%norbit
          np=Orbit%np(io)
          l=Orbit%l(io)
          kappa=Orbit%kappa(io)
          xocc=Orbit%occ(io)
          Orbit%eig(io)=-(zeff/(np))**2
          WRITE(STD_OUT,*) io,np,l,xocc,Orbit%eig(io)
          DO ir=1,Grid%n
           call dirachwfn(np,kappa,zeff,Grid%r(ir),Orbit%eig(io) &
&               ,Orbit%wfn(ir,io),Orbit%lwfn(ir,io))
            IF (ABS(Orbit%wfn(ir,io))<machine_zero) Orbit%wfn(ir,io)=0.d0
            IF (ABS(Orbit%lwfn(ir,io))<machine_zero) Orbit%lwfn(ir,io)=0.d0
          ENDDO
          zeff=zeff-0.5d0*xocc
          zeff=MAX(zeff,1.d0)
       ENDDO
       endif



       ! check charge and rescale
       allocate(dpdr(Grid%n),pbr(Grid%n))
       Orbit%den=0.d0;Orbit%tau=0.d0
       DO io=1,Orbit%norbit
          dpdr=0.d0;pbr=0.d0
          pbr(2:Grid%n)=Orbit%wfn(2:Grid%n,io)/Grid%r(2:Grid%n)
          CALL derivative(Grid,pbr,dpdr,2,Grid%n)
          CALL extrapolate(Grid,pbr)
          CALL extrapolate(Grid,dpdr)
          l=Orbit%l(io)
          fac=l*(l+1)
          Orbit%otau(:,io)=(Grid%r(:)*dpdr(:))**2+fac*(pbr(:))**2
          xocc=Orbit%occ(io)
          DO ir=1,Grid%n
             Orbit%den(ir)=Orbit%den(ir)+xocc*(Orbit%wfn(ir,io)**2)
             Orbit%tau(ir)=Orbit%tau(ir)+xocc*Orbit%otau(ir,io)
             If (diracrelativistic) Orbit%den(ir)=Orbit%den(ir) + &
&                 xocc*((Orbit%lwfn(ir,io))**2)                     
          ENDDO
       ENDDO
!   Note that kinetic energy density (tau) is in Rydberg units
!   Note that kinetic energy density is only correct for non-relativistic
!               formulation
       qcal=integrator(Grid,Orbit%den)
       qf=qcal
       !WRITE(STD_OUT,*) 'qcal electrons = ',qcal, electrons
       !rescale density
       rescale=electrons/qcal
       Orbit%den(1:Grid%n)=Orbit%den(1:Grid%n)*rescale
       Orbit%tau(1:Grid%n)=Orbit%tau(1:Grid%n)*rescale

       deallocate(dpdr,pbr)
       initialconfig=1
        write(std_out,*) 'completed Orbit_Init '; call flush_unit(std_out)

  END SUBROUTINE Orbit_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE NC_Init(Orbit,Pot,skip_reading)
    TYPE(PotentialInfo),INTENT(INOUT) :: Pot
    TYPE(OrbitInfo),INTENT(INOUT) :: Orbit
    LOGICAL,INTENT(IN),OPTIONAL :: skip_reading
    INTEGER  :: i,io,jo,ir,ip,l,nfix,j,norbit_mod,np(5)
    LOGICAL :: read_new_occ
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,xocc,fac
    REAL(8),allocatable :: dpdr(:),pbr(:)
    INTEGER, ALLOCATABLE :: orbit_mod_l(:),orbit_mod_n(:),orbit_mod_k(:)
    REAL(8), ALLOCATABLE :: orbit_mod_occ(:)

    read_new_occ=.true.;if (present(skip_reading)) read_new_occ=.not.skip_reading
    
    !----------New Configuration--AE calculation--------------
       ! readin revision and normalize the density

       IF (read_new_occ) THEN
  
         !Read modified occupations from standard input
         np(1)=Orbit%nps;np(2)=Orbit%npp;np(3)=Orbit%npd;np(4)=Orbit%npf;np(5)=Orbit%npg
         CALL input_dataset_read_occ(norbit_mod,orbit_mod_l,orbit_mod_n,orbit_mod_k,&
  &                                  orbit_mod_occ,np,diracrelativistic)

         DO jo=1,norbit_mod
           nfix=-100
           DO io=1,Orbit%norbit
             IF (orbit_mod_n(jo)==Orbit%np(io).AND.orbit_mod_l(jo)==Orbit%l(io).AND.&
  &              ((.NOT.diracrelativistic).OR.orbit_mod_k(jo)==Orbit%kappa(io))) THEN
               nfix=io
               EXIT
             ENDIF
           ENDDO
           IF (nfix.LE.0.OR.nfix.GT.Orbit%norbit) THEN
             WRITE(STD_OUT,*) 'error in occupations -- ip,l,xocc',&
              orbit_mod_n(jo),orbit_mod_l(jo),orbit_mod_occ(jo),nfix,Orbit%norbit
             STOP
           ENDIF
           Orbit%occ(nfix)=orbit_mod_occ(jo)
         ENDDO
         DEALLOCATE(orbit_mod_l)
         DEALLOCATE(orbit_mod_n)
         DEALLOCATE(orbit_mod_k)
         DEALLOCATE(orbit_mod_occ)
       ENDIF

       electrons=0.d0
       WRITE(STD_OUT,*) ' Corrected occupations are: '
       IF (.NOT.diracrelativistic) THEN
         WRITE(STD_OUT,"(' n  l     occupancy')")
         DO io=1,Orbit%norbit
            WRITE(STD_OUT,'(i2,1x,i2,4x,1p,1e15.7)') Orbit%np(io),Orbit%l(io),Orbit%occ(io)
            electrons=electrons+Orbit%occ(io)
         ENDDO
       ELSE
         WRITE(STD_OUT,"(' n  l kap     occupancy')")
         DO io=1,Orbit%norbit
            WRITE(STD_OUT,'(i2,1x,i2,2x,i2,4x,1p,1e15.7)') Orbit%np(io),Orbit%l(io),Orbit%kappa(io),Orbit%occ(io)
            electrons=electrons+Orbit%occ(io)
         ENDDO
       ENDIF
       Pot%q=electrons
       qf=Pot%nz-electrons
       WRITE(STD_OUT,*)
       WRITE(STD_OUT,*) 'nuclear charge    = ' , Pot%nz
       WRITE(STD_OUT,*) 'electronic charge = ', electrons
       WRITE(STD_OUT,*) 'net charge        = ', qf
       !
       !
       !  calculate initial charge density from stored wavefunctions
       !    also initial energies
       !
       allocate(dpdr(Grid%n),pbr(Grid%n))
       Orbit%den=0.d0;Orbit%tau=0.d0
       DO io=1,Orbit%norbit
          dpdr=0.d0;pbr=0.d0
          pbr(2:Grid%n)=Orbit%wfn(2:Grid%n,io)/Grid%r(2:Grid%n)
          CALL derivative(Grid,pbr,dpdr,2,Grid%n)
          CALL extrapolate(Grid,pbr)
          CALL extrapolate(Grid,dpdr)
          l=Orbit%l(io)
          fac=l*(l+1)
          Orbit%otau(:,io)=(Grid%r(:)*dpdr(:))**2+fac*(pbr(:))**2
          xocc=Orbit%occ(io)
          DO ir=1,Grid%n
             Orbit%den(ir)=Orbit%den(ir)+xocc*(Orbit%wfn(ir,io)**2)
             Orbit%tau(ir)=Orbit%tau(ir)+xocc*Orbit%otau(ir,io)
             If (diracrelativistic) Orbit%den(ir)=Orbit%den(ir) + &
&                 xocc*((Orbit%lwfn(ir,io))**2)                     
          ENDDO
       ENDDO
!   Note that kinetic energy density (tau) is in Rydberg units
!   Note that kinetic energy density is only correct for non-relativistic
!               formulation
       !
       !  check charge
       !
       qcal=integrator(Grid,Orbit%den)
       qf=qcal
       !WRITE(STD_OUT,*) 'qcal electrons = ',qcal, electrons
       !  rescale density
       rescale=electrons/qcal
       Orbit%den=Orbit%den*rescale
       Orbit%tau=Orbit%tau*rescale

       deallocate(dpdr,pbr)
  END SUBROUTINE NC_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!         SC_Init -- set core states in current configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SC_Init(Orbit,Pot,SCF)
    TYPE(PotentialInfo),INTENT(INOUT) :: Pot
    TYPE(OrbitInfo),INTENT(INOUT) :: Orbit
    TYPE (SCFInfo),INTENT(INOUT) :: SCF
    INTEGER :: choosevalence=0

    INTEGER :: io
    LOGICAL :: noalt=.true.
    If (choosevalence/=0) then
     write(std_out,*) 'Error in aeatom -- SC_Init already called'
       stop
    endif

    CALL Set_Valence()
    Call Core_Electron_Report(Orbit,FC,std_out)
    Call Valence_Electron_Report(Orbit,FC,std_out)

    IF (TRIM(Orbit%exctype)=='EXX') THEN
       Call Get_FCKinCoul(Grid,Pot,Orbit,FC,SCF)
       CALL Get_FCEnergy_EXX(Grid,Orbit,FC,SCF)
       CALL Set_Vxref(AEPot)   ! probably not a good idea
    ELSEIF (TRIM(Orbit%exctype)=='EXXKLI') THEN
       Call Get_FCKinCoul(Grid,Pot,Orbit,FC,SCF,noalt)
       CALL Get_FCEnergy_EXX(Grid,Orbit,FC,SCF)
    ELSEIF (TRIM(Orbit%exctype)=='HF'.or.TRIM(Orbit%exctype)=='HFV') THEN
       write(std_out,*) ' Completed Setcore '; call flush_unit(std_out)
       CALL Get_FCEnergy_HF(Grid,Pot,Orbit,FC,SCF)
       CALL Total_FCEnergy_Report(SCF,std_out)
       write(std_out,*) ' Completed Setcore '; call flush_unit(std_out)
    ELSE
       Call Get_FCKinCoul(Grid,Pot,Orbit,FC,SCF)
       CALL Get_FCEXC(SCF)
       CALL Total_FCEnergy_Report(SCF,std_out)
    ENDIF
    !CALL Total_Energy_Report(SCF,std_out)
    !   For EXX, exchange contribution is not known yet.

    choosevalence=1

  END SUBROUTINE SC_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE FC_Init(Orbit,Pot)
    TYPE(PotentialInfo),INTENT(INOUT) :: Pot
    TYPE(OrbitInfo),INTENT(INOUT) :: Orbit
    INTEGER  :: i,io,jo,ir,ip,l,nfix,j,norbit_mod,np(5)
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,xocc
    INTEGER, ALLOCATABLE :: orbit_mod_l(:),orbit_mod_n(:),orbit_mod_k(:)
    REAL(8), ALLOCATABLE :: orbit_mod_occ(:)

       WRITE(STD_OUT,*)' Below are listed the current valence occupations '
       IF (.NOT.diracrelativistic) THEN
         WRITE(STD_OUT,"(' n  l     occupancy')")
         DO io=1,Orbit%norbit
            IF (.NOT.Orbit%iscore(io)) WRITE(STD_OUT,'(i2,1x,i2,4x,1p,1e15.7)') &
&                Orbit%np(io),Orbit%l(io),Orbit%occ(io)
         ENDDO
       ELSE
         WRITE(STD_OUT,"(' n  l kap     occupancy')")
         DO io=1,Orbit%norbit
            IF (.NOT.Orbit%iscore(io)) WRITE(STD_OUT,'(i2,1x,i2,2x,i2,4x,1p,1e15.7)') &
&                Orbit%np(io),Orbit%l(io),Orbit%kappa(io),Orbit%occ(io)
         ENDDO
       ENDIF

       !Read modified occupations from standard input
       np(1)=Orbit%nps;np(2)=Orbit%npp;np(3)=Orbit%npd;np(4)=Orbit%npf;np(5)=Orbit%npg
       CALL input_dataset_read_occ(norbit_mod,orbit_mod_l,orbit_mod_n,orbit_mod_k,&
&                                  orbit_mod_occ,np,diracrelativistic)

       DO jo=1,norbit_mod
         nfix=-100
         DO io=1,Orbit%norbit
           IF (orbit_mod_n(jo)==Orbit%np(io).AND.orbit_mod_l(jo)==Orbit%l(io).AND.&
&              ((.NOT.diracrelativistic).OR.orbit_mod_k(jo)==Orbit%kappa(io))) THEN
             nfix=io
             EXIT
           ENDIF
         ENDDO
         IF (nfix.LE.0.OR.nfix.GT.Orbit%norbit) THEN
           WRITE(STD_OUT,*) 'error in occupations -- ip,l,xocc',&
            orbit_mod_n(jo),orbit_mod_l(jo),orbit_mod_occ(jo),nfix,Orbit%norbit
           STOP
         ENDIF
         Orbit%occ(nfix)=orbit_mod_occ(jo)
       ENDDO
       DEALLOCATE(orbit_mod_l)
       DEALLOCATE(orbit_mod_n)
       DEALLOCATE(orbit_mod_k)
       DEALLOCATE(orbit_mod_occ)

       FC%zvale=0.d0
       WRITE(STD_OUT,*) ' Corrected occupations are: '
       IF (.NOT.diracrelativistic) THEN
         WRITE(STD_OUT,"(' n  l     occupancy')")
         DO io=1,Orbit%norbit
           If (.not.Orbit%iscore(io))then
             WRITE(STD_OUT,'(i2,1x,i2,4x,1p,1e15.7)')  &
&              Orbit%np(io),Orbit%l(io),Orbit%occ(io)
             FC%zvale=FC%zvale+Orbit%occ(io)
           Endif
         ENDDO
       ELSE
         WRITE(STD_OUT,"(' n  l kap     occupancy')")
         DO io=1,Orbit%norbit
           If (.not.Orbit%iscore(io))then
             WRITE(STD_OUT,'(i2,1x,i2,2x,i2,4x,1p,1e15.7)')  &
&              Orbit%np(io),Orbit%l(io),Orbit%kappa(io),Orbit%occ(io)
             FC%zvale=FC%zvale+Orbit%occ(io)
           Endif
         ENDDO
       ENDIF

       electrons=FC%zvale+FC%zcore

       !
       !  calculate initial charge density from stored wavefunctions
       !    also initial energies
       !
       FC%valeden=0; FC%valetau=0
       DO io=1,Orbit%norbit
          IF (.NOT.Orbit%iscore(io)) THEN
             xocc=Orbit%occ(io)
             DO ir=1,Grid%n
                FC%valeden(ir)=FC%valeden(ir)+xocc*(Orbit%wfn(ir,io)**2)
                FC%valetau(ir)=FC%valetau(ir)+xocc*Orbit%otau(ir,io)
              If (diracrelativistic) FC%valeden(ir)=FC%valeden(ir) + &
&                 xocc*((Orbit%lwfn(ir,io))**2)

             ENDDO
          ENDIF
       ENDDO
       !
       !  check charge                                !  rescale density

       !
       qcal=integrator(Grid,FC%valeden)
       if (ABS(qcal)<small0.OR.electrons<small0) then
          FC%valeden=0
       else
          rescale=FC%zvale/qcal
          FC%valeden=FC%valeden*rescale
          FC%valetau=FC%valetau*rescale
       endif

       Orbit%den=FC%valeden+FC%coreden
       WRITE(STD_OUT,*)
       WRITE(STD_OUT,*) 'core charge    = ' , FC%zcore
       WRITE(STD_OUT,*) 'valence charge = ',  FC%zvale
       WRITE(STD_OUT,*) 'net charge     = ',  Pot%nz-FC%zcore-FC%zvale

  END SUBROUTINE FC_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Potential_Init
  !!    Generate Potptr%rv, rvh, rvx given AEOrbit%wfn, AEOrbit%den, AEOrbit%q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Potential_Init(Orbit,Pot)
    IMPLICIT NONE
    TYPE(OrbitInfo),INTENT(IN) :: Orbit
    TYPE(PotentialInfo),INTENT(INOUT) :: Pot
    INTEGER  :: i,io,ir,xocc,ip,l,nfix,j,np
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,ecoul,v0,etxc,eex
    LOGICAL :: success


    CALL poisson(Grid,Pot%q,Orbit%den,Pot%rvh,ecoul,v00=v0)
    write(std_out,*) 'In Potential_Init', Pot%q,ecoul; call flush_unit(std_out)

  END SUBROUTINE Potential_Init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Set Valence Orbitals
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Set_Valence()

    INTEGER :: io,n,ok

    call InitFC(FC,Grid%n)

    FC%zcore=0; FC%zvale=0
    FC%coreden=0; FC%valeden=0
    FC%coretau=0; FC%valetau=0

    DO io=1,FCOrbit%norbit
       WRITE(STD_OUT,'(3i5,1p,2e15.7)') io,FCOrbit%np(io),FCOrbit%l(io),&
&           FCOrbit%occ(io),FCOrbit%eig(io)

       FCOrbit%iscore(io)=input_dataset%orbit_iscore(io)

       IF (FCOrbit%iscore(io)) THEN
          FC%zcore=FC%zcore+FCOrbit%occ(io)
          FC%coreden=FC%coreden+FCOrbit%occ(io)*(FCOrbit%wfn(:,io))**2
          FC%coretau=FC%coretau+FCOrbit%occ(io)*FCOrbit%otau(:,io)
          If (diracrelativistic) FC%coreden=FC%coreden + &
&                 FCOrbit%occ(io)*((FCOrbit%lwfn(:,io))**2)

       ENDIF
       IF (.NOT.FCOrbit%iscore(io)) THEN

          FC%zvale=FC%zvale+FCOrbit%occ(io)
          FC%valeden=FC%valeden+FCOrbit%occ(io)*(FCOrbit%wfn(:,io))**2
          FC%valetau=FC%valetau+FCOrbit%occ(io)*FCOrbit%otau(:,io)
          If (diracrelativistic) FC%valeden=FC%valeden + &
&                 FCOrbit%occ(io)*((FCOrbit%lwfn(:,io))**2)
       ENDIF
    ENDDO

         write(std_out,*) 'Returning from SetValence'
         write(std_out,*) 'Core electrons', FC%zcore,integrator(Grid,FC%coreden)
         write(std_out,*) 'Vale electrons', FC%zvale,integrator(Grid,FC%valeden)

         write(std_out,*) 'Coretau:',FC%coretau(1:20)
         write(std_out,*) 'Valetau:',FC%valetau(1:20)
  END SUBROUTINE Set_Valence

! dump and load subroutines need to be updated!!!!  
  SUBROUTINE dump_aeatom(Fn,Grid,AEOrbit,AEPot,AESCF,FCOrbit,FCPot,FCSCF,FC)
     Character(*), INTENT(IN) :: Fn
     Type(GridInfo), INTENT(IN) :: Grid
     Type(OrbitInfo), INTENT(IN) :: AEOrbit,FCOrbit
     Type(PotentialInfo), INTENT(IN) :: AEPot,FCPot
     Type(SCFInfo), INTENT(IN) :: AESCF,FCSCF
     Type(FCInfo), INTENT(IN) :: FC

     INTEGER, parameter :: ifo=15
     INTEGER :: i,io,n

     Write(STD_OUT,*) 'Note that dump has been called; code not completely checked'

     open(ifo,file=Fn,form='formatted',status='replace')
     write(std_out,*) 'Creating or replacing dump file', Fn

     write(ifo,*) 'ATOMPAW'    ! keyword

    ! grid info
     write(ifo,*) Grid%TYPE,Grid%n,Grid%h
     n=Grid%n
     write(ifo,*)(Grid%r(i),i=1,n)
     if (usingloggrid(Grid)) then
          write(ifo,*)(Grid%drdu(i),Grid%pref(i),Grid%rr02(i),i=1,n)
     endif

    ! atomdata fixed constants
     write(ifo,*)  frozencorecalculation,setupfrozencore,scalarrelativistic,&
&                  finitenucleus,gaussianshapefunction,besselshapefunction,&
&                  ColleSalvetti

    ! Orbit info
     write(ifo,*) AEOrbit%exctype,AEOrbit%nps, AEOrbit%npp, AEOrbit%npd ,&
&      AEOrbit%npf, AEOrbit%npg, AEOrbit%norbit
     write(ifo,*) (AEOrbit%np(io),AEOrbit%l(io),AEOrbit%iscore(io),&
&        AEOrbit%eig(io),AEOrbit%occ(io),io=1,AEOrbit%norbit)
     write(ifo,*) ((AEOrbit%wfn(i,io),i=1,n),io=1,AEOrbit%norbit)
     write(ifo,*) (AEOrbit%den(i),i=1,n)
     write(ifo,*) (AEOrbit%tau(i),i=1,n)

     write(ifo,*) FCOrbit%exctype,FCOrbit%nps, FCOrbit%npp, FCOrbit%npd ,&
&      FCOrbit%npf, FCOrbit%npg, FCOrbit%norbit
     write(ifo,*) (FCOrbit%np(io),FCOrbit%l(io),FCOrbit%iscore(io),&
&        FCOrbit%eig(io),FCOrbit%occ(io),io=1,FCOrbit%norbit)
     write(ifo,*) ((FCOrbit%wfn(i,io),i=1,n),io=1,FCOrbit%norbit)
     write(ifo,*) (FCOrbit%den(i),i=1,n)

   ! Pot info
      write(ifo,*) AEPot%nz,AEPot%sym,AEPot%q,AEPot%v0,AEPot%v0p
      write(ifo,*) (AEPot%rv(i),AEPot%rvn(i),AEPot%rvh(i),AEPot%rvx(i),i=1,n)

      write(ifo,*) FCPot%nz,FCPot%sym,FCPot%q,FCPot%v0,FCPot%v0p
      write(ifo,*) (FCPot%rv(i),FCPot%rvn(i),FCPot%rvh(i),FCPot%rvx(i),i=1,n)

   ! SCFinfo
      write(ifo,*) AESCF%iter,AESCF%delta,AESCF%eone,AESCF%ekin,&
&       AESCF%estatic,AESCF%ecoul,AESCF%eexc,AESCF%oepcs,AESCF%etot,&
&       AESCF%valekin,AESCF%valecoul,AESCF%valeexc,AESCF%corekin,AESCF%evale

      write(ifo,*) FCSCF%iter,FCSCF%delta,FCSCF%eone,FCSCF%ekin,&
&       FCSCF%estatic,FCSCF%ecoul,FCSCF%eexc,FCSCF%oepcs,FCSCF%etot,&
&       FCSCF%valekin,FCSCF%valecoul,FCSCF%valeexc,FCSCF%corekin,FCSCF%evale

   ! FCinfo
      write(ifo,*) FC%zvale,FC%zcore
      write(ifo,*) (FC%coreden(i),FC%valeden(i),i=1,n)
      write(ifo,*) (FC%coretau(i),FC%valetau(i),i=1,n)

   If (TRIM(AEOrbit%exctype)=='EXX'.or.TRIM(AEOrbit%exctype)=='EXXKLI') &
&                  Call EXXdump(Grid,AEOrbit,ifo)
   If (TRIM(AEOrbit%exctype)=='HF'.or.TRIM(AEOrbit%exctype)=='HFV') &
&          Call HFdump(Grid,AEOrbit,ifo)

   close(ifo)

   write(std_out,*) 'Closing dump file'
  END SUBROUTINE dump_aeatom

  SUBROUTINE load_aeatom(Fn,Grid,AEOrbit,AEPot,AESCF,FCOrbit,FCPot,FCSCF,FC)
     Character(*), INTENT(IN) :: Fn
     Type(GridInfo), INTENT(OUT) :: Grid
     Type(OrbitInfo), INTENT(OUT) :: AEOrbit,FCOrbit
     Type(PotentialInfo), INTENT(OUT) :: AEPot,FCPot
     Type(SCFInfo), INTENT(OUT) :: AESCF,FCSCF
     Type(FCInfo), INTENT(OUT) :: FC

     INTEGER, parameter :: ifo=15
     INTEGER :: n,norbit,nps,npp,npd,npf,npg,i,io
     CHARACTER(132) :: inputfile

     Write(STD_OUT,*) 'Note that load has been called; code not completely checked'

     open(ifo,file=Fn,form='formatted',status='old')
     write(std_out,*) 'Reading dump file ', Fn

     read(ifo,*) inputfile(1:7)
     if (inputfile(1:7)=='ATOMPAW') then
       write(std_out,*) 'File seems to be fine'
     else
       write(std_out,*) 'File is not correct -- program will stop'
       stop
     endif

    ! grid info
     read(ifo,*) Grid%TYPE,Grid%n,Grid%h
     n=Grid%n
     Allocate(Grid%r(n))
     read(ifo,*)(Grid%r(i),i=1,n)
     if (usingloggrid(Grid)) then
        Allocate(Grid%drdu(n),Grid%pref(n),Grid%rr02(n))
          read(ifo,*)(Grid%drdu(i),Grid%pref(i),Grid%rr02(i),i=1,n)
     endif

    ! atomdata fixed constants
     read(ifo,*)  frozencorecalculation,setupfrozencore,scalarrelativistic ,&
&     finitenucleus,gaussianshapefunction,besselshapefunction,ColleSalvetti

    ! Orbit info
     read(ifo,*) exctype,nps,npp,npd,npf,npg,norbit
     Call InitOrbit(AEOrbit,norbit,n,exctype)
     AEOrbit%nps=nps;AEOrbit%npp=npp;AEOrbit%npd=npd;AEOrbit%npf=npf;AEOrbit%npg=npg
     read(ifo,*) (AEOrbit%np(io),AEOrbit%l(io),AEOrbit%iscore(io),&
&        AEOrbit%eig(io),AEOrbit%occ(io),io=1,AEOrbit%norbit)
     read(ifo,*) ((AEOrbit%wfn(i,io),i=1,n),io=1,AEOrbit%norbit)
     read(ifo,*) (AEOrbit%den(i),i=1,n)
     read(ifo,*) (AEOrbit%tau(i),i=1,n)

     read(ifo,*) exctype,nps,npp,npd,npf,npg,norbit
     Call InitOrbit(FCOrbit,norbit,n,exctype)
     FCOrbit%nps=nps;FCOrbit%npp=npp;FCOrbit%npd=npd;FCOrbit%npf=npf;FCOrbit%npg=npg
     read(ifo,*) (FCOrbit%np(io),FCOrbit%l(io),FCOrbit%iscore(io),&
&        FCOrbit%eig(io),FCOrbit%occ(io),io=1,FCOrbit%norbit)
     read(ifo,*) ((FCOrbit%wfn(i,io),i=1,n),io=1,FCOrbit%norbit)
     read(ifo,*) (FCOrbit%den(i),i=1,n)

     exctype=AEOrbit%exctype
     IF(TRIM(exctype)/='EXX'.and.TRIM(exctype)/='EXXCS'.and.TRIM(exctype)/='EXXKLI') CALL initexch

   ! Pot info
     Call InitPot(AEPot,n)
     read(ifo,*) AEPot%nz,AEPot%sym,AEPot%q,AEPot%v0,AEPot%v0p
     read(ifo,*) (AEPot%rv(i),AEPot%rvn(i),AEPot%rvh(i),AEPot%rvx(i),i=1,n)

     Call InitPot(FCPot,n)
     read(ifo,*) FCPot%nz,FCPot%sym,FCPot%q,FCPot%v0,FCPot%v0p
     read(ifo,*) (FCPot%rv(i),FCPot%rvn(i),FCPot%rvh(i),FCPot%rvx(i),i=1,n)

   ! SCFinfo
     Call InitSCF(AESCF)
     read(ifo,*) AESCF%iter,AESCF%delta,AESCF%eone,AESCF%ekin,&
&       AESCF%estatic,AESCF%ecoul,AESCF%eexc,AESCF%oepcs,AESCF%etot,&
&       AESCF%valekin,AESCF%valecoul,AESCF%valeexc,AESCF%corekin,AESCF%evale

     Call InitSCF(FCSCF)
     read(ifo,*) FCSCF%iter,FCSCF%delta,FCSCF%eone,FCSCF%ekin,&
&       FCSCF%estatic,FCSCF%ecoul,FCSCF%eexc,FCSCF%oepcs,FCSCF%etot,&
&       FCSCF%valekin,FCSCF%valecoul,FCSCF%valeexc,FCSCF%corekin,FCSCF%evale

   ! FCinfo
     call InitFC(FC,n)
      read(ifo,*) FC%zvale,FC%zcore
      read(ifo,*) (FC%coreden(i),FC%valeden(i),i=1,n)
      read(ifo,*) (FC%coretau(i),FC%valetau(i),i=1,n)

   If (TRIM(AEOrbit%exctype)=='EXX'.or.TRIM(AEOrbit%exctype)=='EXXKLI') &
&                   Call EXXload(Grid,AEOrbit,ifo)
   If (TRIM(AEOrbit%exctype)=='HF'.or.TRIM(AEOrbit%exctype)=='HFV') &
&           Call HFload(Grid,AEOrbit,ifo)

   close(ifo)

   write(std_out,*) 'Closing load file'
   write(std_out,*) 'Load completed normally'

  END SUBROUTINE load_aeatom

END  MODULE AEatom
