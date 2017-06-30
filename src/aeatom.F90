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

  USE Tools
  USE atomdata
  USE excor
  USE exx_mod
  USE general_mod
  USE globalmath
  USE gridmod
  USE hf_mod
  USE ldagga_mod

  IMPLICIT NONE

  REAL(8), PARAMETER, PRIVATE ::small0=1.d-13,rimix=0.5d0,worst=1.d-8
  REAL(8), PARAMETER, PRIVATE ::linrange=50.d0,linh=0.0025d0,mxgridlin=20001
  REAL(8), PARAMETER, PRIVATE ::logrange=80.d0,logh=0.020d0,mxgridlog=2001
!! slightly different grid introduced in version 4 !!
  REAL(8), PARAMETER, PRIVATE :: v4logrange=100.d0,lor00=1.d-5
!!  
  REAL(8), PARAMETER, PRIVATE ::logder_min=-5.d0,logder_max=4.95d0,logder_pts=200
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

  SUBROUTINE SCFatom_Init(ifinput)
    INTEGER, INTENT(IN),OPTIONAL :: ifinput

    !  program to calculate the self-consistent density functional
    !    atom ground state for atom with atomic number nz
    !    for self-consistent potential rv

    REAL(8) :: xocc,qf,small,zeff,hval,gridmatch,gridrange,logdmin,logdmax
    REAL(8) :: qcal,rescale,nzeff,h,hadjusted,r0
    INTEGER :: icount,i,j,k,it,start,np,ierr,gridpoints,logdpoints
    INTEGER :: is,ip,id,jf,ig,io,l,nfix,ir
    INTEGER :: ilin,ilog,inrl,iscl,ipnt,ifin,iend,ilgd,ihfpp,ilcex
    INTEGER :: igrid,irelat,ilogder,ilogv4,ibd
    INTEGER :: nps,npp,npd,npf,npg,finitenucleusmodel
    INTEGER, ALLOCATABLE :: nl(:,:)
    CHARACTER(128) :: exchangecorrelationandgridline,gridkey,relkey
    CHARACTER(132) :: inputline,inputword
    CHARACTER(1) :: CHR

    BDsolve=.false.
    scalarrelativistic=.FALSE.; finitenucleus=.FALSE. ;
    frozencorecalculation=.FALSE.;setupfrozencore=.false.
    gaussianshapefunction=.FALSE.;besselshapefunction=.FALSE.
    ColleSalvetti=.FALSE.  ; HFpostprocess=.FALSE.
    localizedcoreexchange=.false.
    hadjusted=0.d0 ; finitenucleusmodel=-1
    AEPot%sym="";AEPot%nz=0;AEPot%zz=0.d0;AEPot%q=0.d0;
    AEPot%v0=0.d0;AEPot%v0p=0.d0;AEPot%Nv0=0;AEPot%Nv0p=0 

!   First line : Atomic symbol and atomic number
!   ---------------------------------------------
    WRITE(6,*) 'enter atomic symbol and atomic number'
    IF (PRESENT(ifinput)) THEN
       READ(5,'(a)') inputline
       WRITE(ifinput,'(a)') TRIM(inputline)
       READ(inputline,*) AEPot%sym,AEPot%nz
    ELSE
       READ(5,*) AEPot%sym,AEPot%nz
    ENDIF
    AEPot%zz=AEPot%nz

!   2nd line : XC type, grid data, relativistic, point-nucleus, logderiv data, HF data, BDsolve keyword 
!   ----------------------------------------------------------------------------------
    WRITE(6,*) 'exchange-correlation type, among the following:'
    WRITE(6,*) '    * LDA-PW (default), GGA-PBE, GGA-PBESOL'
    WRITE(6,*) '    * LibXC keyword beginning with XC_'
    WRITE(6,*) '    * EXX, EXXOCC, EXXKLI, EXXCS'
    WRITE(6,*) '    * HF, HFV'
    WRITE(6,*) 'further optionally (space) "nonrelativistic/scalarrelativistic" keyword'
    WRITE(6,*) 'further optionally (space) "point-nucleus/finite-nucleus" keyword'
    WRITE(6,*) 'optionally (space) "loggrid/lineargrid" keyword if appropriate'
    WRITE(6,*) 'Note: "loggridv4" puts more points near origin'
    WRITE(6,*) '    further optionally n (number of grid points)'
    WRITE(6,*) '                       r_max (max. grid radius)'
    WRITE(6,*) '                       r_match (exact value of r(n))'
    WRITE(6,*) 'optionally (space) "logderivrange" keyword'
    WRITE(6,*) '    further optionally emin (minimum energy for log. deriv. plot)'
    WRITE(6,*) '                       emax (maximum energy for log. deriv. plot)'
    WRITE(6,*) '                       ne   (#  of energies for log. deriv. plot)'
    WRITE(6,*) 'addition option for "BDsolve" keyword for Block-Davidson solver'

!   Treat line characters
    READ(5,'(a)') exchangecorrelationandgridline
    if(present(ifinput)) WRITE(ifinput,'(a)') TRIM(exchangecorrelationandgridline)
    call Uppercase(exchangecorrelationandgridline)
    exchangecorrelationandgridline=trim(exchangecorrelationandgridline)
    write(6,*) 'writing exchangecorrelationandgridline: ',exchangecorrelationandgridline

!   Retrieve keyword indexes
    ilin=0;ilin=0;ilog=0;ilogv4=0;inrl=0;iscl=0;ipnt=0;ifin=0;ilgd=0
    ihfpp=0;ilcex=0;igrid=0;irelat=0;ilogder=0;ibd=0
    ilin=INDEX(exchangecorrelationandgridline,'LINEARGRID')
    ilog=INDEX(exchangecorrelationandgridline,'LOGGRID')
    ilogv4=INDEX(exchangecorrelationandgridline,'LOGGRIDV4')
    ibd=INDEX(exchangecorrelationandgridline,'BDSOLVE')
    inrl=INDEX(exchangecorrelationandgridline,'NONRELATIVISTIC')
    iscl=INDEX(exchangecorrelationandgridline,'SCALARRELATIVISTIC')
    ipnt=INDEX(exchangecorrelationandgridline,'POINT-NUCLEUS')
    ifin=INDEX(exchangecorrelationandgridline,'FINITE-NUCLEUS')
    ilgd=INDEX(exchangecorrelationandgridline,'LOGDERIVRANGE')
    ihfpp=INDEX(exchangecorrelationandgridline,'HFPOSTPROCESS')
    ilcex=INDEX(exchangecorrelationandgridline,'LOCALIZEDCOREEXCHANGE')
    igrid=max(ilin,ilog)
    irelat=max(inrl,iscl)
    ilogder=ilgd
    if (ifin>0) then
      READ(exchangecorrelationandgridline(ifin+14:ifin+14),'(a)') CHR
      write(6,*) 'CHR ', CHR
      if (CHR==" ") AEPot%finitenucleusmodel=0      
      if (CHR=="2") AEPot%finitenucleusmodel=2      
      if (CHR=="3") AEPot%finitenucleusmodel=3      
      if (CHR=="4") AEPot%finitenucleusmodel=4      
      if (CHR=="5") AEPot%finitenucleusmodel=5      
      write(6,*) 'Finite nucleus model number ', AEPot%finitenucleusmodel
    endif  

!   Treat simple logical variables
    if (iscl>0.and.inrl==0) scalarrelativistic=.true.
    if (ifin>0.and.ipnt==0) finitenucleus=.true.
    if (ihfpp>0) HFpostprocess=.true.
    if (ilcex>0) localizedcoreexchange=.true.
    write(6,'(" Logical variables scalarrelativistic,finitenucleus,HFpostprocess - ", 3L3)') &
&      scalarrelativistic,finitenucleus,HFpostprocess

!   Treat grid data
    gridkey='LINEAR';gridpoints=mxgridlin
    gridrange=linrange;gridmatch=linrange;hadjusted=linh
    if (ilog>0.and.ilin==0.and.ilogv4==0) then
     gridkey='LOGGRID';gridpoints=mxgridlog;
     gridrange=logrange;gridmatch=logrange
    endif
    if (ilog>0.and.ilin==0.and.ilogv4>0) then
     gridkey='LOGGRID';gridpoints=mxgridlog;
     gridrange=v4logrange;gridmatch=v4logrange
    endif
    if (igrid>0) then
     iend=128
     if (irelat >igrid.and.irelat-1 <iend) iend=irelat -1
     if (ilogder>igrid.and.ilogder-1<iend) iend=ilogder-1
     if (ibd>igrid.and.ibd-1<iend) iend=ibd-1
     inputline=""
     if (ilog>0.and.ilogv4==0.and.iend>igrid+7) &
&       inputline=trim(exchangecorrelationandgridline(igrid+7:iend))
     if (ilog>0.and.ilogv4>0.and.iend>igrid+9) &
&       inputline=trim(exchangecorrelationandgridline(igrid+9:iend))
     if (ilin>0.and.iend>igrid+10) &
&       inputline=trim(exchangecorrelationandgridline(igrid+10:iend))
     if (inputline/="") then
      call extractword(1,inputline,inputword);inputword=trim(inputword)
      if (inputword/="") then
       read(inputword,*) gridpoints
       call extractword(2,inputline,inputword);inputword=trim(inputword)
       if (inputword/="") then
        read(inputword,*) gridrange
        gridmatch=gridrange
        call extractword(3,inputline,inputword);inputword=trim(inputword)
        if (inputword/="") read(inputword,*) gridmatch
       endif
      endif
     endif
     if (gridpoints<=0) stop "Number of grid points should be >0 !"
    endif
!   Initialize grid data
    IF (TRIM(gridkey)=='LOGGRID'.and.ilogv4>0) THEN
       hval=logh
       CALL findh_given_r0(AEPot%zz,gridmatch,lor00,gridpoints,hval)
       CALL InitGrid(Grid,hval,gridrange,lor00/AEPot%zz)
    ELSEIF (TRIM(gridkey)=='LOGGRID'.and.ilogv4==0) THEN
       hval=logh
       CALL findh(AEPot%zz,gridmatch,gridpoints,hval,r0)
       CALL InitGrid(Grid,hval,gridrange,r0)
    ELSE
       hadjusted=gridmatch/(gridpoints-1)
       CALL InitGrid(Grid,hadjusted,gridrange)
    ENDIF
     
    CALL InitPot(AEPot,Grid%n)
    CALL Get_Nuclearpotential(Grid,AEPot)

!   Treat logderiv data
    minlogderiv=logder_min;maxlogderiv=logder_max;nlogderiv=logder_pts
    if (ilogder>0) then
     iend=128
     if (igrid >ilogder.and.igrid-1 <iend) iend=igrid -1
     if (irelat>ilogder.and.irelat-1<iend) iend=irelat-1
     inputline=""
     if (iend>ilogder+13) inputline=trim(exchangecorrelationandgridline(ilogder+13:iend))
     if (inputline/="") then
      call extractword(1,inputline,inputword);inputword=trim(inputword)
      if (inputword/="") then
       read(inputword,*) minlogderiv
       call extractword(2,inputline,inputword);inputword=trim(inputword)
       if (inputword/="") then
        read(inputword,*) maxlogderiv
        call extractword(3,inputline,inputword);inputword=trim(inputword)
        if (inputword/="") read(inputword,*) nlogderiv
       endif
      endif
     endif
    endif

! consider bound state solver
   if (ibd>0) BDsolve=.true.
   if (BDsolve.and.gridkey=='LINEAR') then
      write(6,*) &
&      'WARNING:   BlockDavidson solver works very slowly with linear grid'
   endif
   Write(6,*) 'BDSOLVE', BDSOLVE;call flush_unit(6)

!   Treat exchange-correlation/HF keyword
    READ(unit=exchangecorrelationandgridline,fmt=*) exctype
    IF (TRIM(exctype)=='EXX'.or.TRIM(exctype)=='EXXKLI'.or.   &
&       TRIM(exctype)=='EXXOCC') THEN
       CALL EXX_Input_Settings(exchangecorrelationandgridline)
    ELSE IF (TRIM(exctype)=='EXXCS') THEN
       CALL EXX_Input_Settings(exchangecorrelationandgridline)
       ColleSalvetti=.TRUE.
    ELSE IF (TRIM(exctype)=='HF'.or.TRIM(exctype)=='HFV') THEN
    ELSE
       CALL initexch
    ENDIF

!   3rd line and following : electronic configuration of atom
!   ----------------------------------------------------------------------------------
    WRITE(6,'(a,f6.2)') ' Calculation for atomic number = ',AEPot%zz
    call flush_unit(6)
    WRITE(6,*) 'enter maximum principal quantum numbers for s,p,d,f,g'
    IF(PRESENT(ifinput)) THEN
       READ(5,'(a)') inputline
       WRITE(ifinput,'(a)') TRIM(inputline)
       READ(inputline,*) nps,npp,npd,npf,npg
    ELSE
       READ(5,*) nps,npp,npd,npf,npg
    ENDIF
    IF(nps<0)nps=0
    IF(npp<0)npp=0
    IF(npd<0)npd=0
    IF(npf<0)npf=0
    IF(npg<0)npg=0
    WRITE(6,'(5i4)') nps,npp,npd,npf,npg
    i=MAX(nps,npp,npd,npf,npg)
    j=nps
    IF(npp>0) j=j+npp-1
    IF(npd>0) j=j+npd-2
    IF(npf>0) j=j+npf-3
    IF(npg>0) j=j+npg-4

    CALL InitOrbit(AEOrbit,j,Grid%n,exctype)
    AEOrbit%nps=nps;AEOrbit%npp=npp;AEOrbit%npd=npd
    AEOrbit%npf=npf;AEOrbit%npg=npg
    ALLOCATE(nl(i,j));nl=0

    icount=0
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
    WRITE(6,*) AEOrbit%norbit, ' orbitals will be calculated'
    !
    WRITE(6,*)' Below are listed the default occupations '
    WRITE(6,"(' n  l     occupancy')")
    DO io=1,AEOrbit%norbit
       WRITE(6,'(i2,1x,i2,4x,1p,1e15.7)') &
&           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
    ENDDO
    !
    WRITE(6,*)' enter np l occ for all occupations for all revisions'
    WRITE(6,*)'  enter 0 0 0. to end'

    DO
       IF(PRESENT(ifinput)) THEN
          READ(5,'(a)') inputline
          WRITE(ifinput,'(a)') TRIM(inputline)
          READ(inputline,*) ip,l,xocc
       ELSE
          READ(5,*) ip,l,xocc
       ENDIF
       IF (ip.LE.0) EXIT
       nfix=nl(ip,l+1)
       IF (nfix.LE.0.OR.nfix.GT.AEOrbit%norbit) THEN
          WRITE(6,*) 'error in occupations -- ip,l,xocc', &
&              ip,l,xocc,nfix,AEOrbit%norbit
          STOP
       ENDIF
       AEOrbit%occ(nfix)=xocc
    ENDDO

    !
    WRITE(6,*) ' Corrected occupations are: '
    WRITE(6,"(' n  l     occupancy')")
    electrons=0.d0
    DO io=1,AEOrbit%norbit
       WRITE(6,'(i2,1x,i2,4x,1p,1e15.7)')  &
&           AEOrbit%np(io),AEOrbit%l(io),AEOrbit%occ(io)
       electrons=electrons+AEOrbit%occ(io)
    ENDDO
    AEPot%q=electrons
    qf=AEPot%nz-electrons
    WRITE(6,*)
    WRITE(6,*) 'nuclear charge    = ', AEPot%nz
    WRITE(6,*) 'electronic charge = ', electrons
    WRITE(6,*) 'net charge        = ', qf

    CALL InitSCF(AESCF)
    IF (scalarrelativistic) CALL Allocate_Scalar_Relativistic(Grid)

    write(6,*) 'Finish SCFatom_Init' ; call flush_unit(6)

    DEALLOCATE(nl)

  END SUBROUTINE SCFatom_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  SCFatom                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SCFatom(scftype,lotsofoutput,ifinput)
    IMPLICIT NONE
    CHARACTER(2),INTENT(IN)::scftype
    LOGICAL, INTENT(IN) :: lotsofoutput
    INTEGER, INTENT(IN), OPTIONAL :: ifinput
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
       CALL NC_Init(OrbitPtr,PotPtr)
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
       CALL SC_Init(ifinput,OrbitPtr,PotPtr,SCFPtr)

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
       write(6,*) 'Just before HF_SCF'; call flush_unit(6)
       CALL HF_SCF(scftype,lotsofoutput,Grid,OrbitPtr,PotPtr,FC,SCFPtr)
       write(6,*) 'Just after HF_SCF'; call flush_unit(6)

    ELSE
      !LDA or GGA
      CALL LDAGGA_SCF(scftype,lotsofoutput,Grid,OrbitPtr,PotPtr,FC,SCFPtr)
    ENDIF

    If (HFpostprocess) then
       write(6,*) "******PostProcessing HF*********"
       CALL InitSCF(SCFPP)
       call hf_energy_only(Grid,OrbitPtr,PotPtr,SCFPP)
    endif
  END SUBROUTINE SCFatom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Orbit_Init
  !!       From nuclear charge -- generate hydrogenic-like initial wfns
  !!          and densities --
  !!          fill AEOrbit%wfn, AEOrbit%eig, and AEOrbit%den and AEOrbit%q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Orbit_Init(Orbit,Pot)
    IMPLICIT NONE
    TYPE(PotentialInfo),INTENT(INOUT) :: Pot
    TYPE(OrbitInfo),INTENT(INOUT) :: Orbit
    INTEGER  :: i,io,ir,ip,l,nfix,np
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,xocc
    INTEGER :: initialconfig=0

       IF (initialconfig/=0) STOP 'Error in aeatom -- Orbit_Init already called'

       !  calculate initial charge density from hydrogen-like functions
       !  also initial energies
       zeff=Pot%nz
       DO io=1,Orbit%norbit
          np=Orbit%np(io)
          l=Orbit%l(io)
          xocc=Orbit%occ(io)
          Orbit%eig(io)=-(zeff/(np))**2
          WRITE(6,*) io,np,l,xocc,Orbit%eig(io)
          DO ir=1,Grid%n
             Orbit%wfn(ir,io)=hwfn(zeff,np,l,Grid%r(ir))
             IF (ABS(Orbit%wfn(ir,io))<machine_zero) Orbit%wfn(ir,io)=0.d0
          ENDDO
          zeff=zeff-0.5d0*xocc
          zeff=MAX(zeff,1.d0)
       ENDDO

       ! check charge and rescale
       Orbit%den=0.d0
       DO io=1,Orbit%norbit
          xocc=Orbit%occ(io)
          DO ir=1,Grid%n
             Orbit%den(ir)=Orbit%den(ir)+xocc*(Orbit%wfn(ir,io)**2)
          ENDDO
       ENDDO
       qcal=integrator(Grid,Orbit%den)
       qf=qcal
       !WRITE(6,*) 'qcal electrons = ',qcal, electrons
       !rescale density
       rescale=electrons/qcal
       Orbit%den(1:Grid%n)=Orbit%den(1:Grid%n)*rescale

       initialconfig=1
        write(6,*) 'completed Orbit_Init '; call flush_unit(6)

  END SUBROUTINE Orbit_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE NC_Init(Orbit,Pot)
    TYPE(PotentialInfo),INTENT(INOUT) :: Pot
    TYPE(OrbitInfo),INTENT(INOUT) :: Orbit
    INTEGER  :: i,io,ir,ip,l,nfix,j,np
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,xocc
       !----------New Configuration--AE calculation--------------
       ! readin revision and normalize the density

       WRITE(6,*)' enter np l occ for all occupations for all revisions'
       WRITE(6,*)'  enter 0 0 0. to end'

       DO
          READ(5,*) ip,l,xocc
          IF (ip.LE.0) EXIT
          nfix=-100
          DO io=1,Orbit%norbit
             IF (ip==Orbit%np(io).AND.l==Orbit%l(io)) THEN
                nfix=io
                EXIT
             ENDIF
          ENDDO
          IF (nfix.LE.0.OR.nfix.GT.Orbit%norbit) THEN
             WRITE(6,*) 'error in occupations -- ip,l,xocc',ip,l,xocc,nfix,Orbit%norbit
             STOP
          ENDIF
          Orbit%occ(nfix)=xocc
       ENDDO

       WRITE(6,*) ' Corrected occupations are: '
       WRITE(6,"(' n  l     occupancy')")
       electrons=0.d0
       DO io=1,Orbit%norbit
          WRITE(6,'(i2,1x,i2,4x,1p,1e15.7)') Orbit%np(io),Orbit%l(io),Orbit%occ(io)
          electrons=electrons+Orbit%occ(io)
       ENDDO
       Pot%q=electrons
       qf=Pot%nz-electrons
       WRITE(6,*)
       WRITE(6,*) 'nuclear charge    = ' , Pot%nz
       WRITE(6,*) 'electronic charge = ', electrons
       WRITE(6,*) 'net charge        = ', qf
       !
       !
       !  calculate initial charge density from stored wavefunctions
       !    also initial energies
       !
       Orbit%den=0.d0
       DO io=1,Orbit%norbit
          xocc=Orbit%occ(io)
          DO ir=1,Grid%n
             Orbit%den(ir)=Orbit%den(ir)+xocc*(Orbit%wfn(ir,io)**2)
          ENDDO
       ENDDO
       !
       !  check charge
       !
       qcal=integrator(Grid,Orbit%den)
       qf=qcal
       !WRITE(6,*) 'qcal electrons = ',qcal, electrons
       !  rescale density
       rescale=electrons/qcal
       Orbit%den=Orbit%den*rescale

  END SUBROUTINE NC_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!         SC_Init -- set core states in current configuration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SC_Init(ifinput,Orbit,Pot,SCF)
    TYPE(PotentialInfo),INTENT(INOUT) :: Pot
    TYPE(OrbitInfo),INTENT(INOUT) :: Orbit
    TYPE (SCFInfo),INTENT(INOUT) :: SCF
    INTEGER :: choosevalence=0
    INTEGER, INTENT(IN), OPTIONAL :: ifinput

    INTEGER :: io
    LOGICAL :: noalt=.true.
    If (choosevalence/=0) then
     write(6,*) 'Error in aeatom -- SC_Init already called'
       stop
    endif

    CALL Set_Valence(ifinput)
    Call Core_Electron_Report(Orbit,FC,6)
    Call Valence_Electron_Report(Orbit,FC,6)

    IF (TRIM(Orbit%exctype)=='EXX') THEN
       Call Get_FCKinCoul(Grid,Pot,Orbit,FC,SCF)
       CALL Get_FCEnergy_EXX(Grid,Orbit,FC,SCF)
       CALL Set_Vxref(AEPot)   ! probably not a good idea
    ELSEIF (TRIM(Orbit%exctype)=='EXXKLI') THEN
       Call Get_FCKinCoul(Grid,Pot,Orbit,FC,SCF,noalt)
       CALL Get_FCEnergy_EXX(Grid,Orbit,FC,SCF)
    ELSEIF (TRIM(Orbit%exctype)=='HF'.or.TRIM(Orbit%exctype)=='HFV') THEN
       write(6,*) ' Completed Setcore '; call flush_unit(6)
       CALL Get_FCEnergy_HF(Grid,Pot,Orbit,FC,SCF)
       CALL Total_FCEnergy_Report(SCF,6)
       write(6,*) ' Completed Setcore '; call flush_unit(6)
    ELSE
       Call Get_FCKinCoul(Grid,Pot,Orbit,FC,SCF)
       CALL Get_FCEXC(SCF)
       CALL Total_FCEnergy_Report(SCF,6)
    ENDIF
    !CALL Total_Energy_Report(SCF,6)
    !   For EXX, exchange contribution is not known yet.

    choosevalence=1

  END SUBROUTINE SC_Init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE FC_Init(Orbit,Pot)
    TYPE(PotentialInfo),INTENT(INOUT) :: Pot
    TYPE(OrbitInfo),INTENT(INOUT) :: Orbit
    INTEGER  :: i,io,ir,ip,l,nfix,j,np
    REAL(8) :: en0,qcal,qf,rescale,z,small,zeff,xocc

       WRITE(6,*)' Below are listed the current valence occupations '
       WRITE(6,"(' n  l     occupancy')")
       DO io=1,Orbit%norbit
          IF (.NOT.Orbit%iscore(io)) WRITE(6,'(i2,1x,i2,4x,1p,1e15.7)') &
&              Orbit%np(io),Orbit%l(io),Orbit%occ(io)
       ENDDO
       !
       WRITE(6,*)' enter np l occ for all occupations for all revisions'
       WRITE(6,*)'  enter 0 0 0. to end'

       DO
          READ(5,*) ip,l,xocc
          IF (ip.LE.0) EXIT
          nfix=-100
          DO io=1,Orbit%norbit
             IF (ip==Orbit%np(io).AND.l==Orbit%l(io) &
&                 .AND.(.NOT.Orbit%iscore(io))) THEN
                nfix=io
                EXIT
             ENDIF
          ENDDO
          IF (nfix.LE.0.OR.nfix.GT.Orbit%norbit) THEN
             WRITE(6,*) 'error in occupations -- ip,l,xocc',ip,l,xocc,nfix,Orbit%norbit
             STOP
          ENDIF
          Orbit%occ(nfix)=xocc
       ENDDO

       !
       WRITE(6,*) ' Corrected occupations are: '
       WRITE(6,"(' n  l     occupancy')")

       FC%zvale=0.d0
       DO io=1,Orbit%norbit
          If (.not.Orbit%iscore(io))then
            WRITE(6,'(i2,1x,i2,4x,1p,1e15.7)')  &
&             Orbit%np(io),Orbit%l(io),Orbit%occ(io)
            FC%zvale=FC%zvale+Orbit%occ(io)
          Endif
       ENDDO

       electrons=FC%zvale+FC%zcore

       !
       !  calculate initial charge density from stored wavefunctions
       !    also initial energies
       !
       FC%valeden=0
       DO io=1,Orbit%norbit
          IF (.NOT.Orbit%iscore(io)) THEN
             xocc=Orbit%occ(io)
             DO ir=1,Grid%n
                FC%valeden(ir)=FC%valeden(ir)+xocc*(Orbit%wfn(ir,io)**2)
             ENDDO
          ENDIF
       ENDDO
       !
       !  check charge				!  rescale density

       !
       qcal=integrator(Grid,FC%valeden)
       if (ABS(qcal)<small0.OR.electrons<small0) then
          FC%valeden=0
       else
          rescale=FC%zvale/qcal
          FC%valeden=FC%valeden*rescale
       endif

       Orbit%den=FC%valeden+FC%coreden
       WRITE(6,*)
       WRITE(6,*) 'core charge    = ' , FC%zcore
       WRITE(6,*) 'valence charge = ',  FC%zvale
       WRITE(6,*) 'net charge     = ',  Pot%nz-FC%zcore-FC%zvale

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


    CALL poisson(Grid,Pot%q,Orbit%den,Pot%rvh,ecoul,v0)
    write(6,*) 'In Potential_Init', Pot%q,ecoul; call flush_unit(6)

  END SUBROUTINE Potential_Init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Set Valence Orbitals
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE Set_Valence(ifinput)
    INTEGER, INTENT(IN), OPTIONAL :: ifinput

    CHARACTER(1) :: answer
    INTEGER :: io,n,ok
    CHARACTER(132) :: inputline

    call InitFC(FC,Grid%n)

    WRITE(6,*) 'for each state enter c for core or v for valence'
    FC%zcore=0; FC%zvale=0
    FC%coreden=0; FC%valeden=0
    DO io=1,FCOrbit%norbit
       WRITE(6,'(3i5,1p,2e15.7)') io,FCOrbit%np(io),FCOrbit%l(io),&
&           FCOrbit%occ(io),FCOrbit%eig(io)
       DO
          IF (PRESENT(ifinput)) THEN
             READ(5,'(a)') inputline
             WRITE(ifinput,'(a)') TRIM(inputline)
             READ(inputline,*) answer
          ELSE
             READ(5,*) answer
          ENDIF
          IF (answer.NE.'c'.AND.answer.NE.'C'.AND.answer.NE.'v'&
&              .AND.answer.NE.'V')  THEN
             WRITE(6,*) 'Please input c or v'
          ELSE
             EXIT
          ENDIF
       ENDDO
       IF (answer.EQ.'c'.OR.answer.EQ.'C') FCOrbit%iscore(io)=.TRUE.
       IF (answer.EQ.'v'.OR.answer.EQ.'V') FCOrbit%iscore(io)=.FALSE.
       IF (FCOrbit%iscore(io)) THEN
          FC%zcore=FC%zcore+FCOrbit%occ(io)
          FC%coreden=FC%coreden+FCOrbit%occ(io)*(FCOrbit%wfn(:,io))**2
       ENDIF
       IF (.NOT.FCOrbit%iscore(io)) THEN

          FC%zvale=FC%zvale+FCOrbit%occ(io)
          FC%valeden=FC%valeden+FCOrbit%occ(io)*(FCOrbit%wfn(:,io))**2
       ENDIF
    ENDDO

         write(6,*) 'Returning from SetValence'
         write(6,*) 'Core electrons', FC%zcore,integrator(Grid,FC%coreden)
         write(6,*) 'Vale electrons', FC%zvale,integrator(Grid,FC%valeden)

  END SUBROUTINE Set_Valence

  SUBROUTINE dump_aeatom(Fn,Grid,AEOrbit,AEPot,AESCF,FCOrbit,FCPot,FCSCF,FC)
     Character(*), INTENT(IN) :: Fn
     Type(GridInfo), INTENT(IN) :: Grid
     Type(OrbitInfo), INTENT(IN) :: AEOrbit,FCOrbit
     Type(PotentialInfo), INTENT(IN) :: AEPot,FCPot
     Type(SCFInfo), INTENT(IN) :: AESCF,FCSCF
     Type(FCInfo), INTENT(IN) :: FC

     INTEGER, parameter :: ifo=15
     INTEGER :: i,io,n

     Write(6,*) 'Note that dump has been called; code not completely checked'

     open(ifo,file=Fn,form='formatted',status='replace')
     write(6,*) 'Creating or replacing dump file', Fn

     write(ifo,*) 'ATOMPAW'    ! keyword

    ! grid info
     write(ifo,*) Grid%TYPE,Grid%n,Grid%h
     n=Grid%n
     write(ifo,*)(Grid%r(i),i=1,n)
     if (usingloggrid(Grid)) then
          write(ifo,*)(Grid%drdu(i),Grid%pref(i),Grid%rr02(i),i=1,n)
     endif

    ! atomdata fixed constants
     write(ifo,*)  frozencorecalculation,setupfrozencore,scalarrelativistic ,&
&     finitenucleus,gaussianshapefunction,besselshapefunction,ColleSalvetti

    ! Orbit info
     write(ifo,*) AEOrbit%exctype,AEOrbit%nps, AEOrbit%npp, AEOrbit%npd ,&
&      AEOrbit%npf, AEOrbit%npg, AEOrbit%norbit
     write(ifo,*) (AEOrbit%np(io),AEOrbit%l(io),AEOrbit%iscore(io),&
&        AEOrbit%eig(io),AEOrbit%occ(io),io=1,AEOrbit%norbit)
     write(ifo,*) ((AEOrbit%wfn(i,io),i=1,n),io=1,AEOrbit%norbit)
     write(ifo,*) (AEOrbit%den(i),i=1,n)

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

   If (TRIM(AEOrbit%exctype)=='EXX'.or.TRIM(AEOrbit%exctype)=='EXXKLI') &
&                  Call EXXdump(Grid,AEOrbit,ifo)
   If (TRIM(AEOrbit%exctype)=='HF'.or.TRIM(AEOrbit%exctype)=='HFV') &
&          Call HFdump(Grid,AEOrbit,ifo)

   close(ifo)

   write(6,*) 'Closing dump file'
  END SUBROUTINE dump_aeatom

  SUBROUTINE load_aeatom(Fn,Grid,AEOrbit,AEPot,AESCF,FCOrbit,FCPot,FCSCF,FC,ifd)
     Character(*), INTENT(IN) :: Fn
     Type(GridInfo), INTENT(OUT) :: Grid
     Type(OrbitInfo), INTENT(OUT) :: AEOrbit,FCOrbit
     Type(PotentialInfo), INTENT(OUT) :: AEPot,FCPot
     Type(SCFInfo), INTENT(OUT) :: AESCF,FCSCF
     Type(FCInfo), INTENT(OUT) :: FC
     INTEGER, INTENT(IN) :: ifd

     INTEGER, parameter :: ifo=15
     INTEGER :: i,io,n,norbit,nps,npp,npd,npf,npg,j,k
     CHARACTER(132) :: inputfile,checkfile,exctype

     Write(6,*) 'Note that load has been called; code not completely checked'

     open(ifo,file=Fn,form='formatted',status='old')
     write(6,*) 'Reading dump file ', Fn

     read(ifo,*) inputfile(1:7)
     if (inputfile(1:7)=='ATOMPAW') then
       write(6,*) 'File seems to be fine'
     else
       write(6,*) 'File is not correct -- program will stop'
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

   If (TRIM(AEOrbit%exctype)=='EXX'.or.TRIM(AEOrbit%exctype)=='EXXKLI') &
&                   Call EXXload(Grid,AEOrbit,ifo)
   If (TRIM(AEOrbit%exctype)=='HF'.or.TRIM(AEOrbit%exctype)=='HFV') &
&           Call HFload(Grid,AEOrbit,ifo)

   close(ifo)

   write(6,*) 'Closing load file'

  ! crude checking by comparing input file (unit 5) with dummy file (unit ifd)
   rewind(5); rewind(ifd);
    Do i=1,3
      read(5,'(a)') inputfile; read(ifd,'(a)') checkfile
      If (TRIM(inputfile)/=TRIM(checkfile)) then
         write(6,*) 'Inconsistency in files'
          write(6,'(a)') inputfile
          write(6,'(a)') checkfile
          stop
      endif
    enddo

    Bigloop: Do
      read(5,'(a)') inputfile; read(ifd,'(a)') checkfile
      read(inputfile,*) i,j,k
      if (i==0) exit Bigloop
      If (TRIM(inputfile)/=TRIM(checkfile)) then
         write(6,*) 'Inconsistency in files'
          write(6,'(a)') inputfile
          write(6,'(a)') checkfile
          stop
      endif
    enddo Bigloop

    Do io=1,FCOrbit%norbit
      read(5,'(a)') inputfile; read(ifd,'(a)') checkfile
      If (TRIM(inputfile)/=TRIM(checkfile)) then
         write(6,*) 'Inconsistency in files'
          write(6,'(a)') inputfile
          write(6,'(a)') checkfile
          stop
      endif
    enddo

    write(6,*) 'Load completed normally'
  END SUBROUTINE load_aeatom
END  MODULE AEatom
